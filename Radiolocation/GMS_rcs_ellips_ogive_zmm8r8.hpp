

#ifndef __GMS_RCS_ELLIPS_OGIVE_ZMM8R8_HPP__
#define __GMS_RCS_ELLIPS_OGIVE_ZMM8R8_HPP__ 130320230926


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

    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM8R8_FULLVER =
      1000U*GMS_RCS_ELLIPS_OGIVE_ZMM8R8_MAJOR+
      100U*GMS_RCS_ELLIPS_OGIVE_ZMM8R8_MINOR+
      10U*GMS_RCS_ELLIPS_OGIVE_ZMM8R8_MICRO;
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM8R8_CREATION_DATE = "13-03-2023 09:26 AM +00200 (MON 13 MAR 2023 GMT+2)";
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM8R8_DESCRIPTION   = "AVX512 optimized Ellipsoid and Ogive Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimddp.hpp"
#include "GMS_complex_zmm8r8.hpp"

namespace gms {


          namespace radiolocation {


                      /*
                            High-frequency cross-section of perfectly
                            conducting ellipsoid.
                            Bistatic case.
                            Formula 5.1.54
                        */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5154_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d c,
                                            const __m512d th1,
                                            const __m512d phi1,
                                            const __m512d th2,
                                            const __m512d phi2) {

                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,a2,b2,c2,sth1,cphi1,sth2,cphi2,cth1,cth2;
                          register __m512d trm1,trm2,trm3,num,den,x0,strm1,strm2,strm3;
                          a2   = _mm512_mul_pd(a,a);
                          sth1 = xsin(phi1);
                          b2   = _mm512_mul_pd(b,b);
                          cphi1= xcos(phi1);
                          c2   = _mm512_mul_pd(c,c);
                          sth2 = xsin(th2);
                          x0   = _mm512_mul_pd(a2,_mm512_mul_pd(b2,c2));
                          cphi2= xcos(phi2);
                          num  = _mm512_mul_pd(_4pi,x0);
                          cth1 = xcos(th1);
                          trm1 = _mm512_fmadd_pd(sth1,cphi1,_mm512_mul_pd(sth2,cphi2));
                          strm1= _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          cth2 = xcos(th2);
                          trm2 = _mm512_fmadd_pd(sth1,sphi1,_mm512_mul_pd(sth2,sphi2));
                          strm2= _mm512_mul_pd(b2,_mm512_mul_pd(trm2,trm2));
                          trm3 = _mm512_mul_pd(cth1,cth2);
                          strm3= _mm512_mul_pd(c2,_mm512_mul_pd(trm3,trm3));
                          x0   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5154_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pb,
                                            const double * __restrict __ATTR_ALIGN__(64) pc,
                                            const double * __restrict __ATTR_ALIGN__(64) pth1,
                                            const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                            const double * __restrict __ATTR_ALIGN__(64) pth2,
                                            const double * __restrict __ATTR_ALIGN__(64) pphi2) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d b   = _mm512_load_pd(&pb[0]);
                          register __m512d c   = _mm512_load_pd(&pc[0]);
                          register __m512d th1 = _mm512_load_pd(&pth1[0]);
                          register __m512d phi1= _mm512_load_pd(&pphi1[0]);
                          register __m512d th2 = _mm512_load_pd(&pth2[0]);
                          register __m512d phi2= _mm512_load_pd(&pphi2[0]);
                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,a2,b2,c2,sth1,cphi1,sth2,cphi2,cth1,cth2;
                          register __m512d trm1,trm2,trm3,num,den,x0,strm1,strm2,strm3;
                          a2   = _mm512_mul_pd(a,a);
                          sth1 = xsin(phi1);
                          b2   = _mm512_mul_pd(b,b);
                          cphi1= xcos(phi1);
                          c2   = _mm512_mul_pd(c,c);
                          sth2 = xsin(th2);
                          x0   = _mm512_mul_pd(a2,_mm512_mul_pd(b2,c2));
                          cphi2= xcos(phi2);
                          num  = _mm512_mul_pd(_4pi,x0);
                          cth1 = xcos(th1);
                          trm1 = _mm512_fmadd_pd(sth1,cphi1,_mm512_mul_pd(sth2,cphi2));
                          strm1= _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          cth2 = xcos(th2);
                          trm2 = _mm512_fmadd_pd(sth1,sphi1,_mm512_mul_pd(sth2,sphi2));
                          strm2= _mm512_mul_pd(b2,_mm512_mul_pd(trm2,trm2));
                          trm3 = _mm512_mul_pd(cth1,cth2);
                          strm3= _mm512_mul_pd(c2,_mm512_mul_pd(trm3,trm3));
                          x0   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5154_zmm8r8_u(const double * __restrict  pa,
                                            const double * __restrict  pb,
                                            const double * __restrict  pc,
                                            const double * __restrict  pth1,
                                            const double * __restrict  pphi1,
                                            const double * __restrict  pth2,
                                            const double * __restrict  pphi2) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d b   = _mm512_loadu_pd(&pb[0]);
                          register __m512d c   = _mm512_loadu_pd(&pc[0]);
                          register __m512d th1 = _mm512_loadu_pd(&pth1[0]);
                          register __m512d phi1= _mm512_loadu_pd(&pphi1[0]);
                          register __m512d th2 = _mm512_loadu_pd(&pth2[0]);
                          register __m512d phi2= _mm512_loadu_pd(&pphi2[0]);
                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,a2,b2,c2,sth1,cphi1,sth2,cphi2,cth1,cth2;
                          register __m512d trm1,trm2,trm3,num,den,x0,strm1,strm2,strm3;
                          a2   = _mm512_mul_pd(a,a);
                          sth1 = xsin(phi1);
                          b2   = _mm512_mul_pd(b,b);
                          cphi1= xcos(phi1);
                          c2   = _mm512_mul_pd(c,c);
                          sth2 = xsin(th2);
                          x0   = _mm512_mul_pd(a2,_mm512_mul_pd(b2,c2));
                          cphi2= xcos(phi2);
                          num  = _mm512_mul_pd(_4pi,x0);
                          cth1 = xcos(th1);
                          trm1 = _mm512_fmadd_pd(sth1,cphi1,_mm512_mul_pd(sth2,cphi2));
                          strm1= _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          cth2 = xcos(th2);
                          trm2 = _mm512_fmadd_pd(sth1,sphi1,_mm512_mul_pd(sth2,sphi2));
                          strm2= _mm512_mul_pd(b2,_mm512_mul_pd(trm2,trm2));
                          trm3 = _mm512_mul_pd(cth1,cth2);
                          strm3= _mm512_mul_pd(c2,_mm512_mul_pd(trm3,trm3));
                          x0   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                 }


                       /*
                            High-frequency cross-section of perfectly
                            conducting ellipsoid.
                            Backscattering case.
                            Formula 5.1.55
                        */ 


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5155_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d c,
                                            const __m512d th,
                                            const __m512d phi) {

                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a2,b2,c2,sphi,sphis,cth,cths,sth,sths,cphi,cphis;
                          register __m512d num,x0,trm1;
                          a2  = _mm512_mul_pd(a,a);
                          sphi= xsin(phi);
                          b2  = _mm512_mul_pd(b,b);
                          sphis= _mm512_mul_pd(sphi,sphi);
                          cth = xcos(th);
                          c2  = _mm512_mul_pd(c,c);
                          cths= _mm512_mul_pd(cth,cth);
                          sth = xsin(th);
                          num = _mm512_mul_pd(a2,_mm512_mul_pd(b2,c2));
                          sths= _mm512_mul_pd(sth,sth);
                          cphi= xcos(phi);
                          cphis= _mm512_mul_pd(cphi,cphi);
                          trm1 = _mm512_fmadd_pd(_mm512_mul_pd(a2,sths),cphis,
                                                 _mm512_fmadd_pd(_mm512_mul_pd(b2,sths),
                                                                 _mm512_mul_pd(c2,cths)));
                          x0   = _mm512_mul_pd(trm1,trm1);
                          rcs  = _mm512_div_pd(num,x0);
                          return (rcs);
                }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5155_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pc,
                                              const double * __restrict __ATTR_ALIGN__(64) pth,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {
                                           
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d b   = _mm512_load_pd(&pb[0]);
                          register __m512d c   = _mm512_load_pd(&pc[0]);
                          register __m512d th  = _mm512_load_pd(&pth[0]);
                          register __m512d phi = _mm512_load_pd(&pphi[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a2,b2,c2,sphi,sphis,cth,cths,sth,sths,cphi,cphis;
                          register __m512d num,x0,trm1;
                          a2  = _mm512_mul_pd(a,a);
                          sphi= xsin(phi);
                          b2  = _mm512_mul_pd(b,b);
                          sphis= _mm512_mul_pd(sphi,sphi);
                          cth = xcos(th);
                          c2  = _mm512_mul_pd(c,c);
                          cths= _mm512_mul_pd(cth,cth);
                          sth = xsin(th);
                          num = _mm512_mul_pd(a2,_mm512_mul_pd(b2,c2));
                          sths= _mm512_mul_pd(sth,sth);
                          cphi= xcos(phi);
                          cphis= _mm512_mul_pd(cphi,cphi);
                          trm1 = _mm512_fmadd_pd(_mm512_mul_pd(a2,sths),cphis,
                                                 _mm512_fmadd_pd(_mm512_mul_pd(b2,sths),
                                                                 _mm512_mul_pd(c2,cths)));
                          x0   = _mm512_mul_pd(trm1,trm1);
                          rcs  = _mm512_div_pd(num,x0);
                          return (rcs);
                }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5155_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pc,
                                              const double * __restrict  pth,
                                              const double * __restrict  pphi) {
                                           
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d b   = _mm512_loadu_pd(&pb[0]);
                          register __m512d c   = _mm512_loadu_pd(&pc[0]);
                          register __m512d th  = _mm512_loadu_pd(&pth[0]);
                          register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a2,b2,c2,sphi,sphis,cth,cths,sth,sths,cphi,cphis;
                          register __m512d num,x0,trm1;
                          a2  = _mm512_mul_pd(a,a);
                          sphi= xsin(phi);
                          b2  = _mm512_mul_pd(b,b);
                          sphis= _mm512_mul_pd(sphi,sphi);
                          cth = xcos(th);
                          c2  = _mm512_mul_pd(c,c);
                          cths= _mm512_mul_pd(cth,cth);
                          sth = xsin(th);
                          num = _mm512_mul_pd(a2,_mm512_mul_pd(b2,c2));
                          sths= _mm512_mul_pd(sth,sth);
                          cphi= xcos(phi);
                          cphis= _mm512_mul_pd(cphi,cphi);
                          trm1 = _mm512_fmadd_pd(_mm512_mul_pd(a2,sths),cphis,
                                                 _mm512_fmadd_pd(_mm512_mul_pd(b2,sths),
                                                                 _mm512_mul_pd(c2,cths)));
                          x0   = _mm512_mul_pd(trm1,trm1);
                          rcs  = _mm512_div_pd(num,x0);
                          return (rcs);
                }


                    /*
                           High frequency solutions.
                           Perfectly conducting spheroids.
                           Bistatic case.
                           Formula 5.1-67
                        */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5167_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d th1,
                                            const __m512d phi1,
                                            const __m512d th2,
                                            const __m512d phi2) {

                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,a2,b2,b4,sth1,cphi1,sth2,cphi2,cth1,cth2;
                          register __m512d trm1,trm2,trm3,num,den,x0,strm1,strm2,strm3;
                          a2   = _mm512_mul_pd(a,a);
                          sth1 = xsin(phi1);
                          b2   = _mm512_mul_pd(b,b);
                          cphi1= xcos(phi1);
                          b4   = _mm512_mul_pd(b2,b2);
                          sth2 = xsin(th2);
                          x0   = _mm512_mul_pd(a2,b4);
                          cphi2= xcos(phi2);
                          num  = _mm512_mul_pd(_4pi,x0);
                          cth1 = xcos(th1);
                          trm1 = _mm512_fmadd_pd(sth1,cphi1,_mm512_mul_pd(sth2,cphi2));
                          strm1= _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          cth2 = xcos(th2);
                          trm2 = _mm512_fmadd_pd(sth1,sphi1,_mm512_mul_pd(sth2,sphi2));
                          strm2= _mm512_mul_pd(b2,_mm512_mul_pd(trm2,trm2));
                          trm3 = _mm512_mul_pd(cth1,cth2);
                          strm3= _mm512_mul_pd(b2,_mm512_mul_pd(trm3,trm3));
                          x0   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5167_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pth1,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pth2,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi2) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d b   = _mm512_load_pd(&pb[0]);
                          register __m512d th1 = _mm512_load_pd(&pth1[0]);
                          register __m512d phi1= _mm512_load_pd(&pphi1[0]);
                          register __m512d th2 = _mm512_load_pd(&pth2[0]);
                          register __m512d phi2= _mm512_load_pd(&pphi2[0]);
                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,a2,b2,b4,sth1,cphi1,sth2,cphi2,cth1,cth2;
                          register __m512d trm1,trm2,trm3,num,den,x0,strm1,strm2,strm3;
                          a2   = _mm512_mul_pd(a,a);
                          sth1 = xsin(phi1);
                          b2   = _mm512_mul_pd(b,b);
                          cphi1= xcos(phi1);
                          b4   = _mm512_mul_pd(b2,b2);
                          sth2 = xsin(th2);
                          x0   = _mm512_mul_pd(a2,b4);
                          cphi2= xcos(phi2);
                          num  = _mm512_mul_pd(_4pi,x0);
                          cth1 = xcos(th1);
                          trm1 = _mm512_fmadd_pd(sth1,cphi1,_mm512_mul_pd(sth2,cphi2));
                          strm1= _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          cth2 = xcos(th2);
                          trm2 = _mm512_fmadd_pd(sth1,sphi1,_mm512_mul_pd(sth2,sphi2));
                          strm2= _mm512_mul_pd(b2,_mm512_mul_pd(trm2,trm2));
                          trm3 = _mm512_mul_pd(cth1,cth2);
                          strm3= _mm512_mul_pd(b2,_mm512_mul_pd(trm3,trm3));
                          x0   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5167_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pth1,
                                              const double * __restrict  pphi1,
                                              const double * __restrict  pth2,
                                              const double * __restrict  pphi2) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d b   = _mm512_loadu_pd(&pb[0]);
                          register __m512d th1 = _mm512_loadu_pd(&pth1[0]);
                          register __m512d phi1= _mm512_loadu_pd(&pphi1[0]);
                          register __m512d th2 = _mm512_loadu_pd(&pth2[0]);
                          register __m512d phi2= _mm512_loadu_pd(&pphi2[0]);
                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,a2,b2,b4,sth1,cphi1,sth2,cphi2,cth1,cth2;
                          register __m512d trm1,trm2,trm3,num,den,x0,strm1,strm2,strm3;
                          a2   = _mm512_mul_pd(a,a);
                          sth1 = xsin(phi1);
                          b2   = _mm512_mul_pd(b,b);
                          cphi1= xcos(phi1);
                          b4   = _mm512_mul_pd(b2,b2);
                          sth2 = xsin(th2);
                          x0   = _mm512_mul_pd(a2,b4);
                          cphi2= xcos(phi2);
                          num  = _mm512_mul_pd(_4pi,x0);
                          cth1 = xcos(th1);
                          trm1 = _mm512_fmadd_pd(sth1,cphi1,_mm512_mul_pd(sth2,cphi2));
                          strm1= _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          cth2 = xcos(th2);
                          trm2 = _mm512_fmadd_pd(sth1,sphi1,_mm512_mul_pd(sth2,sphi2));
                          strm2= _mm512_mul_pd(b2,_mm512_mul_pd(trm2,trm2));
                          trm3 = _mm512_mul_pd(cth1,cth2);
                          strm3= _mm512_mul_pd(b2,_mm512_mul_pd(trm3,trm3));
                          x0   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                 }


                  /*
                           High frequency solutions.
                           Perfectly conducting spheroids.
                           Backscatter case.
                           Formula 5.1-68
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5168_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d phi) {

                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a2,b2,b4,num,x0;
                          register __m512d cphi,cphis,sphi,sphis,den;
                          a2   = _mm512_mul_pd(a,a);
                          cphi = xcos(phi);
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsin(phi);
                          b4   = _mm512_mul_pd(b2,b2);
                          num  = _mm512_mul_pd(pi,_mm512_mul_pd(a2,b4));
                          cphis= _mm512_mul_pd(cphi,cphi);
                          sphis= _mm512_mul_pd(sphi,sphi);
                          x0   = _mm512_fmadd_pd(a2,cphis,_mm512_mul_pd(b2,sphis));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5168_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pb,
                                            const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d b   = _mm512_load_pd(&pb[0]);
                          register __m512d phi = _mm512_load_pd(&pph1[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a2,b2,b4,num,x0;
                          register __m512d cphi,cphis,sphi,sphis,den;
                          a2   = _mm512_mul_pd(a,a);
                          cphi = xcos(phi);
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsin(phi);
                          b4   = _mm512_mul_pd(b2,b2);
                          num  = _mm512_mul_pd(pi,_mm512_mul_pd(a2,b4));
                          cphis= _mm512_mul_pd(cphi,cphi);
                          sphis= _mm512_mul_pd(sphi,sphi);
                          x0   = _mm512_fmadd_pd(a2,cphis,_mm512_mul_pd(b2,sphis));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5168_zmm8r8_u(const double * __restrict  pa,
                                            const double * __restrict  pb,
                                            const double * __restrict  pphi) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d b   = _mm512_loadu_pd(&pb[0]);
                          register __m512d phi = _mm512_loadu_pd(&pph1[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a2,b2,b4,num,x0;
                          register __m512d cphi,cphis,sphi,sphis,den;
                          a2   = _mm512_mul_pd(a,a);
                          cphi = xcos(phi);
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsin(phi);
                          b4   = _mm512_mul_pd(b2,b2);
                          num  = _mm512_mul_pd(pi,_mm512_mul_pd(a2,b4));
                          cphis= _mm512_mul_pd(cphi,cphi);
                          sphis= _mm512_mul_pd(sphi,sphi);
                          x0   = _mm512_fmadd_pd(a2,cphis,_mm512_mul_pd(b2,sphis));
                          den  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_div_pd(num,den);
                          return (rcs);
                }


                    /*
                           High frequency solutions.
                           Perfectly conducting spheroids.
                           Incidence along a the major axis -- the backscatter RCS.
                           Formula 5.1-69
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5169_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d k0) {

                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,a2,b2,b4,k0a,k02a2,arg,sarg;
                          register __m512d trm1,trm2,trm3,x0;  
                          a2   = _mm512_mul_pd(a,a);
                          k0a  = _mm512_mul_pd(k0,a);
                          arg  = _mm512_add_pd(k0a,k0a);
                          b2   = _mm512_mul_pd(b,b);
                          k02a2= _mm512_mul_pd(a2,_mm512_mul_pd(k0,k0));
                          sarg = xsin(arg);
                          b4   = _mm512_mul_pd(b2,b2);
                          trm1 = _mm512_div_pd(_mm512_mul_pd(pi,b4),a2);
                          x0   = _mm512_mul_pd(sarg,sarg);
                          trm2 = _mm512_sub_pd(_1,_mm512_div_pd(sarg,arg));
                          trm3 = _mm512_div_pd(x0,k02a2);
                          rcs  = _mm512_mul_pd(trm1,_mm512_add_pd(trm2,trm3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5169_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pb,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d b   = _mm512_load_pd(&pb[0]);
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,a2,b2,b4,k0a,k02a2,arg,sarg;
                          register __m512d trm1,trm2,trm3,x0;  
                          a2   = _mm512_mul_pd(a,a);
                          k0a  = _mm512_mul_pd(k0,a);
                          arg  = _mm512_add_pd(k0a,k0a);
                          b2   = _mm512_mul_pd(b,b);
                          k02a2= _mm512_mul_pd(a2,_mm512_mul_pd(k0,k0));
                          sarg = xsin(arg);
                          b4   = _mm512_mul_pd(b2,b2);
                          trm1 = _mm512_div_pd(_mm512_mul_pd(pi,b4),a2);
                          x0   = _mm512_mul_pd(sarg,sarg);
                          trm2 = _mm512_sub_pd(_1,_mm512_div_pd(sarg,arg));
                          trm3 = _mm512_div_pd(x0,k02a2);
                          rcs  = _mm512_mul_pd(trm1,_mm512_add_pd(trm2,trm3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5169_zmm8r8_u(const double * __restrict  pa,
                                            const double * __restrict pb,
                                            const double * __restrict  pk0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d b   = _mm512_loadu_pd(&pb[0]);
                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,a2,b2,b4,k0a,k02a2,arg,sarg;
                          register __m512d trm1,trm2,trm3,x0;  
                          a2   = _mm512_mul_pd(a,a);
                          k0a  = _mm512_mul_pd(k0,a);
                          arg  = _mm512_add_pd(k0a,k0a);
                          b2   = _mm512_mul_pd(b,b);
                          k02a2= _mm512_mul_pd(a2,_mm512_mul_pd(k0,k0));
                          sarg = xsin(arg);
                          b4   = _mm512_mul_pd(b2,b2);
                          trm1 = _mm512_div_pd(_mm512_mul_pd(pi,b4),a2);
                          x0   = _mm512_mul_pd(sarg,sarg);
                          trm2 = _mm512_sub_pd(_1,_mm512_div_pd(sarg,arg));
                          trm3 = _mm512_div_pd(x0,k02a2);
                          rcs  = _mm512_mul_pd(trm1,_mm512_add_pd(trm2,trm3));
                          return (rcs);
                }


                  /*
                           High frequency solutions.
                           Perfectly conducting spheroids.
                           Incidence along a the minor axis -- the backscatter RCS.
                           Formula 5.1-70

                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5170_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d k0) {

                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,a2,b2,k0b,k02b2,arg,sarg;
                          register __m512d trm1,trm2,trm3,x0;  
                          a2   = _mm512_mul_pd(a,a);
                          k0b  = _mm512_mul_pd(k0,b);
                          arg  = _mm512_add_pd(k0b,k0b);
                          b2   = _mm512_mul_pd(b,b);
                          k02b2= _mm512_mul_pd(b2,_mm512_mul_pd(k0,k0));
                          sarg = xsin(arg);
                          trm1 = _mm512_mul_pd(pi,a2);
                          x0   = _mm512_mul_pd(sarg,sarg);
                          trm2 = _mm512_sub_pd(_1,_mm512_div_pd(sarg,arg));
                          trm3 = _mm512_div_pd(x0,k02b2);
                          rcs  = _mm512_mul_pd(trm1,_mm512_add_pd(trm2,trm3));
                          return (rcs);
                }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5170_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pb,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d b   = _mm512_load_pd(&pb[0]);
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,a2,b2,k0b,k02b2,arg,sarg;
                          register __m512d trm1,trm2,trm3,x0;  
                          a2   = _mm512_mul_pd(a,a);
                          k0b  = _mm512_mul_pd(k0,b);
                          arg  = _mm512_add_pd(k0b,k0b);
                          b2   = _mm512_mul_pd(b,b);
                          k02b2= _mm512_mul_pd(b2,_mm512_mul_pd(k0,k0));
                          sarg = xsin(arg);
                          trm1 = _mm512_mul_pd(pi,a2);
                          x0   = _mm512_mul_pd(sarg,sarg);
                          trm2 = _mm512_sub_pd(_1,_mm512_div_pd(sarg,arg));
                          trm3 = _mm512_div_pd(x0,k02b2);
                          rcs  = _mm512_mul_pd(trm1,_mm512_add_pd(trm2,trm3));
                          return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5170_zmm8r8_u(const double * __restrict  pa,
                                            const double * __restrict  pb,
                                            const double * __restrict  pk0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d b   = _mm512_loadu_pd(&pb[0]);
                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,a2,b2,k0b,k02b2,arg,sarg;
                          register __m512d trm1,trm2,trm3,x0;  
                          a2   = _mm512_mul_pd(a,a);
                          k0b  = _mm512_mul_pd(k0,b);
                          arg  = _mm512_add_pd(k0b,k0b);
                          b2   = _mm512_mul_pd(b,b);
                          k02b2= _mm512_mul_pd(b2,_mm512_mul_pd(k0,k0));
                          sarg = xsin(arg);
                          trm1 = _mm512_mul_pd(pi,a2);
                          x0   = _mm512_mul_pd(sarg,sarg);
                          trm2 = _mm512_sub_pd(_1,_mm512_div_pd(sarg,arg));
                          trm3 = _mm512_div_pd(x0,k02b2);
                          rcs  = _mm512_mul_pd(trm1,_mm512_add_pd(trm2,trm3));
                          return (rcs);
                }


                   /*
                         Oblate spheroids.
                         Low frequency solutions.
                         Helper parameters: Ia,Ib
                         Formula 5.1-91
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d IaIb_f5191_zmm8r8(const __m512d a,
                                             const __m512d c) {

                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d IaIb,e,a2c,arg,asarg,e2m1s;
                          register __m512d trm1,trm2,trm3,e2m1,x0,x1;
                          e    = _mm512_div_pd(a,c);
                          a2c  = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                          e2   = _mm512_mul_pd(e,e);
                          trm1 = _mm512_div_pd(_2,a2c);
                          e2m1 = _mm512_sub_pd(e2,_1);
                          trm2 = _mm512_div_pd(_1,_mm512_add_pd(e2m1,e2m1));
                          x0   = _mm512_mul_pd(e2,e2m1);
                          e2m1s= _mm512_sqrt_pd(e2m1);
                          arg  = _mm512_div_pd(e2m1s,e);
                          x1   = _mm512_sqrt_pd(x0);
                          asarg= _mm512_asin_pd(arg);
                          trm3 = _mm512_rcp14_pd(x1);
                          asarg= _mm512_sub_pd(asarg,_1);
                          IaIb = _mm512_mul_pd(_mm512_mul_pd(trm1,trm2),
                                               _mm512_mul_pd(trm3,asarg));
                          return (IaIb);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d IaIb_f5191_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                               const double * __restrict __ATTR_ALIGN__(64) pc) {

                          register __m512d a = _mm512_load_pd(&pa[0]);
                          register __m512d c = _mm512_load_pd(&pc[0]);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d IaIb,e,a2c,arg,asarg,e2m1s;
                          register __m512d trm1,trm2,trm3,e2m1,x0,x1;
                          e    = _mm512_div_pd(a,c);
                          a2c  = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                          e2   = _mm512_mul_pd(e,e);
                          trm1 = _mm512_div_pd(_2,a2c);
                          e2m1 = _mm512_sub_pd(e2,_1);
                          trm2 = _mm512_div_pd(_1,_mm512_add_pd(e2m1,e2m1));
                          x0   = _mm512_mul_pd(e2,e2m1);
                          e2m1s= _mm512_sqrt_pd(e2m1);
                          arg  = _mm512_div_pd(e2m1s,e);
                          x1   = _mm512_sqrt_pd(x0);
                          asarg= _mm512_asin_pd(arg);
                          trm3 = _mm512_rcp14_pd(x1);
                          asarg= _mm512_sub_pd(asarg,_1);
                          IaIb = _mm512_mul_pd(_mm512_mul_pd(trm1,trm2),
                                               _mm512_mul_pd(trm3,asarg));
                          return (IaIb);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d IaIb_f5191_zmm8r8_u(const double * __restrict  pa,
                                               const double * __restrict  pc) {

                          register __m512d a = _mm512_loadu_pd(&pa[0]);
                          register __m512d c = _mm512_loadu_pd(&pc[0]);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d IaIb,e,a2c,arg,asarg,e2m1s;
                          register __m512d trm1,trm2,trm3,e2m1,x0,x1;
                          e    = _mm512_div_pd(a,c);
                          a2c  = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                          e2   = _mm512_mul_pd(e,e);
                          trm1 = _mm512_div_pd(_2,a2c);
                          e2m1 = _mm512_sub_pd(e2,_1);
                          trm2 = _mm512_div_pd(_1,_mm512_add_pd(e2m1,e2m1));
                          x0   = _mm512_mul_pd(e2,e2m1);
                          e2m1s= _mm512_sqrt_pd(e2m1);
                          arg  = _mm512_div_pd(e2m1s,e);
                          x1   = _mm512_sqrt_pd(x0);
                          asarg= _mm512_asin_pd(arg);
                          trm3 = _mm512_rcp14_pd(x1);
                          asarg= _mm512_sub_pd(asarg,_1);
                          IaIb = _mm512_mul_pd(_mm512_mul_pd(trm1,trm2),
                                               _mm512_mul_pd(trm3,asarg));
                          return (IaIb);
                 }


                     /*
                         Oblate spheroids.
                         Low frequency solutions.
                         Helper parameters: Ia,Ib
                         Formula 5.1-92
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d Ic_f5192_zmm8r8(const __m512d a,
                                             const __m512d c) {

                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d Ic,e,a2c,arg,asarg,e2m1s;
                          register __m512d trm1,trm2,trm3,e2m1,x0,x1;
                          e    = _mm512_div_pd(a,c);
                          a2c  = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                          e2   = _mm512_mul_pd(e,e);
                          trm1 = _mm512_div_pd(_2,a2c);
                          e2m1 = _mm512_sub_pd(e2,_1);
                          trm2 = _mm512_div_pd(e2,e2m1,e2m1);
                          x0   = _mm512_sub_pd(_1,e2m1);
                          e2m1s= _mm512_sqrt_pd(e2m1);
                          arg  = _mm512_div_pd(e2m1s,e);
                          x1   = _mm512_sqrt_pd(x0);
                          asarg= _mm512_asin_pd(arg);
                          trm3 = _mm512_rcp14_pd(x1);
                          Ic   = _mm512_mul_pd(_mm512_mul_pd(trm1,trm2),
                                               _mm512_mul_pd(trm3,asarg));
                          return (Ic);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d Ic_f5192_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pc) {

                          register __m512d a = _mm512_load_pd(&pa[0]);
                          register __m512d c = _mm512_load_pd(&pc[0]);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d Ic,e,a2c,arg,asarg,e2m1s;
                          register __m512d trm1,trm2,trm3,e2m1,x0,x1;
                          e    = _mm512_div_pd(a,c);
                          a2c  = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                          e2   = _mm512_mul_pd(e,e);
                          trm1 = _mm512_div_pd(_2,a2c);
                          e2m1 = _mm512_sub_pd(e2,_1);
                          trm2 = _mm512_div_pd(e2,e2m1,e2m1);
                          x0   = _mm512_sub_pd(_1,e2m1);
                          e2m1s= _mm512_sqrt_pd(e2m1);
                          arg  = _mm512_div_pd(e2m1s,e);
                          x1   = _mm512_sqrt_pd(x0);
                          asarg= _mm512_asin_pd(arg);
                          trm3 = _mm512_rcp14_pd(x1);
                          Ic   = _mm512_mul_pd(_mm512_mul_pd(trm1,trm2),
                                               _mm512_mul_pd(trm3,asarg));
                          return (Ic);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d Ic_f5192_zmm8r8_u(const double * __restrict  pa,
                                           const double * __restrict  pc) {

                          register __m512d a = _mm512_loadu_pd(&pa[0]);
                          register __m512d c = _mm512_loadu_pd(&pc[0]);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d Ic,e,a2c,arg,asarg,e2m1s;
                          register __m512d trm1,trm2,trm3,e2m1,x0,x1;
                          e    = _mm512_div_pd(a,c);
                          a2c  = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                          e2   = _mm512_mul_pd(e,e);
                          trm1 = _mm512_div_pd(_2,a2c);
                          e2m1 = _mm512_sub_pd(e2,_1);
                          trm2 = _mm512_div_pd(e2,e2m1,e2m1);
                          x0   = _mm512_sub_pd(_1,e2m1);
                          e2m1s= _mm512_sqrt_pd(e2m1);
                          arg  = _mm512_div_pd(e2m1s,e);
                          x1   = _mm512_sqrt_pd(x0);
                          asarg= _mm512_asin_pd(arg);
                          trm3 = _mm512_rcp14_pd(x1);
                          Ic   = _mm512_mul_pd(_mm512_mul_pd(trm1,trm2),
                                               _mm512_mul_pd(trm3,asarg));
                          return (Ic);
                 }


                   /*
                         Oblate spheroids (perfectly conducting);
                         Low frequency solutions.
                         Backscatter RCS.
                         Formula 5.1-89
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5189_zmm8r8(const __m512d a,
                                            const __m512d c,
                                            const __m512d tht,
                                            const __m512d k0) {

                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d ptrm= _mm512_set1_pd(5.585053606381854646155810459164);
                          register __m512d rcs,Ia,Ic,k04,x0,stht,ctht,stht2,ctht2,trm1;
                          register __m512d trm2,trm3,trm23,trm1,IaIc,x1;
                          x0   = _mm512_mul_pd(k0,k0);
                          Ia   = IaIb_f5191_zmm8r8(a,c);
                          k04  = _mm512_mul_pd(x0,x0);
                          Ic   = Ic_f5192_zmm8r8(a,c);
                          trm1 = _mm512_mul_pd(k04,ptrm);
                          stht = xsin(tht); 
                          IaIc = _mm512_add_pd(Ia,Ic);
                          ctht = xcos(tht);
                          stht2= _mm512_mul_pd(stht,stht);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          x0   = _mm512_fmadd_pd(stht2,hlf,_1);
                          trm3 = _mm512_div_pd(ctht2,IaIc);
                          trm2 = _mm512_div_pd(x0,Ia);
                          x1   = _mm512_add_pd(trm2,trm3);
                          trm23= _mm512_mul_pd(x1,x1);
                          rcs  = _mm512_mul_pd(trm1,trm23);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5189_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pc,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0 ) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d c   = _mm512_load_pd(&pc[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d ptrm= _mm512_set1_pd(5.585053606381854646155810459164);
                          register __m512d rcs,Ia,Ic,k04,x0,stht,ctht,stht2,ctht2,trm1;
                          register __m512d trm2,trm3,trm23,trm1,IaIc,x1;
                          x0   = _mm512_mul_pd(k0,k0);
                          Ia   = IaIb_f5191_zmm8r8(a,c);
                          k04  = _mm512_mul_pd(x0,x0);
                          Ic   = Ic_f5192_zmm8r8(a,c);
                          trm1 = _mm512_mul_pd(k04,ptrm);
                          stht = xsin(tht); 
                          IaIc = _mm512_add_pd(Ia,Ic);
                          ctht = xcos(tht);
                          stht2= _mm512_mul_pd(stht,stht);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          x0   = _mm512_fmadd_pd(stht2,hlf,_1);
                          trm3 = _mm512_div_pd(ctht2,IaIc);
                          trm2 = _mm512_div_pd(x0,Ia);
                          x1   = _mm512_add_pd(trm2,trm3);
                          trm23= _mm512_mul_pd(x1,x1);
                          rcs  = _mm512_mul_pd(trm1,trm23);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5189_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pc,
                                              const double * __restrict  ptht,
                                              const double * __restrict  pk0 ) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d c   = _mm512_loadu_pd(&pc[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d ptrm= _mm512_set1_pd(5.585053606381854646155810459164);
                          register __m512d rcs,Ia,Ic,k04,x0,stht,ctht,stht2,ctht2,trm1;
                          register __m512d trm2,trm3,trm23,trm1,IaIc,x1;
                          x0   = _mm512_mul_pd(k0,k0);
                          Ia   = IaIb_f5191_zmm8r8(a,c);
                          k04  = _mm512_mul_pd(x0,x0);
                          Ic   = Ic_f5192_zmm8r8(a,c);
                          trm1 = _mm512_mul_pd(k04,ptrm);
                          stht = xsin(tht); 
                          IaIc = _mm512_add_pd(Ia,Ic);
                          ctht = xcos(tht);
                          stht2= _mm512_mul_pd(stht,stht);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          x0   = _mm512_fmadd_pd(stht2,hlf,_1);
                          trm3 = _mm512_div_pd(ctht2,IaIc);
                          trm2 = _mm512_div_pd(x0,Ia);
                          x1   = _mm512_add_pd(trm2,trm3);
                          trm23= _mm512_mul_pd(x1,x1);
                          rcs  = _mm512_mul_pd(trm1,trm23);
                          return (rcs);
                 }


                     /*
                         Oblate spheroids (perfectly conducting);
                         Low frequency solutions.
                         Backscatter RCS.
                         Formula 5.1-90
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5190_zmm8r8(const __m512d a,
                                            const __m512d c,
                                            const __m512d tht,
                                            const __m512d k0) {

                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d ptrm= _mm512_set1_pd(5.585053606381854646155810459164);
                          register __m512d rcs,Ia,Ic,k04,x0,stht,ctht,stht2,ctht2,trm1;
                          register __m512d trm2,trm3,strm,trm1,IaIc,x1,inv;
                          x0   = _mm512_mul_pd(k0,k0);
                          Ia   = IaIb_f5191_zmm8r8(a,c);
                          k04  = _mm512_mul_pd(x0,x0);
                          Ic   = Ic_f5192_zmm8r8(a,c);
                          trm1 = _mm512_mul_pd(k04,ptrm);
                          stht = xsin(tht); 
                          IaIc = _mm512_add_pd(Ia,Ic);
                          inv  = _mm512_rcp14_pd(IaIc);
                          ctht = xcos(tht);
                          stht2= _mm512_mul_pd(stht,stht);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          trm2 = _mm512_div_pd(ctht2,Ia);
                          trm3 = _mm512_div_pd(stht3,Ic);
                          strm = _mm512_add_pd(trm2,trm3);
                          x1   = _mm512_mul_pd(strm,strm);
                          rcs  = _mm512_mul_pd(trm1,x1);
                          return (rcs);
                 } 


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5190_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pc,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d c   = _mm512_load_pd(&pc[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d ptrm= _mm512_set1_pd(5.585053606381854646155810459164);
                          register __m512d rcs,Ia,Ic,k04,x0,stht,ctht,stht2,ctht2,trm1;
                          register __m512d trm2,trm3,strm,trm1,IaIc,x1,inv;
                          x0   = _mm512_mul_pd(k0,k0);
                          Ia   = IaIb_f5191_zmm8r8(a,c);
                          k04  = _mm512_mul_pd(x0,x0);
                          Ic   = Ic_f5192_zmm8r8(a,c);
                          trm1 = _mm512_mul_pd(k04,ptrm);
                          stht = xsin(tht); 
                          IaIc = _mm512_add_pd(Ia,Ic);
                          inv  = _mm512_rcp14_pd(IaIc);
                          ctht = xcos(tht);
                          stht2= _mm512_mul_pd(stht,stht);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          trm2 = _mm512_div_pd(ctht2,Ia);
                          trm3 = _mm512_div_pd(stht3,Ic);
                          strm = _mm512_add_pd(trm2,trm3);
                          x1   = _mm512_mul_pd(strm,strm);
                          rcs  = _mm512_mul_pd(trm1,x1);
                          return (rcs);
                 } 


                   /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          E-field (theta), perpendicular, formula 5.1-83
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5183_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict ESr,
                                          __m512d * __restrict ESi) {

                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sph2,ear,eai,cer,cei,epsrm1,epsim1;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(ctht2,stht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(ctht1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       t0r = _mm512_sub_pd(mul1r,mul2r);
                       t0i = _mm512_sub_pd(mul1i,mul2i);
                       cmul_zmm8r8(facr,faci,t0r,t0i,*ESr,*ESi);
                }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5183_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          double * __restrict __ATTR_ALIGN__(64) ESr,
                                          double * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512d  k0    = _mm512_load_pd(&pk0[0]);
                       register __m512d  r     = _mm512_load_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_load_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_load_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_load_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_load_pd(&pmui[0]);
                       register __m512d  a     = _mm512_load_pd(&pa[0]);
                       register __m512d  c     = _mm512_load_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_load_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_load_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_load_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sph2,ear,eai,cer,cei,epsrm1,epsim1,resr,resi;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(ctht2,stht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(ctht1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       t0r = _mm512_sub_pd(mul1r,mul2r);
                       t0i = _mm512_sub_pd(mul1i,mul2i);
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_store_pd(&ESr[0], resr);
                       _mm512_store_pd(&ESi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5183_zmm8r8_u(const double * __restrict  pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht1,
                                          const double * __restrict  ptht2,
                                          const double * __restrict  pphi2,
                                          double * __restrict  ESr,
                                          double * __restrict  ESi) {

                       register __m512d  k0    = _mm512_loadu_pd(&pk0[0]);
                       register __m512d  r     = _mm512_loadu_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_loadu_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_loadu_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_loadu_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_loadu_pd(&pmui[0]);
                       register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                       register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_loadu_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_loadu_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_loadu_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sph2,ear,eai,cer,cei,epsrm1,epsim1,resr,resi;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(ctht2,stht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(ctht1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       t0r = _mm512_sub_pd(mul1r,mul2r);
                       t0i = _mm512_sub_pd(mul1i,mul2i);
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_storeu_pd(&ESr[0], resr);
                       _mm512_storeu_pd(&ESi[0], resi);
                }


                  /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          H-field (phi), perpendicular formula 5.1-83
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5183_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict HSr,
                                          __m512d * __restrict HSi) {

                    ESth_f5183_zmm8r8(k0,r,epsr,epsi,mur,mui,a,
                                       c,tht1,tht2,phi2,*HSr,*HSi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5183_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          double * __restrict __ATTR_ALIGN__(64) HSr,
                                          double * __restrict __ATTR_ALIGN__(64) HSi) {

                         ESth_f5183_zmm8r8_a(pk0,pr,pepsr,pepsi,pmur,pmui,
                                              pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
                }


                 __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5183_zmm8r8_u(const double * __restrict  pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht1,
                                          const double * __restrict  ptht2,
                                          const double * __restrict  pphi2,
                                          double * __restrict  HSr,
                                          double * __restrict  HSi) {

                         ESth_f5183_zmm8r8_u(pk0,pr,pepsr,pepsi,pmur,pmui,
                                              pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
                }


                 /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          E-field (phi), perpendicular formula 5.1-84

                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5184_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict ESr,
                                          __m512d * __restrict ESi) {

                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sth2,sth1,ear,eai,cer,cei,epsrm1,epsim1,cphi2,t2r,t2i;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c,den3r,den3i;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       cphi2  = xcos(phi2);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       sth2   = xsin(tht2);
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(cth1,_mm512_mul_pd(ctht2,cphi2));
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(stht1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(cphi2,den1r);
                       t1r    = _mm512_div_pd(x0,den2r);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       den3r  = _mm512_fmadd_pd(murm1,Ic,_2a2c);
                       t2r    = _mm512_div_pd(x1,den3r);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0i    = _mm512_div_pd(cphi2,den1i);
                       t1i    = _mm512_div_pd(x0,den2i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       den3i  = _mm512_fmadd_pd(muim1,Ic,_2a2c);
                       t2i    = _mm512_div_pd(x1,den3i);
                       cmul_zmm8r8(murm1,muim1,t2r,t2i,&mul3r,&mul3i);
                       t0r = _mm512_sub_pd(mul1r,_mm512_sub_pd(mul2r,mul3r));
                       t0i = _mm512_sub_pd(mul1i,_mm512_sub_pd(mul2i,mul3i));
                       cmul_zmm8r8(facr,faci,t0r,t0i,*ESr,*ESi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5184_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          double * __restrict __ATTR_ALIGN__(64) ESr,
                                          double * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512d  k0    = _mm512_load_pd(&pk0[0]);
                       register __m512d  r     = _mm512_load_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_load_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_load_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_load_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_load_pd(&pmui[0]);
                       register __m512d  a     = _mm512_load_pd(&pa[0]);
                       register __m512d  c     = _mm512_load_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_load_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_load_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_load_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht1,ctht2,facr,faci,resr,resi;
                       register __m512d sth2,sth1,ear,eai,cer,cei,epsrm1,epsim1,cphi2,t2r,t2i;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c,den3r,den3i;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       cphi2  = xcos(phi2);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       sth2   = xsin(tht2);
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(cth1,_mm512_mul_pd(ctht2,cphi2));
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(stht1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(cphi2,den1r);
                       t1r    = _mm512_div_pd(x0,den2r);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       den3r  = _mm512_fmadd_pd(murm1,Ic,_2a2c);
                       t2r    = _mm512_div_pd(x1,den3r);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0i    = _mm512_div_pd(cphi2,den1i);
                       t1i    = _mm512_div_pd(x0,den2i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       den3i  = _mm512_fmadd_pd(muim1,Ic,_2a2c);
                       t2i    = _mm512_div_pd(x1,den3i);
                       cmul_zmm8r8(murm1,muim1,t2r,t2i,&mul3r,&mul3i);
                       t0r = _mm512_sub_pd(mul1r,_mm512_sub_pd(mul2r,mul3r));
                       t0i = _mm512_sub_pd(mul1i,_mm512_sub_pd(mul2i,mul3i));
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_store_pd(&ESr[0], resr);
                       _mm512_store_pd(&ESi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5184_zmm8r8_u(const double * __restrict  pk0,
                                          const double * __restrict pr,
                                          const double * __restrict pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict ptht1,
                                          const double * __restrict  ptht2,
                                          const double * __restrict  pphi2,
                                          double * __restrict  ESr,
                                          double * __restrict  ESi) {

                       register __m512d  k0    = _mm512_loadu_pd(&pk0[0]);
                       register __m512d  r     = _mm512_loadu_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_loadu_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_loadu_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_loadu_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_loadu_pd(&pmui[0]);
                       register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                       register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_loadu_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_loadu_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_load_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht1,ctht2,facr,faci,resr,resi;
                       register __m512d sth2,sth1,ear,eai,cer,cei,epsrm1,epsim1,cphi2,t2r,t2i;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c,den3r,den3i;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       cphi2  = xcos(phi2);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       sth2   = xsin(tht2);
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(cth1,_mm512_mul_pd(ctht2,cphi2));
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(stht1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(cphi2,den1r);
                       t1r    = _mm512_div_pd(x0,den2r);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       den3r  = _mm512_fmadd_pd(murm1,Ic,_2a2c);
                       t2r    = _mm512_div_pd(x1,den3r);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0i    = _mm512_div_pd(cphi2,den1i);
                       t1i    = _mm512_div_pd(x0,den2i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       den3i  = _mm512_fmadd_pd(muim1,Ic,_2a2c);
                       t2i    = _mm512_div_pd(x1,den3i);
                       cmul_zmm8r8(murm1,muim1,t2r,t2i,&mul3r,&mul3i);
                       t0r = _mm512_sub_pd(mul1r,_mm512_sub_pd(mul2r,mul3r));
                       t0i = _mm512_sub_pd(mul1i,_mm512_sub_pd(mul2i,mul3i));
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_storeu_pd(&ESr[0], resr);
                       _mm512_storeu_pd(&ESi[0], resi);
                }


                  /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          H-field (theta), perpendicular formula 5.1-84
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSth_f5184_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict HSr,
                                          __m512d * __restrict HSi) {

                        ESph_f5184_zmm8r8(k0,r,epsr,epsi,mur,mui,
                                           a,c,tht1,tht2,phi2,*HSr,*HSi);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5184_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          double * __restrict __ATTR_ALIGN__(64) HSr,
                                          double * __restrict __ATTR_ALIGN__(64) HSi) {

                         ESph_f5184_zmm8r8_a(pk0,pr,pepsr,pepsi,pmur,pmui,
                                              pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5184_zmm8r8_u(const double * __restrict pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht1,
                                          const double * __restrict  ptht2,
                                          const double * __restrict  pphi2,
                                          double * __restrict  HSr,
                                          double * __restrict  HSi) {

                         ESph_f5184_zmm8r8_u(pk0,pr,pepsr,pepsi,pmur,pmui,
                                              pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
                }


                     /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          E-field (theta), parallel, formula 5.1-85

                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5185_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict ESr,
                                          __m512d * __restrict ESi) {

                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sth2,sth1,ear,eai,cer,cei,epsrm1,epsim1,cphi2,t2r,t2i;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c,den3r,den3i;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       cphi2  = xcos(phi2);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       sth2   = xsin(tht2);
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(cth1,_mm512_mul_pd(ctht2,cphi2));
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(stht1,stht2);
                       den2r  = _mm512_fmadd_pd(epsrm1,Ic,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       den3r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       t2r    = _mm512_div_pd(cphi2,den3r);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(epsim1,Ic,_2a2c);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t1r,t1i,&mul2r,&mul2i);
                       den3i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t2i    = _mm512_div_pd(cphi2,den3i);
                       cmul_zmm8r8(murm1,muim1,t2r,t2i,&mul3r,&mul3i);
                       t0r = _mm512_add_pd(mul1r,_mm512_sub_pd(mul2r,mul3r));
                       t0i = _mm512_add_pd(mul1i,_mm512_sub_pd(mul2i,mul3i));
                       cmul_zmm8r8(facr,faci,t0r,t0i,*ESr,*ESi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5185_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          double * __restrict __ATTR_ALIGN__(64) ESr,
                                          double * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512d  k0    = _mm512_load_pd(&pk0[0]);
                       register __m512d  r     = _mm512_load_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_load_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_load_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_load_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_load_pd(&pmui[0]);
                       register __m512d  a     = _mm512_load_pd(&pa[0]);
                       register __m512d  c     = _mm512_load_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_load_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_load_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_load_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht1,ctht2,facr,faci,resr,resi;
                       register __m512d sth2,sth1,ear,eai,cer,cei,epsrm1,epsim1,cphi2,t2r,t2i;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c,den3r,den3i;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       cphi2  = xcos(phi2);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       sth2   = xsin(tht2);
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(cth1,_mm512_mul_pd(ctht2,cphi2));
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(stht1,stht2);
                       den2r  = _mm512_fmadd_pd(epsrm1,Ic,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       den3r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       t2r    = _mm512_div_pd(cphi2,den3r);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(epsim1,Ic,_2a2c);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t1r,t1i,&mul2r,&mul2i);
                       den3i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t2i    = _mm512_div_pd(cphi2,den3i);
                       cmul_zmm8r8(murm1,muim1,t2r,t2i,&mul3r,&mul3i);
                       t0r = _mm512_add_pd(mul1r,_mm512_sub_pd(mul2r,mul3r));
                       t0i = _mm512_add_pd(mul1i,_mm512_sub_pd(mul2i,mul3i));
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_store_pd(&ESr[0], resr);
                       _mm512_store_pd(&ESi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5185_zmm8r8_u(const double * __restrict pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht1,
                                          const double * __restrict  ptht2,
                                          const double * __restrict  pphi2,
                                          double * __restrict  ESr,
                                          double * __restrict  ESi) {

                       register __m512d  k0    = _mm512_loadu_pd(&pk0[0]);
                       register __m512d  r     = _mm512_loadu_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_loadu_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_loadu_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_loadu_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_loadu_pd(&pmui[0]);
                       register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                       register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_loadu_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_loadu_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_loadu_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht1,ctht2,facr,faci,resr,resi;
                       register __m512d sth2,sth1,ear,eai,cer,cei,epsrm1,epsim1,cphi2,t2r,t2i;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c,den3r,den3i;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       cphi2  = xcos(phi2);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       sth2   = xsin(tht2);
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(cth1,_mm512_mul_pd(ctht2,cphi2));
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(stht1,stht2);
                       den2r  = _mm512_fmadd_pd(epsrm1,Ic,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       den3r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       t2r    = _mm512_div_pd(cphi2,den3r);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(epsim1,Ic,_2a2c);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t1r,t1i,&mul2r,&mul2i);
                       den3i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t2i    = _mm512_div_pd(cphi2,den3i);
                       cmul_zmm8r8(murm1,muim1,t2r,t2i,&mul3r,&mul3i);
                       t0r = _mm512_add_pd(mul1r,_mm512_sub_pd(mul2r,mul3r));
                       t0i = _mm512_add_pd(mul1i,_mm512_sub_pd(mul2i,mul3i));
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_store_pd(&ESr[0], resr);
                       _mm512_store_pd(&ESi[0], resi);
                }


                   /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          H-field (phi), parallel, formula 5.1-85

                    */



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5185_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict HSr,
                                          __m512d * __restrict HSi) {

                        ESth_f5185_zmm8r8(k0,r,epsr,epsi,mur,mui,
                                           a,c,tht1,tht2,phi2,*HSr,*HSi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5185_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                             const double * __restrict __ATTR_ALIGN__(64) pmur,
                                             const double * __restrict __ATTR_ALIGN__(64) pmui,
                                             const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) pc,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                             const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                             double * __restrict __ATTR_ALIGN__(64) HSr,
                                             double * __restrict __ATTR_ALIGN__(64) HSi) {

                       ESth_f5185_zmm8r8_a(pk0,pr,pepsr,pepsi,pmur,pmui,
                                            pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSph_f5185_zmm8r8_u(const double * __restrict pk0,
                                             const double * __restrict  pr,
                                             const double * __restrict  pepsr,
                                             const double * __restrict  pepsi,
                                             const double * __restrict  pmur,
                                             const double * __restrict  pmui,
                                             const double * __restrict  pa,
                                             const double * __restrict  pc,
                                             const double * __restrict  ptht1,
                                             const double * __restrict  ptht2,
                                             const double * __restrict  pphi2,
                                             double * __restrict  HSr,
                                             double * __restrict  HSi) {

                       ESth_f5185_zmm8r8_u(pk0,pr,pepsr,pepsi,pmur,pmui,
                                            pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
                }


                  /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          E-field (phi), parallel, formula 5.1-86


                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5186_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict ESr,
                                          __m512d * __restrict ESi) {

                       const __m512d _23 = _mm512_set1_pd(-0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sph2,ear,eai,cer,cei,epsrm1,epsim1;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(ctht1,stht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(ctht2,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       t0r = _mm512_sub_pd(mul1r,mul2r);
                       t0i = _mm512_sub_pd(mul1i,mul2i);
                       cmul_zmm8r8(facr,faci,t0r,t0i,*ESr,*ESi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5186_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          double * __restrict __ATTR_ALIGN__(64) ESr,
                                          double * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512d  k0    = _mm512_load_pd(&pk0[0]);
                       register __m512d  r     = _mm512_load_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_load_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_load_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_load_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_load_pd(&pmui[0]);
                       register __m512d  a     = _mm512_load_pd(&pa[0]);
                       register __m512d  c     = _mm512_load_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_load_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_load_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_load_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(-0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sph2,ear,eai,cer,cei,epsrm1,epsim1,resr,resi;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(ctht1,stht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(ctht2,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       t0r = _mm512_sub_pd(mul1r,mul2r);
                       t0i = _mm512_sub_pd(mul1i,mul2i);
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_store_pd(&ESr[0], resr);
                       _mm512_store_pd(&ESi[0], resi);
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5186_zmm8r8_u(const double * __restrict  pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht1,
                                          const double * __restrict  ptht2,
                                          const double * __restrict  pphi2,
                                          double * __restrict  ESr,
                                          double * __restrict _ ESi) {

                       register __m512d  k0    = _mm512_loadu_pd(&pk0[0]);
                       register __m512d  r     = _mm512_loadu_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_loadu_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_loadu_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_loadu_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_loadu_pd(&pmui[0]);
                       register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                       register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_loadu_pd(&ptht1[0]); 
                       register __m512d  tht2  = _mm512_loadu_pd(&ptht2[0]);
                       register __m512d  phi2  = _mm512_loadu_pd(&pphi2[0]);
                       const __m512d _23 = _mm512_set1_pd(-0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,a2c,ctht1,ctht2,facr,faci;
                       register __m512d sph2,ear,eai,cer,cei,epsrm1,epsim1,resr,resi;
                       register __m512d murm1,muim1,x0,x1,t0r,t0i,t1r,t1i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,mul1r,mul1i,mul2r,mul2i;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht2  = xcos(tht2);
                       murm1  = _mm512_sub_pd(mur,_1)
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht2  = xsin(tht2);
                       k0r    = _mm512_mul_pd(k0,r);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ctht1  = xcos(tht1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       x0     = _mm512_mul_pd(ctht1,stht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       x1     = _mm512_mul_pd(ctht2,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       t0r    = _mm512_div_pd(x0,den1r);
                       t1r    = _mm512_div_pd(x1,den2r);
                       t0i    = _mm512_div_pd(x0,den1i);
                       t1i    = _mm512_div_pd(x1,den2i);
                       cmul_zmm8r8(epsrm1,epsim1,t0r,t0i,&mul1r,&mul1i);
                       cmul_zmm8r8(murm1,muim1,t1r,t1i,&mul2r,&mul2i);
                       t0r = _mm512_sub_pd(mul1r,mul2r);
                       t0i = _mm512_sub_pd(mul1i,mul2i);
                       cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_storeu_pd(&ESr[0], resr);
                       _mm512_storeu_pd(&ESi[0], resi);
                }


                      /*
                          Low-frequency oblate spheroid Rayleigh
                          bistatic scattered fields.
                          H-field (theta), parallel, formula 5.1-86


                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSth_f5186_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d phi2,
                                          __m512d * __restrict HSr,
                                          __m512d * __restrict HSi) {

                        ESph_f5186_zmm8r8(k0,r,epsr,epsi,mur,mui,
                                           a,c,tht1,tht2,phi3,*HSr,*HSi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSth_f5186_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          double * __restrict __ATTR_ALIGN__(64) HSr,
                                          double * __restrict __ATTR_ALIGN__(64) HSi) {

                        ESph_f5186_zmm8r8_a(pk0,pr,pepsr,pepsi,pmur,pmui,
                                             pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HSth_f5186_zmm8r8_u(const double * __restrict pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht1,
                                          const double * __restrict  ptht2,
                                          const double * __restrict  pphi2,
                                          double * __restrict  HSr,
                                          double * __restrict  HSi) {

                        ESph_f5186_zmm8r8_u(pk0,pr,pepsr,pepsi,pmur,pmui,
                                             pa,pc,ptht1,ptht2,pphi2,HSr,HSi);
               }


                    /*
                          Low-frequency oblate spheroid Rayleigh
                          Backscattered fields.
                          E-field (phi), perpendicular, formula 5.1-87         
 
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5187_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht,
                                          __m512d * __restrict ESr,
                                          __m512d * __restrict ESi) {

                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht,stht,facr,faci;
                       register __m512d ear,eai,cer,cei,epsrm1,epsim1,t2r,t2i,den3r,den3i;
                       register __m512d murm1,muim1,x0,t0r,t0i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,div1r,div1i,div2r,div2i,div3r,div3i;
                       register __m512d num1r,num1i,num2r,num2i,ctht2,stht2;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht   = xcos(tht);
                       murm1  = _mm512_sub_pd(mur,_1)
                       ctht2  = _mm512_mul_pd(ctht,ctht);
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht   = xsin(tht);
                       k0r    = _mm512_mul_pd(k0,r);
                       stht2  = _mm512_mul_pd(stht,stht);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       num1r  = _mm512_mul_pd(murm1,ctht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       num2r  = _mm512_mul_pd(murm1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den3r  = _mm512_fmadd_pd(murm1,Ic,_2a2c);
                       num1i  = _mm512_mul_pd(muim1,ctht2);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       num2i  = _mm512_mul_pd(muim1,stht2);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       den3i  = _mm512_fmadd_pd(muim1,Ic,_2a2c);
                       cdiv_zmm8r8(epsrm1,epsim1,den1r,den1i,&div1r,&div1i);
                       cdiv_zmm8r8(num1r,num1i,den2r,den2i,&div2r,&div2i);
                       cdiv_zmm8r8(num2r,num2i,den3r,den3i,&div2r,&div2i);
                       t0r = _mm512_sub_pd(div1r,_mm512_sub_pd(div2r,div3r));
                       t0r = _mm512_sub_pd(div1i,_mm512_sub_pd(div2i,div3i));
                       cexp_zmm8r8(facr,faci,t0r,t0i,*ESr,*ESi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5187_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht,
                                          double * __restrict __ATTR_ALIGN__(64) ESr,
                                          double * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512d  k0    = _mm512_load_pd(&pk0[0]);
                       register __m512d  r     = _mm512_load_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_load_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_load_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_load_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_load_pd(&pmui[0]);
                       register __m512d  a     = _mm512_load_pd(&pa[0]);
                       register __m512d  c     = _mm512_load_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_load_pd(&ptht[0]); 
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht,stht,facr,faci;
                       register __m512d ear,eai,cer,cei,epsrm1,epsim1,t2r,t2i,den3r,den3i;
                       register __m512d murm1,muim1,x0,t0r,t0i,_2a2c,resr,resi;
                       register __m512d den1r,den1i,den2r,den2i,div1r,div1i,div2r,div2i,div3r,div3i;
                       register __m512d num1r,num1i,num2r,num2i,ctht2,stht2;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht   = xcos(tht);
                       murm1  = _mm512_sub_pd(mur,_1)
                       ctht2  = _mm512_mul_pd(ctht,ctht);
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht   = xsin(tht);
                       k0r    = _mm512_mul_pd(k0,r);
                       stht2  = _mm512_mul_pd(stht,stht);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       num1r  = _mm512_mul_pd(murm1,ctht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       num2r  = _mm512_mul_pd(murm1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den3r  = _mm512_fmadd_pd(murm1,Ic,_2a2c);
                       num1i  = _mm512_mul_pd(muim1,ctht2);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       num2i  = _mm512_mul_pd(muim1,stht2);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       den3i  = _mm512_fmadd_pd(muim1,Ic,_2a2c);
                       cdiv_zmm8r8(epsrm1,epsim1,den1r,den1i,&div1r,&div1i);
                       cdiv_zmm8r8(num1r,num1i,den2r,den2i,&div2r,&div2i);
                       cdiv_zmm8r8(num2r,num2i,den3r,den3i,&div2r,&div2i);
                       t0r = _mm512_sub_pd(div1r,_mm512_sub_pd(div2r,div3r));
                       t0r = _mm512_sub_pd(div1i,_mm512_sub_pd(div2i,div3i));
                       cexp_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_store_pd(&ESr[0], resr);
                       _mm512_store_pd(&ESi[0], resi);
               }

              
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f5187_zmm8r8_u(const double * __restrict  pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict pmur,
                                          const double * __restrict  pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht,
                                          double * __restrict  ESr,
                                          double * __restrict  ESi) {

                       register __m512d  k0    = _mm512_loadu_pd(&pk0[0]);
                       register __m512d  r     = _mm512_loadu_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_loadu_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_loadu_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_loadu_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_loadu_pd(&pmui[0]);
                       register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                       register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_loadu_pd(&ptht[0]); 
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht,stht,facr,faci;
                       register __m512d ear,eai,cer,cei,epsrm1,epsim1,t2r,t2i,den3r,den3i;
                       register __m512d murm1,muim1,x0,t0r,t0i,_2a2c,resr,resi;
                       register __m512d den1r,den1i,den2r,den2i,div1r,div1i,div2r,div2i,div3r,div3i;
                       register __m512d num1r,num1i,num2r,num2i,ctht2,stht2;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht   = xcos(tht);
                       murm1  = _mm512_sub_pd(mur,_1)
                       ctht2  = _mm512_mul_pd(ctht,ctht);
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht   = xsin(tht);
                       k0r    = _mm512_mul_pd(k0,r);
                       stht2  = _mm512_mul_pd(stht,stht);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       num1r  = _mm512_mul_pd(murm1,ctht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       num2r  = _mm512_mul_pd(murm1,stht2);
                       den2r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       den3r  = _mm512_fmadd_pd(murm1,Ic,_2a2c);
                       num1i  = _mm512_mul_pd(muim1,ctht2);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       num2i  = _mm512_mul_pd(muim1,stht2);
                       den2i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       den3i  = _mm512_fmadd_pd(muim1,Ic,_2a2c);
                       cdiv_zmm8r8(epsrm1,epsim1,den1r,den1i,&div1r,&div1i);
                       cdiv_zmm8r8(num1r,num1i,den2r,den2i,&div2r,&div2i);
                       cdiv_zmm8r8(num2r,num2i,den3r,den3i,&div2r,&div2i);
                       t0r = _mm512_sub_pd(div1r,_mm512_sub_pd(div2r,div3r));
                       t0r = _mm512_sub_pd(div1i,_mm512_sub_pd(div2i,div3i));
                       cexp_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_storeu_pd(&ESr[0], resr);
                       _mm512_storeu_pd(&ESi[0], resi);
               }



                   /*
                          Low-frequency oblate spheroid Rayleigh
                          Backscattered fields.
                          E-field (theta), parallel, formula 5.1-88         
 
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5188_zmm8r8(const __m512d k0,
                                          const __m512d r,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          const __m512d a,
                                          const __m512d c,
                                          const __m512d tht,
                                          __m512d * __restrict ESr,
                                          __m512d * __restrict ESi) {

                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht,stht,facr,faci;
                       register __m512d ear,eai,cer,cei,epsrm1,epsim1,t2r,t2i,den3r,den3i;
                       register __m512d murm1,muim1,x0,t0r,t0i,_2a2c;
                       register __m512d den1r,den1i,den2r,den2i,div1r,div1i,div2r,div2i,div3r,div3i;
                       register __m512d num1r,num1i,num2r,num2i,ctht2,stht2;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht   = xcos(tht);
                       murm1  = _mm512_sub_pd(mur,_1)
                       ctht2  = _mm512_mul_pd(ctht,ctht);
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht   = xsin(tht);
                       k0r    = _mm512_mul_pd(k0,r);
                       stht2  = _mm512_mul_pd(stht,stht);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       num1r  = _mm512_mul_pd(epsrm1,ctht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       num2r  = _mm512_mul_pd(epsrm1,stht2);
                       den2r  = _mm512_fmadd_pd(epsrm1,Ic,_2a2c);
                       den3r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       num1i  = _mm512_mul_pd(epsim1,ctht2);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       num2i  = _mm512_mul_pd(epsim1,stht2);
                       den2i  = _mm512_fmadd_pd(epsim1,Ic,_2a2c);
                       den3i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       cdiv_zmm8r8(num1r,num1i,den1r,den1i,&div1r,&div1i);
                       cdiv_zmm8r8(num2r,num2i,den2r,den2i,&div2r,&div2i);
                       cdiv_zmm8r8(murm1,muim1,den3r,den3i,&div2r,&div2i);
                       t0r = _mm512_add_pd(div1r,_mm512_sub_pd(div2r,div3r));
                       t0r = _mm512_add_pd(div1i,_mm512_sub_pd(div2i,div3i));
                       cexp_zmm8r8(facr,faci,t0r,t0i,*ESr,*ESi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5188_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pc,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht,
                                          double * __restrict __ATTR_ALIGN__(64) ESr,
                                          double * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512d  k0    = _mm512_load_pd(&pk0[0]);
                       register __m512d  r     = _mm512_load_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_load_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_load_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_load_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_load_pd(&pmui[0]);
                       register __m512d  a     = _mm512_load_pd(&pa[0]);
                       register __m512d  c     = _mm512_load_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_load_pd(&ptht[0]); 
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht,stht,facr,faci;
                       register __m512d ear,eai,cer,cei,epsrm1,epsim1,t2r,t2i,den3r,den3i;
                       register __m512d murm1,muim1,x0,t0r,t0i,_2a2c,resr,resi;
                       register __m512d den1r,den1i,den2r,den2i,div1r,div1i,div2r,div2i,div3r,div3i;
                       register __m512d num1r,num1i,num2r,num2i,ctht2,stht2;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht   = xcos(tht);
                       murm1  = _mm512_sub_pd(mur,_1)
                       ctht2  = _mm512_mul_pd(ctht,ctht);
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht   = xsin(tht);
                       k0r    = _mm512_mul_pd(k0,r);
                       stht2  = _mm512_mul_pd(stht,stht);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       num1r  = _mm512_mul_pd(epsrm1,ctht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       num2r  = _mm512_mul_pd(epsrm1,stht2);
                       den2r  = _mm512_fmadd_pd(epsrm1,Ic,_2a2c);
                       den3r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       num1i  = _mm512_mul_pd(epsim1,ctht2);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       num2i  = _mm512_mul_pd(epsim1,stht2);
                       den2i  = _mm512_fmadd_pd(epsim1,Ic,_2a2c);
                       den3i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       cdiv_zmm8r8(num1r,num1i,den1r,den1i,&div1r,&div1i);
                       cdiv_zmm8r8(num2r,num2i,den2r,den2i,&div2r,&div2i);
                       cdiv_zmm8r8(murm1,muim1,den3r,den3i,&div2r,&div2i);
                       t0r = _mm512_add_pd(div1r,_mm512_sub_pd(div2r,div3r));
                       t0r = _mm512_add_pd(div1i,_mm512_sub_pd(div2i,div3i));
                       cexp_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_store_pd(&ESr[0], resr);
                       _mm512_store_pd(&ESi[0], resi);
               }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f5188_zmm8r8_u(const double * __restrict  pk0,
                                          const double * __restrict  pr,
                                          const double * __restrict  pepsr,
                                          const double * __restrict pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict pmui,
                                          const double * __restrict  pa,
                                          const double * __restrict  pc,
                                          const double * __restrict  ptht,
                                          double * __restrict  ESr,
                                          double * __restrict  ESi) {

                       register __m512d  k0    = _mm512_loadu_pd(&pk0[0]);
                       register __m512d  r     = _mm512_loadu_pd(&pr[0]);
                       register __m512d  epsr  = _mm512_loadu_pd(&pepsr[0]);
                       register __m512d  epsi  = _mm512_loadu_pd(&pepsi[0]);
                       register __m512d  mur   = _mm512_loadu_pd(&pmur[0]);
                       register __m512d  mui   = _mm512_loadu_pd(&pmui[0]);
                       register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                       register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                       register __m512d  tht1  = _mm512_loadu_pd(&ptht[0]); 
                       const __m512d _23 = _mm512_set1_pd(0.666666666666666666666666666667);
                       const __m512d _1  = _mm512_set1_pd(1.0f);
                       const __m512d _2  = _mm512_set1_pd(2.0f);
                       const __m512d _0  = _mm512_setzero_pd();
                       register __m512d k02,a2,invr,Ia,Ic,a2c,ctht,stht,facr,faci;
                       register __m512d ear,eai,cer,cei,epsrm1,epsim1,t2r,t2i,den3r,den3i;
                       register __m512d murm1,muim1,x0,t0r,t0i,_2a2c,resr,resi;
                       register __m512d den1r,den1i,den2r,den2i,div1r,div1i,div2r,div2i,div3r,div3i;
                       register __m512d num1r,num1i,num2r,num2i,ctht2,stht2;
                       epsrm1 = _mm512_sub_pd(epsr,_1);
                       k02    = _mm512_mul_pd(k0,k0);
                       ctht   = xcos(tht);
                       murm1  = _mm512_sub_pd(mur,_1)
                       ctht2  = _mm512_mul_pd(ctht,ctht);
                       invr   = _mm512_rcp14_pd(r);
                       epsim1 = _mm512_sub_pd(epsi,_1);
                       stht   = xsin(tht);
                       k0r    = _mm512_mul_pd(k0,r);
                       stht2  = _mm512_mul_pd(stht,stht);
                       muim1  = _mm512_sub_pd(mui,_1);
                       ear    = _0;
                       x0     = _mm512_mul_pd(_23,k02);
                       eai    = k0r;
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       Ia     = IaIb_f5191_zmm8r8(a,c);
                       facr   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cer));
                       a2c    = _mm512_mul_pd(c,_mm512_mul_pd(a,a));
                       faci   = _mm512_mul_pd(x0,_mm512_mul_pd(invr,cei));
                       _2a2c  = _mm512_div_pd(_2,a2c);
                       Ic     = Ic_f5192_zmm8r8(a,c);
                       num1r  = _mm512_mul_pd(epsrm1,ctht2);
                       den1r  = _mm512_fmadd_pd(epsrm1,Ia,_2a2c);
                       num2r  = _mm512_mul_pd(epsrm1,stht2);
                       den2r  = _mm512_fmadd_pd(epsrm1,Ic,_2a2c);
                       den3r  = _mm512_fmadd_pd(murm1,Ia,_2a2c);
                       num1i  = _mm512_mul_pd(epsim1,ctht2);
                       den1i  = _mm512_fmadd_pd(epsim1,Ia,_2a2c);
                       num2i  = _mm512_mul_pd(epsim1,stht2);
                       den2i  = _mm512_fmadd_pd(epsim1,Ic,_2a2c);
                       den3i  = _mm512_fmadd_pd(muim1,Ia,_2a2c);
                       cdiv_zmm8r8(num1r,num1i,den1r,den1i,&div1r,&div1i);
                       cdiv_zmm8r8(num2r,num2i,den2r,den2i,&div2r,&div2i);
                       cdiv_zmm8r8(murm1,muim1,den3r,den3i,&div2r,&div2i);
                       t0r = _mm512_add_pd(div1r,_mm512_sub_pd(div2r,div3r));
                       t0r = _mm512_add_pd(div1i,_mm512_sub_pd(div2i,div3i));
                       cexp_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                       _mm512_storeu_pd(&ESr[0], resr);
                       _mm512_storeu_pd(&ESi[0], resi);
               }


                    /*
                          High-frequency solutions.
                          Bistatic case RCS of oblate spheroid.
                          Formula 5.1-93
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5193_zmm8r8(const __m512d a,
                                            const __m512d c,
                                            const __m512d tht1,
                                            const __m512d tht2,
                                            const __m512d phi2) {

                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,num,a4,a2,c2,den,sden,trm1,trm2,trm3;
                          register __m512d stht1,stht2,cphi2,ctht1,ctht2,sphi2,strm1,strm2,strm3;
                          stht1 = xsin(tht1);
                          a2    = _mm512_mul_pd(a,a);
                          cphi2 = xcos(phi2);
                          c2    = _mm512_mul_pd(c,c);
                          ctht1 = xcos(tht1);
                          a4    = _mm512_mul_pd(a2,a2);
                          ctht2 = xcos(tht2);
                          num   = _mm512_mul_pd(_4pi,_mm512_mul_pd(a4,c2));
                          sphi2 = xsin(phi2);
                          stht2 = xsin(tht2);
                          trm1  = _mm512_fmadd_pd(stht2,cphi2,stht1);
                          strm1 = _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          trm2  = _mm512_mul_pd(stht2,sphi2);
                          strm2 = _mm512_mul_pd(a2,_mm512_mul_pd(trm2,trm2));
                          trm3  = _mm512_mul_pd(ctht1,ctht2);
                          strm3 = _mm512_mul_pd(c2,_mm512_mul_pd(trm3,trm3));
                          den   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          sden  = _mm512_mul_pd(den,den);
                          rcs   = _mm512_div_pd(num,sden);
                          return (rcs);
                } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5193_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) pc,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                             const double * __restrict __ATTR_ALIGN__(64) pphi2) {

                          register __m512d  a     = _mm512_load_pd(&pa[0]);
                          register __m512d  c     = _mm512_load_pd(&pc[0]);
                          register __m512d  tht1  = _mm512_load_pd(&ptht1[0]); 
                          register __m512d  tht2  = _mm512_load_pd(&ptht2[0]); 
                          register __m512d  pphi2 = _mm512_load_pd(&pphi2[0]); 
                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,num,a4,a2,c2,den,sden,trm1,trm2,trm3;
                          register __m512d stht1,stht2,cphi2,ctht1,ctht2,sphi2,strm1,strm2,strm3;
                          stht1 = xsin(tht1);
                          a2    = _mm512_mul_pd(a,a);
                          cphi2 = xcos(phi2);
                          c2    = _mm512_mul_pd(c,c);
                          ctht1 = xcos(tht1);
                          a4    = _mm512_mul_pd(a2,a2);
                          ctht2 = xcos(tht2);
                          num   = _mm512_mul_pd(_4pi,_mm512_mul_pd(a4,c2));
                          sphi2 = xsin(phi2);
                          stht2 = xsin(tht2);
                          trm1  = _mm512_fmadd_pd(stht2,cphi2,stht1);
                          strm1 = _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          trm2  = _mm512_mul_pd(stht2,sphi2);
                          strm2 = _mm512_mul_pd(a2,_mm512_mul_pd(trm2,trm2));
                          trm3  = _mm512_mul_pd(ctht1,ctht2);
                          strm3 = _mm512_mul_pd(c2,_mm512_mul_pd(trm3,trm3));
                          den   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          sden  = _mm512_mul_pd(den,den);
                          rcs   = _mm512_div_pd(num,sden);
                          return (rcs);
                } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5193_zmm8r8_u( const double * __restrict  pa,
                                             const double * __restrict  pc,
                                             const double * __restrict  ptht1,
                                             const double * __restrict  ptht2,
                                             const double * __restrict  pphi2) {

                          register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                          register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                          register __m512d  tht1  = _mm512_loadu_pd(&ptht1[0]); 
                          register __m512d  tht2  = _mm512_loadu_pd(&ptht2[0]); 
                          register __m512d  pphi2 = _mm512_loadu_pd(&pphi2[0]); 
                          const __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,num,a4,a2,c2,den,sden,trm1,trm2,trm3;
                          register __m512d stht1,stht2,cphi2,ctht1,ctht2,sphi2,strm1,strm2,strm3;
                          stht1 = xsin(tht1);
                          a2    = _mm512_mul_pd(a,a);
                          cphi2 = xcos(phi2);
                          c2    = _mm512_mul_pd(c,c);
                          ctht1 = xcos(tht1);
                          a4    = _mm512_mul_pd(a2,a2);
                          ctht2 = xcos(tht2);
                          num   = _mm512_mul_pd(_4pi,_mm512_mul_pd(a4,c2));
                          sphi2 = xsin(phi2);
                          stht2 = xsin(tht2);
                          trm1  = _mm512_fmadd_pd(stht2,cphi2,stht1);
                          strm1 = _mm512_mul_pd(a2,_mm512_mul_pd(trm1,trm1));
                          trm2  = _mm512_mul_pd(stht2,sphi2);
                          strm2 = _mm512_mul_pd(a2,_mm512_mul_pd(trm2,trm2));
                          trm3  = _mm512_mul_pd(ctht1,ctht2);
                          strm3 = _mm512_mul_pd(c2,_mm512_mul_pd(trm3,trm3));
                          den   = _mm512_add_pd(strm1,_mm512_add_pd(strm2,strm3));
                          sden  = _mm512_mul_pd(den,den);
                          rcs   = _mm512_div_pd(num,sden);
                          return (rcs);
                } 


                   /*
                          High-frequency solutions.
                          Backscatter case RCS of oblate spheroid.
                          Formula 5.1-94
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5194_zmm8r8(const __m512d a,
                                            const __m512d c,
                                            const __m512d tht) {

                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a4,a2,c2,num,den,sden;
                          register __m512d stht,ctht,stht2,ctht2;
                          stht = xsin(tht);
                          a2   = _mm512_mul_pd(a,a);
                          ctht = xcos(tht);
                          c2   = _mm512_mul_pd(c,c);
                          stht2= _mm512_mul_pd(stht,stht);
                          a4   = _mm512_mul_pd(a2,a2);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          num  = _mm512_mul_pd(pi,_mm512_mul_pd(a4,c2));
                          den  = _mm512_fmadd_pd(a2,stht2,_mm512_mul_pd(c2,ctht2));
                          sden = _mm512_mul_pd(den,den);
                          rcs  = _mm512_div_pd(num,sden);
                          return (rcs):
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5194_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) pc,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d  a     = _mm512_load_pd(&pa[0]);
                          register __m512d  c     = _mm512_load_pd(&pc[0]);
                          register __m512d  tht  = _mm512_load_pd(&ptht[0]); 
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a4,a2,c2,num,den,sden;
                          register __m512d stht,ctht,stht2,ctht2;
                          stht = xsin(tht);
                          a2   = _mm512_mul_pd(a,a);
                          ctht = xcos(tht);
                          c2   = _mm512_mul_pd(c,c);
                          stht2= _mm512_mul_pd(stht,stht);
                          a4   = _mm512_mul_pd(a2,a2);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          num  = _mm512_mul_pd(pi,_mm512_mul_pd(a4,c2));
                          den  = _mm512_fmadd_pd(a2,stht2,_mm512_mul_pd(c2,ctht2));
                          sden = _mm512_mul_pd(den,den);
                          rcs  = _mm512_div_pd(num,sden);
                          return (rcs):
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5194_zmm8r8_u(const double * __restrict  pa,
                                             const double * __restrict  pc,
                                             const double * __restrict  ptht) {

                          register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                          register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                          register __m512d  tht  = _mm512_loadu_pd(&ptht[0]); 
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          register __m512d rcs,a4,a2,c2,num,den,sden;
                          register __m512d stht,ctht,stht2,ctht2;
                          stht = xsin(tht);
                          a2   = _mm512_mul_pd(a,a);
                          ctht = xcos(tht);
                          c2   = _mm512_mul_pd(c,c);
                          stht2= _mm512_mul_pd(stht,stht);
                          a4   = _mm512_mul_pd(a2,a2);
                          ctht2= _mm512_mul_pd(ctht,ctht);
                          num  = _mm512_mul_pd(pi,_mm512_mul_pd(a4,c2));
                          den  = _mm512_fmadd_pd(a2,stht2,_mm512_mul_pd(c2,ctht2));
                          sden = _mm512_mul_pd(den,den);
                          rcs  = _mm512_div_pd(num,sden);
                          return (rcs):
                 }


                    /*
                           Prolate spheroid.
                           Axial incidence "Physical Optics" backscatter RCS.
                           Formula 5.1-95
                      */

 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5195_zmm8r8(const __m512d a,
                                            const __m512d c,
                                            const __m512d k0) {

                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,k0c,a4,a2,c2,_2k0c,k0cs;
                          register __m512d sarg,trm1,trm2,trm3,x0;
                          k0c  = _mm512_mul_pd(k0,c);
                          a2   = _mm512_mul_pd(a,a);
                          _2k0c= _mm512_add_pd(k0c,k0c);
                          c2   = _mm512_mul_pd(c,c);
                          sarg = xsin(_2k0c);
                          trm2 = _mm512_div_pd(sarg,k0c);
                          a4   = _mm512_mul_pd(a2,a2);
                          k0cs = _mm512_mul_pd(k0c,k0c);
                          trm1 = _mm512_div_pd(_mm512_mul_pd(pi,ar),c2);
                          trm3 = _mm512_div_pd(sarg,k0cs);
                          x0   = _mm512_add_pd(_mm512_sub_pd(_1,trm2),trm3);
                          rcs  = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5195_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pc,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512d  a     = _mm512_load_pd(&pa[0]);
                          register __m512d  c     = _mm512_load_pd(&pc[0]);
                          register __m512d  k0    = _mm512_load_pd(&pk0[0]); 
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,k0c,a4,a2,c2,_2k0c,k0cs;
                          register __m512d sarg,trm1,trm2,trm3,x0;
                          k0c  = _mm512_mul_pd(k0,c);
                          a2   = _mm512_mul_pd(a,a);
                          _2k0c= _mm512_add_pd(k0c,k0c);
                          c2   = _mm512_mul_pd(c,c);
                          sarg = xsin(_2k0c);
                          trm2 = _mm512_div_pd(sarg,k0c);
                          a4   = _mm512_mul_pd(a2,a2);
                          k0cs = _mm512_mul_pd(k0c,k0c);
                          trm1 = _mm512_div_pd(_mm512_mul_pd(pi,ar),c2);
                          trm3 = _mm512_div_pd(sarg,k0cs);
                          x0   = _mm512_add_pd(_mm512_sub_pd(_1,trm2),trm3);
                          rcs  = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5195_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pc,
                                              const double * __restrict  pk0) {

                          register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                          register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                          register __m512d  k0    = _mm512_loadu_pd(&pk0[0]); 
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,k0c,a4,a2,c2,_2k0c,k0cs;
                          register __m512d sarg,trm1,trm2,trm3,x0;
                          k0c  = _mm512_mul_pd(k0,c);
                          a2   = _mm512_mul_pd(a,a);
                          _2k0c= _mm512_add_pd(k0c,k0c);
                          c2   = _mm512_mul_pd(c,c);
                          sarg = xsin(_2k0c);
                          trm2 = _mm512_div_pd(sarg,k0c);
                          a4   = _mm512_mul_pd(a2,a2);
                          k0cs = _mm512_mul_pd(k0c,k0c);
                          trm1 = _mm512_div_pd(_mm512_mul_pd(pi,ar),c2);
                          trm3 = _mm512_div_pd(sarg,k0cs);
                          x0   = _mm512_add_pd(_mm512_sub_pd(_1,trm2),trm3);
                          rcs  = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                 }


                  /*
                           Prolate spheroid.
                           Incidence perpendicular to the axis. 
                           "Physical Optics" backscatter RCS.
                           Formula 5.1-96
                      */


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5196_zmm8r8(const __m512d a,
                                            const __m512d c,
                                            const __m512d k0) {

                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,k0a,c2,_2k0a,k0as;
                          register __m512d sarg,trm1,trm2,trm3,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          c2   = _mm512_mul_pd(c,c);
                          _2k0a= _mm512_add_pd(k0a,k0a);
                          sarg = xsin(_2k0a);
                          trm2 = _mm512_div_pd(sarg,k0a);
                          k0as = _mm512_mul_pd(k0a,k0a);
                          trm1 = _mm512_mul_pd(pi,c2);
                          x1   = _mm512_mul_pd(sarg,sarg);
                          trm3 = _mm512_div_pd(x1,k0as);
                          x0   = _mm512_add_pd(_mm512_sub_pd(_1,trm2),trm3);
                          rcs  = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5196_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pc,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512d  a     = _mm512_load_pd(&pa[0]);
                          register __m512d  c     = _mm512_load_pd(&pc[0]);
                          register __m512d  k0    = _mm512_load_pd(&pk0[0]); 
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,k0a,c2,_2k0a,k0as;
                          register __m512d sarg,trm1,trm2,trm3,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          c2   = _mm512_mul_pd(c,c);
                          _2k0a= _mm512_add_pd(k0a,k0a);
                          sarg = xsin(_2k0a);
                          trm2 = _mm512_div_pd(sarg,k0a);
                          k0as = _mm512_mul_pd(k0a,k0a);
                          trm1 = _mm512_mul_pd(pi,c2);
                          x1   = _mm512_mul_pd(sarg,sarg);
                          trm3 = _mm512_div_pd(x1,k0as);
                          x0   = _mm512_add_pd(_mm512_sub_pd(_1,trm2),trm3);
                          rcs  = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                 }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5196_zmm8r8_u(  const double * __restrict  pa,
                                              const double * __restrict  pc,
                                              const double * __restrict  pk0) {

                          register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                          register __m512d  c     = _mm512_loadu_pd(&pc[0]);
                          register __m512d  k0    = _mm512_loadu_pd(&pk0[0]); 
                          const __m512d pi = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d rcs,k0a,c2,_2k0a,k0as;
                          register __m512d sarg,trm1,trm2,trm3,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          c2   = _mm512_mul_pd(c,c);
                          _2k0a= _mm512_add_pd(k0a,k0a);
                          sarg = xsin(_2k0a);
                          trm2 = _mm512_div_pd(sarg,k0a);
                          k0as = _mm512_mul_pd(k0a,k0a);
                          trm1 = _mm512_mul_pd(pi,c2);
                          x1   = _mm512_mul_pd(sarg,sarg);
                          trm3 = _mm512_div_pd(x1,k0as);
                          x0   = _mm512_add_pd(_mm512_sub_pd(_1,trm2),trm3);
                          rcs  = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                 }


                   /*
                        Perfectly conducting ogive.
                        Backscatter RCS, axial incidence.
                        Parabolic ogive.
                        Formula 5.2-1
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f521_zmm8r8(const __m512d a,
                                           const __m512d b,
                                           const __m512d k0) {

                          const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _4pi = _mm512_set1_pd(1.27323954473516268615107010698);
                          const __m512d _1   = _mm512_set1_pd(1.0f);
                          const __m512d c0   = _mm512_set1_pd(0.8f);
                          const __m512d c1   = _mm512_set1_pd(1.066666666666666666666666666667);
                          register __m512d rcs,V,rho,V2,k04,k02,expr;
                          register __m512d trm1,trm2,b2,invr,x0,inve,strm2;
                          k02  = _mm512_mul_pd(k0,k0);
                          rho  = _mm512_mul_pd(c0,_mm512_div_pd(a,b));
                          k04  = _mm512_mul_pd(k02,k02);
                          b2   = _mm512_mul_pd(b,b);
                          x0   = _mm512_mul_pd(pi,rho);
                          expr = xexpf(rho);
                          invr = _mm512_rcp14_pd(x0);
                          V    = _mm512_mul_pd(_mm512_set1_pd(3.351032163829112787693486275498),
                                               _mm512_mul_pd(a,b2));
                          inve = _mm512_rcp14_pd(expr);
                          V2   = _mm512_mul_pd(V,V);
                          trm2 = _mm512_fmadd_pd(invr,inve,_1);
                          trm1 = _mm512_mul_pd(_4pi,_mm512_mul_pd(k04,V2));
                          strm2= _mm512_mul_pd(trm2,trm2);
                          rcs  = _mm512_mul_pd(trm1,strm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f521_zmm8r8_a(   const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512d  a     = _mm512_load_pd(&pa[0]);
                          register __m512d  b     = _mm512_load_pd(&pb[0]);
                          register __m512d  k0    = _mm512_load_pd(&pk0[0]); 
                          const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _4pi = _mm512_set1_pd(1.27323954473516268615107010698);
                          const __m512d _1   = _mm512_set1_pd(1.0f);
                          const __m512d c0   = _mm512_set1_pd(0.8f);
                          const __m512d c1   = _mm512_set1_pd(1.066666666666666666666666666667);
                          register __m512d rcs,V,rho,V2,k04,k02,expr;
                          register __m512d trm1,trm2,b2,invr,x0,inve,strm2;
                          k02  = _mm512_mul_pd(k0,k0);
                          rho  = _mm512_mul_pd(c0,_mm512_div_pd(a,b));
                          k04  = _mm512_mul_pd(k02,k02);
                          b2   = _mm512_mul_pd(b,b);
                          x0   = _mm512_mul_pd(pi,rho);
                          expr = xexpf(rho);
                          invr = _mm512_rcp14_pd(x0);
                          V    = _mm512_mul_pd(_mm512_set1_pd(3.351032163829112787693486275498),
                                               _mm512_mul_pd(a,b2));
                          inve = _mm512_rcp14_pd(expr);
                          V2   = _mm512_mul_pd(V,V);
                          trm2 = _mm512_fmadd_pd(invr,inve,_1);
                          trm1 = _mm512_mul_pd(_4pi,_mm512_mul_pd(k04,V2));
                          strm2= _mm512_mul_pd(trm2,trm2);
                          rcs  = _mm512_mul_pd(trm1,strm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f521_zmm8r8_u(   const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pk0) {

                          register __m512d  a     = _mm512_loadu_pd(&pa[0]);
                          register __m512d  b     = _mm512_loadu_pd(&pb[0]);
                          register __m512d  k0    = _mm512_loadu_pd(&pk0[0]); 
                          const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _4pi = _mm512_set1_pd(1.27323954473516268615107010698);
                          const __m512d _1   = _mm512_set1_pd(1.0f);
                          const __m512d c0   = _mm512_set1_pd(0.8f);
                          const __m512d c1   = _mm512_set1_pd(1.066666666666666666666666666667);
                          register __m512d rcs,V,rho,V2,k04,k02,expr;
                          register __m512d trm1,trm2,b2,invr,x0,inve,strm2;
                          k02  = _mm512_mul_pd(k0,k0);
                          rho  = _mm512_mul_pd(c0,_mm512_div_pd(a,b));
                          k04  = _mm512_mul_pd(k02,k02);
                          b2   = _mm512_mul_pd(b,b);
                          x0   = _mm512_mul_pd(pi,rho);
                          expr = xexpf(rho);
                          invr = _mm512_rcp14_pd(x0);
                          V    = _mm512_mul_pd(_mm512_set1_pd(3.351032163829112787693486275498),
                                               _mm512_mul_pd(a,b2));
                          inve = _mm512_rcp14_pd(expr);
                          V2   = _mm512_mul_pd(V,V);
                          trm2 = _mm512_fmadd_pd(invr,inve,_1);
                          trm1 = _mm512_mul_pd(_4pi,_mm512_mul_pd(k04,V2));
                          strm2= _mm512_mul_pd(trm2,trm2);
                          rcs  = _mm512_mul_pd(trm1,strm2);
                          return (rcs);
                }


                  /*
                        Perfectly conducting ogive.
                        Backscatter RCS, axial incidence.
                        Circular ogive.
                        Formula 5.2-4
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f524_zmm8r8(const __m512d alp,
                                           const __m512d r0,
                                           const __m512d k0) {

                          const __m512d tpi   = _mm512_set1_pd(2.0f*3.14159265358979323846264338328);
                          const __m512d _4pi  = _mm512_set1_pd(1.27323954473516268615107010698);
                          const __m512d c0    = _mm512_set1_pd(0.33333333333333333333333333333333333);
                          const __m512d c1    = _mm512_set1_pd(1.5f);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d rcs,rho,V,V2,k04,salp,calp,trm2,trm3;
                          register __m512d num,den,expr,invr,inve,x0,x1,x2,trm1,x3;
                          x0  = _mm512_mul_pd(k0,k0);
                          x1  = _mm512_mul_pd(tpi,_mm512_mul_pd(r0,r0));
                          salp= xsin(alp);
                          k04 = _mm512_mul_pd(x0,x0);
                          calp= xcos(alp);
                          trm1= _mm512_mul_pd(_4pi,k04);
                          x0  = _mm512_mul_pd(salp,_mm512_mul_pd(salp,salp));
                          x3  = _mm512_sub_pd(_1,calp);
                          x2  = _mm512_mul_pd(c0,x0);
                          den = _mm512_mul_pd(x3,_mm512_mul_pd(x3,x3));
                          num = _mm512_sub_pd(salp,_mm512_fmsub_pd(alp,calp,x2));
                          rho = _mm512_mul_pd(c1,_mm512_div_pd(num,den));
                          invr= _mm512_rcp14_pd(rho);
                          V   = _mm512_mul_pd(x1,num);
                          expr= xexpf(rho);
                          V2  = _mm512_mul_pd(V,V);
                          inve= _mm512_rcp14_pd(expr);
                          trm2= _mm512_mul_pd(trm1,V2);
                          trm3= _mm512_fmadd_pd(inve,invr,_1);
                          rcs = _mm512_mul_pd(trm2,trm3);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f524_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) palp,
                                              const double * __restrict __ATTR_ALIGN__(64) pr0,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512d  alp   = _mm512_load_pd(&palp[0]);
                          register __m512d  r0    = _mm512_load_pd(&pr0[0]);
                          register __m512d  k0    = _mm512_load_pd(&pk0[0]); 
                          const __m512d tpi   = _mm512_set1_pd(2.0f*3.14159265358979323846264338328);
                          const __m512d _4pi  = _mm512_set1_pd(1.27323954473516268615107010698);
                          const __m512d c0    = _mm512_set1_pd(0.33333333333333333333333333333333333);
                          const __m512d c1    = _mm512_set1_pd(1.5f);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d rcs,rho,V,V2,k04,salp,calp,trm2,trm3;
                          register __m512d num,den,expr,invr,inve,x0,x1,x2,trm1,x3;
                          x0  = _mm512_mul_pd(k0,k0);
                          x1  = _mm512_mul_pd(tpi,_mm512_mul_pd(r0,r0));
                          salp= xsin(alp);
                          k04 = _mm512_mul_pd(x0,x0);
                          calp= xcos(alp);
                          trm1= _mm512_mul_pd(_4pi,k04);
                          x0  = _mm512_mul_pd(salp,_mm512_mul_pd(salp,salp));
                          x3  = _mm512_sub_pd(_1,calp);
                          x2  = _mm512_mul_pd(c0,x0);
                          den = _mm512_mul_pd(x3,_mm512_mul_pd(x3,x3));
                          num = _mm512_sub_pd(salp,_mm512_fmsub_pd(alp,calp,x2));
                          rho = _mm512_mul_pd(c1,_mm512_div_pd(num,den));
                          invr= _mm512_rcp14_pd(rho);
                          V   = _mm512_mul_pd(x1,num);
                          expr= xexpf(rho);
                          V2  = _mm512_mul_pd(V,V);
                          inve= _mm512_rcp14_pd(expr);
                          trm2= _mm512_mul_pd(trm1,V2);
                          trm3= _mm512_fmadd_pd(inve,invr,_1);
                          rcs = _mm512_mul_pd(trm2,trm3);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f524_zmm8r8_u( const double * __restrict  palp,
                                              const double * __restrict  pr0,
                                              const double * __restrict  pk0) {

                          register __m512d  alp   = _mm512_loadu_pd(&palp[0]);
                          register __m512d  r0    = _mm512_loadu_pd(&pr0[0]);
                          register __m512d  k0    = _mm512_loadu_pd(&pk0[0]); 
                          const __m512d tpi   = _mm512_set1_pd(2.0f*3.14159265358979323846264338328);
                          const __m512d _4pi  = _mm512_set1_pd(1.27323954473516268615107010698);
                          const __m512d c0    = _mm512_set1_pd(0.33333333333333333333333333333333333);
                          const __m512d c1    = _mm512_set1_pd(1.5f);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d rcs,rho,V,V2,k04,salp,calp,trm2,trm3;
                          register __m512d num,den,expr,invr,inve,x0,x1,x2,trm1,x3;
                          x0  = _mm512_mul_pd(k0,k0);
                          x1  = _mm512_mul_pd(tpi,_mm512_mul_pd(r0,r0));
                          salp= xsin(alp);
                          k04 = _mm512_mul_pd(x0,x0);
                          calp= xcos(alp);
                          trm1= _mm512_mul_pd(_4pi,k04);
                          x0  = _mm512_mul_pd(salp,_mm512_mul_pd(salp,salp));
                          x3  = _mm512_sub_pd(_1,calp);
                          x2  = _mm512_mul_pd(c0,x0);
                          den = _mm512_mul_pd(x3,_mm512_mul_pd(x3,x3));
                          num = _mm512_sub_pd(salp,_mm512_fmsub_pd(alp,calp,x2));
                          rho = _mm512_mul_pd(c1,_mm512_div_pd(num,den));
                          invr= _mm512_rcp14_pd(rho);
                          V   = _mm512_mul_pd(x1,num);
                          expr= xexpf(rho);
                          V2  = _mm512_mul_pd(V,V);
                          inve= _mm512_rcp14_pd(expr);
                          trm2= _mm512_mul_pd(trm1,V2);
                          trm3= _mm512_fmadd_pd(inve,invr,_1);
                          rcs = _mm512_mul_pd(trm2,trm3);
                          return (rcs);
                }


                 /*
                        Circular ogive.
                        Axial incidence -- RCS.
                        Formula 5.2-9
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f529_zmm8r8(const __m512d gam0,
                                           const __m512d alp,
                                           const __m512d k0a) {

                          const __m512d _16pi = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d rcs,trm1,trm2,trm3,talp,calp,_2k0a,tan2,x0;
                          register __m512d gam2,tan4,calp4,den1,num,carg,den2;
                          talp   = xtan(alp);
                          gam2   = _mm512_mul_pd(gam0,gam0);
                          calp   = xcos(alp);
                          _2k0a  = _mm512_add_pd(k0a,k0a);
                          den1   = _mm512_add_pd(_1,calp);
                          tan2   = _mm512_mul_pd(talp,talp);
                          den2   = _mm512_mul_pd(den1,den1);
                          tan4   = _mm512_mul_pd(tan2,tan2);
                          x0     = _mm512_mul_pd(calp,calp);
                          trm1   = _mm512_div_pd(_mm512_mul_pd(gam2,tan4),_16pi);
                          calp4  = _mm512_mul_pd(x0,x0);
                          carg   = xcos(_2k0a);
                          trm3   = _mm512_div_pd(calp4,den2);
                          num    = _mm512_mul_pd(x0,carg);
                          trm2   = _mm512_div_pd(num,den1);
                          x0     = _mm512_add_pd(_1,_mm512_add_pd(trm2,trm3));
                          rcs    = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f529_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const double * __restrict __ATTR_ALIGN__(64) palp,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d  gam0   = _mm512_load_pd(&pgam0[0]);
                          register __m512d  alp    = _mm512_load_pd(&palp[0]);
                          register __m512d  k0a    = _mm512_load_pd(&pk0a[0]); 
                          const __m512d _16pi = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d rcs,trm1,trm2,trm3,talp,calp,_2k0a,tan2,x0;
                          register __m512d gam2,tan4,calp4,den1,num,carg,den2;
                          talp   = xtan(alp);
                          gam2   = _mm512_mul_pd(gam0,gam0);
                          calp   = xcos(alp);
                          _2k0a  = _mm512_add_pd(k0a,k0a);
                          den1   = _mm512_add_pd(_1,calp);
                          tan2   = _mm512_mul_pd(talp,talp);
                          den2   = _mm512_mul_pd(den1,den1);
                          tan4   = _mm512_mul_pd(tan2,tan2);
                          x0     = _mm512_mul_pd(calp,calp);
                          trm1   = _mm512_div_pd(_mm512_mul_pd(gam2,tan4),_16pi);
                          calp4  = _mm512_mul_pd(x0,x0);
                          carg   = xcos(_2k0a);
                          trm3   = _mm512_div_pd(calp4,den2);
                          num    = _mm512_mul_pd(x0,carg);
                          trm2   = _mm512_div_pd(num,den1);
                          x0     = _mm512_add_pd(_1,_mm512_add_pd(trm2,trm3));
                          rcs    = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f529_zmm8r8_u( const double * __restrict pgam0,
                                              const double * __restrict  palp,
                                              const double * __restrict  pk0a) {

                          register __m512d  gam0   = _mm512_loadu_pd(&pgam0[0]);
                          register __m512d  alp    = _mm512_loadu_pd(&palp[0]);
                          register __m512d  k0a    = _mm512_loadu_pd(&pk0a[0]); 
                          const __m512d _16pi = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d rcs,trm1,trm2,trm3,talp,calp,_2k0a,tan2,x0;
                          register __m512d gam2,tan4,calp4,den1,num,carg,den2;
                          talp   = xtan(alp);
                          gam2   = _mm512_mul_pd(gam0,gam0);
                          calp   = xcos(alp);
                          _2k0a  = _mm512_add_pd(k0a,k0a);
                          den1   = _mm512_add_pd(_1,calp);
                          tan2   = _mm512_mul_pd(talp,talp);
                          den2   = _mm512_mul_pd(den1,den1);
                          tan4   = _mm512_mul_pd(tan2,tan2);
                          x0     = _mm512_mul_pd(calp,calp);
                          trm1   = _mm512_div_pd(_mm512_mul_pd(gam2,tan4),_16pi);
                          calp4  = _mm512_mul_pd(x0,x0);
                          carg   = xcos(_2k0a);
                          trm3   = _mm512_div_pd(calp4,den2);
                          num    = _mm512_mul_pd(x0,carg);
                          trm2   = _mm512_div_pd(num,den1);
                          x0     = _mm512_add_pd(_1,_mm512_add_pd(trm2,trm3));
                          rcs    = _mm512_mul_pd(trm1,x0);
                          return (rcs);
                }


                   /*
                           Circular ogive.
                           RCS as function of theta angle, i.e. (0<<theta<<(90-alpha))
                           Formula 5.2-6
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f526_zmm8r8(const __m512d gam0,
                                           const __m512d alp,
                                           const __m512d tht) { // in degrees.
                                           
                          const __m512d _16pi  = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _1     = _mm512_set1_pd(1.0f);
                          register __m512d rcs,talp,talp2,talp4,tath,tath2;
                          register __m512d num,den,ctht,ctht2,ctht6,_1mt2;
                          register __m512d x0,x1,gam2;
                          gam2  = _mm512_mul_pd(gam0,gam0);
                          talp  = xtan(alp);
                          ctht  = xcos(tht);
                          tath  = xtan(tht);
                          talp2 = _mm512_mul_pd(talp,talp);
                          ctht2 = _mm512_mul_pd(ctht,ctht);
                          _1mt2 = _mm512_sub_pd(_1,talp2);
                          ctht6 = _mm512_mul_pd(ctht2,_mm512_mul_pd(ctht2,ctht2));
                          tath2 = _mm512_mul_pd(tath,tath);
                          x1    = _mm512_mul_pd(_1mt2,tath2);
                          num   = _mm512_mul_pd(gam2,_mm512_mul_pd(talp2,talp2));
                          x0    = _mm512_mul_pd(_16pi,ctht6);
                          den   = _mm512_mul_pd(x1,_mm512_mul_pd(x1,x1));
                          den   = _mm512_mul_pd(x0,den);
                          rcs   = _mm512_div_pd(num,den);
                          return (rcs);
                }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f526_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const double * __restrict __ATTR_ALIGN__(64) palp,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht) { // in degrees.
                           
                          register __m512d  gam0   = _mm512_load_pd(&pgam0[0]);
                          register __m512d  alp    = _mm512_load_pd(&palp[0]);
                          register __m512d  tht    = _mm512_load_pd(&ptht[0]);                 
                          const __m512d _16pi  = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _1     = _mm512_set1_pd(1.0f);
                          register __m512d rcs,talp,talp2,talp4,tath,tath2;
                          register __m512d num,den,ctht,ctht2,ctht6,_1mt2;
                          register __m512d x0,x1,gam2;
                          gam2  = _mm512_mul_pd(gam0,gam0);
                          talp  = xtan(alp);
                          ctht  = xcos(tht);
                          tath  = xtan(tht);
                          talp2 = _mm512_mul_pd(talp,talp);
                          ctht2 = _mm512_mul_pd(ctht,ctht);
                          _1mt2 = _mm512_sub_pd(_1,talp2);
                          ctht6 = _mm512_mul_pd(ctht2,_mm512_mul_pd(ctht2,ctht2));
                          tath2 = _mm512_mul_pd(tath,tath);
                          x1    = _mm512_mul_pd(_1mt2,tath2);
                          num   = _mm512_mul_pd(gam2,_mm512_mul_pd(talp2,talp2));
                          x0    = _mm512_mul_pd(_16pi,ctht6);
                          den   = _mm512_mul_pd(x1,_mm512_mul_pd(x1,x1));
                          den   = _mm512_mul_pd(x0,den);
                          rcs   = _mm512_div_pd(num,den);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f526_zmm8r8_u(const double * __restrict  pgam0,
                                              const double * __restrict  palp,
                                              const double * __restrict  ptht) { // in degrees.
                           
                          register __m512d  gam0   = _mm512_loadu_pd(&pgam0[0]);
                          register __m512d  alp    = _mm512_loadu_pd(&palp[0]);
                          register __m512d  tht    = _mm512_loadu_pd(&ptht[0]);                 
                          const __m512d _16pi  = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _1     = _mm512_set1_pd(1.0f);
                          register __m512d rcs,talp,talp2,talp4,tath,tath2;
                          register __m512d num,den,ctht,ctht2,ctht6,_1mt2;
                          register __m512d x0,x1,gam2;
                          gam2  = _mm512_mul_pd(gam0,gam0);
                          talp  = xtan(alp);
                          ctht  = xcos(tht);
                          tath  = xtan(tht);
                          talp2 = _mm512_mul_pd(talp,talp);
                          ctht2 = _mm512_mul_pd(ctht,ctht);
                          _1mt2 = _mm512_sub_pd(_1,talp2);
                          ctht6 = _mm512_mul_pd(ctht2,_mm512_mul_pd(ctht2,ctht2));
                          tath2 = _mm512_mul_pd(tath,tath);
                          x1    = _mm512_mul_pd(_1mt2,tath2);
                          num   = _mm512_mul_pd(gam2,_mm512_mul_pd(talp2,talp2));
                          x0    = _mm512_mul_pd(_16pi,ctht6);
                          den   = _mm512_mul_pd(x1,_mm512_mul_pd(x1,x1));
                          den   = _mm512_mul_pd(x0,den);
                          rcs   = _mm512_div_pd(num,den);
                          return (rcs);
                }


                    /*
                           Circular ogive.
                           RCS as function of theta angle, i.e. theta = 90-alpha
                           Formula 5.2-7
                     */


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f527_zmm8r8(const __m512d b,
                                           const __m512d alp) {

                          const __m512d _4pi  = _mm512_set1_pd(12.566370614359172953850573533118);
                          const __m512d hlf   = _mm512_set1_pd(0.5f);
                          register __m512d rcs,b2,alp2,talp,x0;
                          b2  = _mm512_mul_pd(b,b);
                          alp2= _mm512_mul_pd(alp,hlf);
                          talp= xtan(alp2);
                          x0  = _mm512_mul_pd(_4pi,_mm512_mul_pd(talp,talp));
                          rcs = _mm512_div_pd(b2,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f527_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) palp) {

                          register __m512d  b    = _mm512_load_pd(&pb[0]);
                          register __m512d  alp  = _mm512_load_pd(&palp[0]);    
                          const __m512d _4pi  = _mm512_set1_pd(12.566370614359172953850573533118);
                          const __m512d hlf   = _mm512_set1_pd(0.5f);
                          register __m512d rcs,b2,alp2,talp,x0;
                          b2  = _mm512_mul_pd(b,b);
                          alp2= _mm512_mul_pd(alp,hlf);
                          talp= xtan(alp2);
                          x0  = _mm512_mul_pd(_4pi,_mm512_mul_pd(talp,talp));
                          rcs = _mm512_div_pd(b2,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f527_zmm8r8_u(const double * __restrict  pb,
                                              const double * __restrict  palp) {

                          register __m512d  b    = _mm512_loadu_pd(&pb[0]);
                          register __m512d  alp  = _mm512_loadu_pd(&palp[0]);    
                          const __m512d _4pi  = _mm512_set1_pd(12.566370614359172953850573533118);
                          const __m512d hlf   = _mm512_set1_pd(0.5f);
                          register __m512d rcs,b2,alp2,talp,x0;
                          b2  = _mm512_mul_pd(b,b);
                          alp2= _mm512_mul_pd(alp,hlf);
                          talp= xtan(alp2);
                          x0  = _mm512_mul_pd(_4pi,_mm512_mul_pd(talp,talp));
                          rcs = _mm512_div_pd(b2,x0);
                          return (rcs);
                 }


                   /*
                           Circular ogive.
                           RCS as function of theta angle, i.e. theta = (90-alpha) < theta << 90
                           Formula 5.2-8
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f528_zmm8r8(const __m512d r0,
                                           const __m512d b,
                                           const __m512d tht) {

                          const __m512d pi  = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          register __m512d rcs,num,den,r02,stht,x0,rat;
                          r02  = _mm512_mul_pd(r0,r0);
                          stht = xsin(tht);
                          x0   = _mm512_mul_pd(pi,r02);
                          num  = _mm512_sub_pd(r0,b);
                          den  = _mm512_mul_pd(r0,stht);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512_mul_pd(x0,_mm512_sub_pd(_1,rat));
                          return (rat);
                 }

                    
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f528_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pb,
                                           const double * __restrict __ATTR_ALIGN__(64) pr0,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht ) {

                          register __m512d  b     = _mm512_load_pd(&pb[0]);
                          register __m512d  r0    = _mm512_load_pd(&pr0[0]);
                          register __m512d  tht    = _mm512_load_pd(&ptht[0]);    
                          const __m512d pi  = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          register __m512d rcs,num,den,r02,stht,x0,rat;
                          r02  = _mm512_mul_pd(r0,r0);
                          stht = xsin(tht);
                          x0   = _mm512_mul_pd(pi,r02);
                          num  = _mm512_sub_pd(r0,b);
                          den  = _mm512_mul_pd(r0,stht);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512_mul_pd(x0,_mm512_sub_pd(_1,rat));
                          return (rat);
                 }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f528_zmm8r8_u(const double * __restrict  pb,
                                           const double * __restrict  pr0,
                                           const double * __restrict  ptht ) {

                          register __m512d  b     = _mm512_loadu_pd(&pb[0]);
                          register __m512d  r0    = _mm512_loadu_pd(&pr0[0]);
                          register __m512d  tht   = _mm512_loadu_pd(&ptht[0]);    
                          const __m512d pi  = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          register __m512d rcs,num,den,r02,stht,x0,rat;
                          r02  = _mm512_mul_pd(r0,r0);
                          stht = xsin(tht);
                          x0   = _mm512_mul_pd(pi,r02);
                          num  = _mm512_sub_pd(r0,b);
                          den  = _mm512_mul_pd(r0,stht);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512_mul_pd(x0,_mm512_sub_pd(_1,rat));
                          return (rat);
                 }


                  /*
                       Dispatch kernel for Circular ogive RCS.
                       Formulae: 5.2-6, 5.2-7, 5.2-8
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f52678_zmm8r8(const __m512d gam0,
                                             const __m512d alp,
                                             const __m512d tht,
                                             const __m512d b,
                                             const __m512d r0) {

                          const __m512d _90        = _mm512_set1_pd(90.0f);
                          const __m512d _0         = _mm512_setzero_pd();
                          const __m512d dif        = _mm512_sub_pd(_90,alp);
                          const __mmask16 m1      = _mm512_cmp_pd_mask(_0,tht,_CMP_LE_OQ);
                          const __mmask16 m2      = _mm512_cmp_pd_mask(tht,dif,_CMP_LT_OQ);
                          const __mmask16 m3      = _mm512_cmp_pd_mask(tht,dif,_CMP_EQ_OQ);
                          const __mmask16 m4      = _mm512_cmp_pd_mask(dif,tht,_CMP_LT_OQ);
                          const __mmask16 m5      = _mm512_cmp_pd_mask(tht,_90,_CMP_LE_OQ);
                          if(m1 && m2) {
                              return (rcs_f526_zmm8r8(gam0,alp,tht));
                          }
                          else if(m3) {
                              return (rcs_f527_zmm8r8(b,alp));
                          }
                          else(m4 && m5) {
                               return (rcs_f528_zmm8r8(b,r0,tht)); 
                          }
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f52678_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                               const double * __restrict __ATTR_ALIGN__(64) palp,
                                               const double * __restrict __ATTR_ALIGN__(64) ptht,
                                               const double * __restrict __ATTR_ALIGN__(64) pb,
                                               const double * __restrict __ATTR_ALIGN__(64) pr0) {

                          const __m512d _90        = _mm512_set1_pd(90.0f);
                          const __m512d _0         = _mm512_setzero_pd();
                          const __m512d dif        = _mm512_sub_pd(_90,alp);
                          const __mmask16 m1      = _mm512_cmp_pd_mask(_0,tht,_CMP_LE_OQ);
                          const __mmask16 m2      = _mm512_cmp_pd_mask(tht,dif,_CMP_LT_OQ);
                          const __mmask16 m3      = _mm512_cmp_pd_mask(tht,dif,_CMP_EQ_OQ);
                          const __mmask16 m4      = _mm512_cmp_pd_mask(dif,tht,_CMP_LT_OQ);
                          const __mmask16 m5      = _mm512_cmp_pd_mask(tht,_90,_CMP_LE_OQ);
                          if(m1 && m2) {
                              return (rcs_f526_zmm8r8_a(pgam0,palp,ptht));
                          }
                          else if(m3) {
                              return (rcs_f527_zmm8r8_a(pb,palp));
                          }
                          else(m4 && m5) {
                               return (rcs_f528_zmm8r8_a(pb,pr0,ptht)); 
                          }
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f52678_zmm8r8_u(const double * __restrict  pgam0,
                                               const double * __restrict  palp,
                                               const double * __restrict  ptht,
                                               const double * __restrict  pb,
                                               const double * __restrict  pr0) {

                          const __m512d _90        = _mm512_set1_pd(90.0f);
                          const __m512d _0         = _mm512_setzero_pd();
                          const __m512d dif        = _mm512_sub_pd(_90,alp);
                          const __mmask16 m1      = _mm512_cmp_pd_mask(_0,tht,_CMP_LE_OQ);
                          const __mmask16 m2      = _mm512_cmp_pd_mask(tht,dif,_CMP_LT_OQ);
                          const __mmask16 m3      = _mm512_cmp_pd_mask(tht,dif,_CMP_EQ_OQ);
                          const __mmask16 m4      = _mm512_cmp_pd_mask(dif,tht,_CMP_LT_OQ);
                          const __mmask16 m5      = _mm512_cmp_pd_mask(tht,_90,_CMP_LE_OQ);
                          if(m1 && m2) {
                              return (rcs_f526_zmm8r8_u(pgam0,palp,ptht));
                          }
                          else if(m3) {
                              return (rcs_f527_zmm8r8_u(pb,palp));
                          }
                          else(m4 && m5) {
                               return (rcs_f528_zmm8r8_u(pb,pr0,ptht)); 
                          }
                }


                  /*
                      Long thin bodies of revolution.
                      Axial-incidence backscatter RCS.
                      Long thin parabolic ogive - MISSILE NOISE CONE.
                      Formula 5.2-10
                  */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5210_zmm8r8(const __m512d gam0,
                                            const __m512d alp,
                                            const __m512d k0a) {

                        const __m512d _16pi  = _mm512_set1_pd(50.265482457436691815402294132472);
                        const __m512d _1     = _mm512_set1_pd(1.0f);  
                        const __m512d _3     = _mm512_set1_pd(3.0f);
                        const __m512d _4     = _mm512_set1_pd(4.0f);
                        register __m512d rcs,gam2,talp,talp2,talp4,_2k0a;
                        register __m512d k0as,ear,eai,cer,cei,cabs,x0,x1;
                        register __m512d trm1,t0r,t0i,t1r,t1i,t2r,t2i,trm2r,trm2i;
                        gam2  = _mm512_mul_pd(gam0,gam0);
                        talp  = xtan(alp);
                        _2k0a = _mm512_add_pd(k0a,k0a);
                        k0as  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a));
                        talp2 = _mm512_mul_pd(talp,talp);
                        ear   = _mm512_setzero_pd();
                        talp4 = _mm512_mul_pd(talp2,talp2);
                        x0    = _mm512_div_pd(_3,_2k0a);
                        eai   = _mm512_add_pd(_2k0a,_2k0a);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        trm1  = _mm512_div_pd(_mm512_mul_pd(gam2,talp4),_16pi);
                        t0r   = _mm512_sub_pd(_1,cer);
                        t0i   = _mm512_sub_pd(_1,cei);
                        x1    = _mm512_div_pd(_3,k0as);
                        t2r   = _mm512_mul_pd(x1,t0r);
                        t1r   = _mm512_mul_pd(x0,_mm512_add_pd(_1,cer));
                        t2i   = _mm512_mul_pd(x1,t0i);
                        t1i   = _mm512_mul_pd(x0,_mm512_add_pd(_1,cei));
                        trm2r = _mm512_sub_pd(t0r,_mm512_sub_pd(t1r,t2r));
                        trm2i = _mm512_sub_pd(t0i,_mm512_sub_pd(t1i,t2i));
                        cabs  = cabs_zmm8r8(trm2r,trm2i);
                        rcs   = _mm512_mul_pd(trm1,cabs);
                        return (rcs);                 
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5210_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                               const double * __restrict __ATTR_ALIGN__(64) palp,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                        register __m512d  gam0     = _mm512_load_pd(&pgam0[0]);
                        register __m512d  alp      = _mm512_load_pd(&palp[0]);
                        register __m512d  k0s      = _mm512_load_pd(&pk0a[0]); 
                        const __m512d _16pi  = _mm512_set1_pd(50.265482457436691815402294132472);
                        const __m512d _1     = _mm512_set1_pd(1.0f);  
                        const __m512d _3     = _mm512_set1_pd(3.0f);
                        const __m512d _4     = _mm512_set1_pd(4.0f);
                        register __m512d rcs,gam2,talp,talp2,talp4,_2k0a;
                        register __m512d k0as,ear,eai,cer,cei,cabs,x0,x1;
                        register __m512d trm1,t0r,t0i,t1r,t1i,t2r,t2i,trm2r,trm2i;
                        gam2  = _mm512_mul_pd(gam0,gam0);
                        talp  = xtan(alp);
                        _2k0a = _mm512_add_pd(k0a,k0a);
                        k0as  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a));
                        talp2 = _mm512_mul_pd(talp,talp);
                        ear   = _mm512_setzero_pd();
                        talp4 = _mm512_mul_pd(talp2,talp2);
                        x0    = _mm512_div_pd(_3,_2k0a);
                        eai   = _mm512_add_pd(_2k0a,_2k0a);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        trm1  = _mm512_div_pd(_mm512_mul_pd(gam2,talp4),_16pi);
                        t0r   = _mm512_sub_pd(_1,cer);
                        t0i   = _mm512_sub_pd(_1,cei);
                        x1    = _mm512_div_pd(_3,k0as);
                        t2r   = _mm512_mul_pd(x1,t0r);
                        t1r   = _mm512_mul_pd(x0,_mm512_add_pd(_1,cer));
                        t2i   = _mm512_mul_pd(x1,t0i);
                        t1i   = _mm512_mul_pd(x0,_mm512_add_pd(_1,cei));
                        trm2r = _mm512_sub_pd(t0r,_mm512_sub_pd(t1r,t2r));
                        trm2i = _mm512_sub_pd(t0i,_mm512_sub_pd(t1i,t2i));
                        cabs  = cabs_zmm8r8(trm2r,trm2i);
                        rcs   = _mm512_mul_pd(trm1,cabs);
                        return (rcs);                 
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5210_zmm8r8_u( const double * __restrict  pgam0,
                                               const double * __restrict  palp,
                                               const double * __restrict  pk0a,) {

                        register __m512d  gam0     = _mm512_loadu_pd(&pgam0[0]);
                        register __m512d  alp      = _mm512_loadu_pd(&palp[0]);
                        register __m512d  k0s      = _mm512_loadu_pd(&pk0a[0]); 
                        const __m512d _16pi  = _mm512_set1_pd(50.265482457436691815402294132472);
                        const __m512d _1     = _mm512_set1_pd(1.0f);  
                        const __m512d _3     = _mm512_set1_pd(3.0f);
                        const __m512d _4     = _mm512_set1_pd(4.0f);
                        register __m512d rcs,gam2,talp,talp2,talp4,_2k0a;
                        register __m512d k0as,ear,eai,cer,cei,cabs,x0,x1;
                        register __m512d trm1,t0r,t0i,t1r,t1i,t2r,t2i,trm2r,trm2i;
                        gam2  = _mm512_mul_pd(gam0,gam0);
                        talp  = xtan(alp);
                        _2k0a = _mm512_add_pd(k0a,k0a);
                        k0as  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a));
                        talp2 = _mm512_mul_pd(talp,talp);
                        ear   = _mm512_setzero_pd();
                        talp4 = _mm512_mul_pd(talp2,talp2);
                        x0    = _mm512_div_pd(_3,_2k0a);
                        eai   = _mm512_add_pd(_2k0a,_2k0a);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        trm1  = _mm512_div_pd(_mm512_mul_pd(gam2,talp4),_16pi);
                        t0r   = _mm512_sub_pd(_1,cer);
                        t0i   = _mm512_sub_pd(_1,cei);
                        x1    = _mm512_div_pd(_3,k0as);
                        t2r   = _mm512_mul_pd(x1,t0r);
                        t1r   = _mm512_mul_pd(x0,_mm512_add_pd(_1,cer));
                        t2i   = _mm512_mul_pd(x1,t0i);
                        t1i   = _mm512_mul_pd(x0,_mm512_add_pd(_1,cei));
                        trm2r = _mm512_sub_pd(t0r,_mm512_sub_pd(t1r,t2r));
                        trm2i = _mm512_sub_pd(t0i,_mm512_sub_pd(t1i,t2i));
                        cabs  = cabs_zmm8r8(trm2r,trm2i);
                        rcs   = _mm512_mul_pd(trm1,cabs);
                        return (rcs);                 
                }


                 /*
                      Long thin bodies of revolution.
                      Axial-incidence backscatter RCS.
                      High-frequency limit.
                      Long thin parabolic ogive - MISSILE NOISE CONE.
                      Formula 5.2-12
                  */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5212_zmm8r8(const __m512d gam0,
                                            const __m512d alp,
                                            const __m512d k0a) {

                          const __m512d _1     = _mm512_set1_pd(1.0f);  
                          const __m512d _4     = _mm512_set1_pd(4.0f);
                          const __m512d _8pi   = _mm512_set1_pd(25.132741228718345907701147066236);
                          register __m512d rcs,gam2,talp,talp2,talp4,arg,carg,x0,rat;
                          gam2 = _mm512_mul_pd(gam0,gam0);
                          arg = _mm512_mul_pd(_4,k0a);
                          talp = xtan(alp);
                          carg = xcos(arg);
                          talp2= _mm512_mul_pd(talp,talp);
                          x0   = _mm512_add_pd(_1,carg);
                          talp4= _mm512_mul_pd(talp2,talp2);
                          rat  = _mm512_div_pd(_mm512_mul_pd(gam2,talp4),_8pi);
                          rcs  = _mm512_mul_pd(rat,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5212_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                               const double * __restrict __ATTR_ALIGN__(64) palp,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d  gam0     = _mm512_load_pd(&pgam0[0]);
                          register __m512d  alp      = _mm512_load_pd(&palp[0]);
                          register __m512d  k0s      = _mm512_load_pd(&pk0a[0]); 
                          const __m512d _1     = _mm512_set1_pd(1.0f);  
                          const __m512d _4     = _mm512_set1_pd(4.0f);
                          const __m512d _8pi   = _mm512_set1_pd(25.132741228718345907701147066236);
                          register __m512d rcs,gam2,talp,talp2,talp4,arg,carg,x0,rat;
                          gam2 = _mm512_mul_pd(gam0,gam0);
                          arg = _mm512_mul_pd(_4,k0a);
                          talp = xtan(alp);
                          carg = xcos(arg);
                          talp2= _mm512_mul_pd(talp,talp);
                          x0   = _mm512_add_pd(_1,carg);
                          talp4= _mm512_mul_pd(talp2,talp2);
                          rat  = _mm512_div_pd(_mm512_mul_pd(gam2,talp4),_8pi);
                          rcs  = _mm512_mul_pd(rat,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5212_zmm8r8_u(const double * __restrict  pgam0,
                                               const double * __restrict palp,
                                               const double * __restrict  pk0a) {

                          register __m512d  gam0     = _mm512_loadu_pd(&pgam0[0]);
                          register __m512d  alp      = _mm512_loadu_pd(&palp[0]);
                          register __m512d  k0s      = _mm512_loadu_pd(&pk0a[0]); 
                          const __m512d _1     = _mm512_set1_pd(1.0f);  
                          const __m512d _4     = _mm512_set1_pd(4.0f);
                          const __m512d _8pi   = _mm512_set1_pd(25.132741228718345907701147066236);
                          register __m512d rcs,gam2,talp,talp2,talp4,arg,carg,x0,rat;
                          gam2 = _mm512_mul_pd(gam0,gam0);
                          arg = _mm512_mul_pd(_4,k0a);
                          talp = xtan(alp);
                          carg = xcos(arg);
                          talp2= _mm512_mul_pd(talp,talp);
                          x0   = _mm512_add_pd(_1,carg);
                          talp4= _mm512_mul_pd(talp2,talp2);
                          rat  = _mm512_div_pd(_mm512_mul_pd(gam2,talp4),_8pi);
                          rcs  = _mm512_mul_pd(rat,x0);
                          return (rcs);
                 }


                    /*
                         Circular ogive nose-on RCS (en emprical solution).
                         Formula 5.2-14
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5214_zmm8r8(const __m512d gam0,
                                            const __m512d alp,
                                            const __m512d k0a) {

                          const __m512d _1     = _mm512_set1_pd(1.0f);  
                          const __m512d _4     = _mm512_set1_pd(4.0f);
                          const __m512d _2pi   = _mm512_set1_pd(6.283185307179586476925286766559);
                          const __m512d pi     = _mm512_set1_pd(3.14159265358979323846264338328f;
                          const __m512d c0     = _mm512_set1_pd(20000000000.0);
                          const __m512d c1     = _mm512_set1_pd(4.3);
                          const __m512d c2     = _mm512_set1_pd(1.25f);
                          register __m512d rcs,gam2,alp2,alp2p,s2alp,arg,sarg;
                          register __m512d num,den,rat,num1,den1;
                          gam2  = _mm512_mul_pd(gam0,gam0);
                          alp2  = _mm512_add_pd(alp,alp);
                          s2alp = xsin(alp2);
                          num   = _mm512_mul_pd(_mm512_mul_pd(_4,k0a),alp);
                          den   = _mm512_mul_pd(pi,s2alp);
                          alp2p = _mm512_pow_pd(alp2,c1);
                          rat   = _mm512_sub_pd(_mm512_div_pd(num,den),c2);
                          num1  = _mm512_mul_pd(pi,_mm512_mul_pd(gam2,alp2p));
                          arg   = _mm512_mul_pd(_2pi,rat);
                          den1  = _mm512_div_pd(num1,c0);
                          sarg  = xsin(arg);
                          rcs   = _mm512_mul_pd(den1,_mm512_add_pd(_1,sarg));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5214_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const double * __restrict __ATTR_ALIGN__(64) palp,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d  gam0     = _mm512_load_pd(&pgam0[0]);
                          register __m512d  alp      = _mm512_load_pd(&palp[0]);
                          register __m512d  k0s      = _mm512_load_pd(&pk0a[0]); 
                          const __m512d _1     = _mm512_set1_pd(1.0f);  
                          const __m512d _4     = _mm512_set1_pd(4.0f);
                          const __m512d _2pi   = _mm512_set1_pd(6.283185307179586476925286766559);
                          const __m512d pi     = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d c0     = _mm512_set1_pd(20000000000.0);
                          const __m512d c1     = _mm512_set1_pd(4.3);
                          const __m512d c2     = _mm512_set1_pd(1.25f);
                          register __m512d rcs,gam2,alp2,alp2p,s2alp,arg,sarg;
                          register __m512d num,den,rat,num1,den1;
                          gam2  = _mm512_mul_pd(gam0,gam0);
                          alp2  = _mm512_add_pd(alp,alp);
                          s2alp = xsin(alp2);
                          num   = _mm512_mul_pd(_mm512_mul_pd(_4,k0a),alp);
                          den   = _mm512_mul_pd(pi,s2alp);
                          alp2p = _mm512_pow_pd(alp2,c1);
                          rat   = _mm512_sub_pd(_mm512_div_pd(num,den),c2);
                          num1  = _mm512_mul_pd(pi,_mm512_mul_pd(gam2,alp2p));
                          arg   = _mm512_mul_pd(_2pi,rat);
                          den1  = _mm512_div_pd(num1,c0);
                          sarg  = xsin(arg);
                          rcs   = _mm512_mul_pd(den1,_mm512_add_pd(_1,sarg));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f5214_zmm8r8_u(const double * __restrict pgam0,
                                              const double * __restrict  palp,
                                              const double * __restrict  pk0a) {

                          register __m512d  gam0     = _mm512_loadu_pd(&pgam0[0]);
                          register __m512d  alp      = _mm512_loadu_pd(&palp[0]);
                          register __m512d  k0s      = _mm512_loadu_pd(&pk0a[0]); 
                          const __m512d _1     = _mm512_set1_pd(1.0f);  
                          const __m512d _4     = _mm512_set1_pd(4.0f);
                          const __m512d _2pi   = _mm512_set1_pd(6.283185307179586476925286766559);
                          const __m512d pi     = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d c0     = _mm512_set1_pd(20000000000.0);
                          const __m512d c1     = _mm512_set1_pd(4.3);
                          const __m512d c2     = _mm512_set1_pd(1.25f);
                          register __m512d rcs,gam2,alp2,alp2p,s2alp,arg,sarg;
                          register __m512d num,den,rat,num1,den1;
                          gam2  = _mm512_mul_pd(gam0,gam0);
                          alp2  = _mm512_add_pd(alp,alp);
                          s2alp = xsin(alp2);
                          num   = _mm512_mul_pd(_mm512_mul_pd(_4,k0a),alp);
                          den   = _mm512_mul_pd(pi,s2alp);
                          alp2p = _mm512_pow_pd(alp2,c1);
                          rat   = _mm512_sub_pd(_mm512_div_pd(num,den),c2);
                          num1  = _mm512_mul_pd(pi,_mm512_mul_pd(gam2,alp2p));
                          arg   = _mm512_mul_pd(_2pi,rat);
                          den1  = _mm512_div_pd(num1,c0);
                          sarg  = xsin(arg);
                          rcs   = _mm512_mul_pd(den1,_mm512_add_pd(_1,sarg));
                          return (rcs);
                 }



                  



       }



}














#endif /*__GMS_RCS_ELLIPS_OGIVE_ZMM8R8_HPP__*/
