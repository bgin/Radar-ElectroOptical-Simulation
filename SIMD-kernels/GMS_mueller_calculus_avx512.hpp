

#ifndef __GMS_MUELLER_CALCULUS_AVX512_HPP__
#define __GMS_MUELLER_CALCULUS_AVX512_HPP__
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


namespace file_info {

 const unsigned int gGMS_MUELLER_CALCULUS_AVX512_MAJOR = 1U;
 const unsigned int gGMS_MUELLER_CALCULUS_AVX512_MINOR = 0U;
 const unsigned int gGMS_MUELLER_CALCULUS_AVX512_MICRO = 0U;
 const unsigned int gGMS_MUELLER_CALCULUS_AVX512_FULLVER =
  1000U*gGMS_MUELLER_CALCULUS_AVX512_MAJOR+100U*gGMS_MUELLER_CALCULUS_AVX512_MINOR+10U*gGMS_MUELLER_CALCULUS_AVX512_MICRO;
 const char * const pgGMS_MUELLER_CALCULUS_AVX512_CREATION_DATE = "22-12-2021 12:50 +00200 (WED 22 DEC 2021 12:50 GMT+2)";
 const char * const pgGMS_MUELLER_CALCULUS_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_MUELLER_CALCULUS_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_MUELLER_CALCULUS_AVX512_SYNOPSIS      = "AVX512 based Mueller calculus implementation.";


}


#include <complex>
#include "GMS_mueller_types_avx512.hpp"
#include "GMS_config.h"


namespace gms {

        namespace  math {


	                 // Based on: R.A. Chipman, "Polarimetry" chapter in Handbook of Optics Volume 2
                         // (McGraw-Hill, New York, 1995).

	               /*
                                Jones-Vector components -- packed single-precision
                                16-tuple vectors of complex numbers [deinterleaved]
                           */
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x16c16
                      JVec2x16c16_set_1() {

		           JVec2x16c16 v;
			   v.p = ZMM16c4(); // 'p' wave component
			   v.s = ZMM16c4(); // 's' wave component
			   return (v);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x16c16
                      JVec2x16c16_set_2(const std::complex<float> c1,
		                        const std::complex<float> c2) {

                          JVec2x16c16 v;
			  v.p = ZMM16c4(c1);
			  v.s = ZMM16c4(c2);
			  return (v);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x16c16
                      JVec2x16c16_set_3(const float * __restrict __ATTR_ALIGN__(64) re1,
		                        const float * __restrict __ATTR_ALIGN__(64) im1,
					const float * __restrict __ATTR_ALIGN__(64) re2,
					const float * __restrict __ATTR_ALIGN__(64) im2) {

                          JVec2x16c16 v;
			  v.p = ZMM16c4(re1,im1);
			  v.s = ZMM16c4(re2,im2);
			  return (v);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x16c16
                      JVec2x16c16_set_4(const float re1,
		                        const float im1,
					const float re2,
					const float im2) {

                           JVec2x16c16 v;
			   v.p = ZMM16c4(re1,im1);
			   v.s = ZMM16c4(re2,im2);
			   return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x16c16
                    JVec2x16c16_set_5(const __m512 re1,
		                      const __m512 im1,
				      const __m512 re2,
				      const __m512 im2) {

                           JVec2x16c16 v;
			   v.p = ZMM16c4(re1,im1);
			   v.s = ZMM16c4(re2,im2);
			   return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x16c16
                    JVec2x16c16_set_6(const JVec2x16c16 x) {

                          JVec2x16c16 v;
			  v.p = x.j0;
			  v.s = x.j1;
			  return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x16c16
                    JVec2x16c16_set_7(const ZMM16c4 s,
		                      const ZMM16c4 p) {

                         JVec2x16c16 v;
			 v.p = s;
			 v.s = p;
			 return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x16c16
                    JVec2x16c16_mul_ZMM16c4(const JVec2x16c16 x
		                            const ZMM16c4 y) {

                         JVec2x16c16 v;
			 v.p = x.p*y;
			 v.s = x.s*y;
			 return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    void
		    JVec2x16c16_mul_ZMM16c4(JVec2x16c16 & x,
		                            const ZMM16c4 y) {

                         x.p *= y;
			 x.s *= y;
		  }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x16c16
                    JVec2x16c16_div_ZMM16c4(const JVec2x16c16 x,
		                            const ZMM16c4 y) {

                         JVec2x16c16 v;
			 v.p = x.p/y;
			 v.s = x.s/y;
			 return (v);
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
                   void
		   JVec2x16c16_div_ZMM16c4(JVec2x16c16 & x,
		                           const ZMM16c4 y) {

                         x.p /= y;
			 x.s /= y;
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
		   ZMM16c4
		   JVec2x16c16_mul_JVec2x16c16(const JVec2x16c16 x,
		                               const JVec2x16c16 y) {

                         ZMM16c4 c;
			 c = x.p*y.p+x.s*y.s;
			 return (c);
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
		   JVec2x16c16
		   JVec2x16c16_conj(const JVec2x16c16 x) {

		          JVec2x16c16 v;
			  v.p = conj(x.p);
			  v.s = conj(x.s);
			  return (v);
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
		   JVec2x16c16
		   JVec2x16c16_add_JVec2x16c16(const JVec2x16c16 x,
		                               const JVec2x16c16 y) {

                          JVec2x16c16 v;
			  v.p = x.p+y.p;
			  v.s = x.s+y.s;
			  return (v);
		 }


		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  void
		  JVec2x16c16_add_JVec2x16c16(JVec2x16c16 & x,
		                              const JVec2x16c16 y) {

                         x.p += y.p;
			 x.s += y.s;
		 }


		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  JVec2x16c16
		  JVec2x16c16_sub_JVec2x16c16(const JVec2x16c16 x,
		                               const JVec2x16c16 y) {

                          JVec2x16c16 v;
			  v.p = x.p-y.p;
			  v.s = x.s-y.s;
			  return (v);
		 }


		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  void
		  JVec2x16c16_sub_JVec2x16c16(JVec2x16c16 & x,
		                              const JVec2x16c16 y) {

                         x.p -= y.p;
			 x.s -= y.s;
		 }


		 // An intensity (vector norm)
		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  __m512
		  JVec2x16c16_norm(const JVec2x16c16 x) {
 
                     __m512 t0 = _mm512_fmadd_ps(x.p,x.p,
		                          _mm512_mul_ps(x.s,x.s));
		     return (_mm512_sqrt_ps(t0));
		 }


		  __ATTR_REGCALL__
		  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  JVec2x16c16
		  JVec2x16c16_negative(const JVec2x16c16 x) {

                       JVec2x16c16 v;
		       v.p = ~x.p;
		       v.s = ~x.s;
		       return (v);
		 }


		  __ATTR_REGCALL__
		  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  __m512
		  JVec2x16c16_psi(const JVec2x16c16 x) {

                      __m512 t0 = _mm512_div_ps(cabs(x.p),cabs(x.s));
		      return (_mm512_atan_ps(t0));
		 }


		  __ATTR_REGCALL__
		  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  __m512
		  JVec2x16c16_delta(const JVec2x16c16 x) {

                       return (_mm512_sub_ps(carg(x.s),carg(x.p)));
		 }


	        /*
                                Jones-Vector components -- packed double-precision
                                16-8 vectors of complex numbers [deinterleaved]
                               
                */


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x8c8
                      JVec2x8c8_set_1() {

		           JVec2x8c8 v;
			   v.p = ZMM8c8(); // 'p' wave component
			   v.s = ZMM8c8(); // 's' wave component
			   return (v);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x8c8
                      JVec2x8c8_set_2(const std::complex<double> c1,
		                      const std::complex<double> c2) {

                          JVec2x8c8 v;
			  v.p = ZMM8c8(c1);
			  v.s = ZMM8c8(c2);
			  return (v);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x8c9
                      JVec2x8c8_set_3(  const double * __restrict __ATTR_ALIGN__(64) re1,
		                        const double * __restrict __ATTR_ALIGN__(64) im1,
					const double * __restrict __ATTR_ALIGN__(64) re2,
					const double * __restrict __ATTR_ALIGN__(64) im2) {

                          JVec2x8c8 v;
			  v.p = ZMM8c8(re1,im1);
			  v.s = ZMM8c8(re2,im2);
			  return (v);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_ALIGN__(32)
		      static inline
		      JVec2x8c8
                      JVec2x8c8_set_4(  const double re1,
		                        const double im1,
					const double re2,
					const double im2) {

                           JVec2x8c8 v;
			   v.p = ZMM8c8(re1,im1);
			   v.s = ZMM8c8(re2,im2);
			   return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x8c8
                    JVec2x8c8_set_5(  const __m512d re1,
		                      const __m512d im1,
				      const __m512d re2,
				      const __m512d im2) {

                           JVec2x8c8 v;
			   v.p = ZMM8c8(re1,im1);
			   v.s = ZMM8c8(re2,im2);
			   return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x8c8
                    JVec2x8c8_set_6(const JVec2x8c8 x) {

                          JVec2x16c16 v;
			  v.p = x.j0;
			  v.s = x.j1;
			  return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x16c16
                    JVec2x16c16_set_7(const ZMM16c4 s,
		                      const ZMM16c4 p) {

                         JVec2x8c8 v;
			 v.p = s;
			 v.s = p;
			 return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x8c8
                    JVec2x8c8_mul_ZMM8c8(const JVec2x8c8 x
		                         const ZMM8c8 y) {

                         JVec2x8c8 v;
			 v.p = x.p*y;
			 v.s = x.s*y;
			 return (v);
		   }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    void
		    JVec2x8c8_mul_ZMM8c8(JVec2x8c8 & x,
		                         const ZMM8c8 y) {

                         x.p *= y;
			 x.s *= y;
		  }


		    __ATTR_REGCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_ALIGN__(32)
		    static inline
		    JVec2x8c8
                    JVec2x8c8_div_ZMM8c8(const JVec2x8c8 x,
		                            const ZMM8c8 y) {

                         JVec2x8c8 v;
			 v.p = x.p/y;
			 v.s = x.s/y;
			 return (v);
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
                   void
		   JVec2x8c8_div_ZMM8c8(JVec2x8c8 & x,
		                           const ZMM8c8 y) {

                         x.p /= y;
			 x.s /= y;
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
		   ZMM8c8
		   JVec2x8c8_mul_JVec2x8c8(const JVec2x8c8 x,
		                               const JVec2x8c8 y) {

                         ZMM8c8 c;
			 c = x.p*y.p+x.s*y.s;
			 return (c);
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
		   JVec2x8c8
		   JVec2x8c8_conj(const JVec2x8c8 x) {

		          JVec2x8c8 v;
			  v.p = conj(x.p);
			  v.s = conj(x.s);
			  return (v);
		  }


		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
		   JVec2x8c8
		   JVec2x8c8_add_JVec2x8c8(const JVec2x8c8 x,
		                               const JVec2x8c8 y) {

                          JVec2x8c8 v;
			  v.p = x.p+y.p;
			  v.s = x.s+y.s;
			  return (v);
		 }


		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  void
		  JVec2x8c8_add_JVec2x8c8(JVec2x8c8 & x,
		                              const JVec2x8c8 y) {

                         x.p += y.p;
			 x.s += y.s;
		 }


		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  JVec2x8c8
		  JVec2x8c8_sub_JVec2x8c8(const JVec2x8c8 x,
		                               const JVec2x8c8 y) {

                          JVec2x8c8 v;
			  v.p = x.p-y.p;
			  v.s = x.s-y.s;
			  return (v);
		 }


		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  void
		  JVec2x8c8_sub_JVec2x8c8(JVec2x8c8 & x,
		                          const JVec2x8c8 y) {

                         x.p -= y.p;
			 x.s -= y.s;
		 }


		 // An intensity (vector norm)
		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  __m512
		  JVec2x8c8_norm(const JVec2x8c8 x) {
 
                     __m512 t0 = _mm512_fmadd_ps(x.p,x.p,
		                          _mm512_mul_ps(x.s,x.s));
		     return (_mm512_sqrt_ps(t0));
		 }


		  __ATTR_REGCALL__
		  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  JVec2x8c8
		  JVec2x8c8_negative(const JVec2x8c8 x) {

                       JVec2x8c8 v;
		       v.p = ~x.p;
		       v.s = ~x.s;
		       return (v);
		 }


		  __ATTR_REGCALL__
		  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  __m512
		  JVec2x8c8_psi(const JVec2x8c8 x) {

                      __m512 t0 = _mm512_div_ps(cabs(x.p),cabs(x.s));
		      return (_mm512_atan_ps(t0));
		 }


		  __ATTR_REGCALL__
		  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  __m512
		  JVec2x8c8_delta(const JVec2x8c8 x) {

                       return (_mm512_sub_ps(carg(x.s),carg(x.p)));
		 }


		 /*

                         Jones Matrix implementation based on SIMD 
                         16-tuple complex vector [deinterleaved]
                         single-precision.
                         @Reference

                           
                                  typedef struct  __ATTR_ALIGN__(64) JMat4x16c16 {

                                      ZMM16c4 pp;
				      ZMM16c4 ss;
				      ZMM16c4 ps;
				      ZMM16c4 sp;
		                 }JMat4x16c16;
                                 typedef struct __ATTR_ALIGN__(64)  JMat4x8c8 {

                                       ZMM8c8 pp;
				       ZMM8c8 ss;
				       ZMM8c8 ps;
				       ZMM8c8 sp;
		                }JMat4x8c8;
                  */

		   __ATTR_REGCALL__
                   __ATTR_ALWAYS_INLINE__
		   __ATTR_ALIGN__(32)
		   static inline
		   JMat4x16c16
		   JMat4x16c16_set_1() {
 
                      JMat4x16c16 v;
		      v.pp = ZMM16c4(); // 'pp' component
		      v.ss = ZMM16c4(); // 'ss' component
		      v.ps = ZMM16c4(); // 'ps' component
		      v.sp = ZMM16c4(); // 'sp' component
		      return (v);
		 }


		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  JMat4x16c16
		  JMat4x16c16_set_2(const std::complex<float> pp,
		                    const std::complex<float> ss,
				    const std::complex<float> ps,
				    const std::complex<float> sp) {

                       JMat4x16c16 v;
		       v.pp = ZMM16c4(pp);
		       v.ss = ZMM16c4(ss);
		       v.ps = ZMM16c4(ps);
		       v.sp = ZMM16c4(sp);
		       return (v);
		}
		

		  __ATTR_REGCALL__
                  __ATTR_ALWAYS_INLINE__
		  __ATTR_ALIGN__(32)
		  static inline
		  JMat4x16c16
		  JMat4x16c16_set_3(const float * __restrict __ATTR_ALIGN__(64) re1,
		                    const float * __restrict __ATTR_ALIGN__(64) im1,
				    const float * __restrict __ATTR_ALIGN__(64) re2,
				    const float * __restrict __ATTR_ALIGN__(64) im2,
				    const float * __restrict __ATTR_ALIGN__(64) re3,
				    const float * __restrict __ATTR_ALIGN__(64) im3,
				    const float * __restrict __ATTR_ALIGN__(64) re4,
				    const float * __restrict __ATTR_ALIGN__(64) im4) {

                       JMatrix4x16c16 v;
		       v.pp = ZMM16c4(re1,im1);
		       v.ss = ZMM16c4(re2,im2);
		       v.ps = ZMM16c4(re3,im3);
		       v.sp = ZMM16c4(re4,im4);
		       return (v);
		}


		 __ATTR_REGCALL__
                 __ATTR_ALWAYS_INLINE__
		 __ATTR_ALIGN__(32)
		 static inline
		 JMat4x16c16
		 JMat4x16c16_set_4(const float re1,
		                   const float im1,
				   const float re2,
				   const float im2,
				   const float re3,
				   const float im3,
				   const float re4,
				   const float im4) {

                       JMatrix4x16c16 v;
		       v.pp = ZMM16c4(re1,im1);
		       v.ss = ZMM16c4(re2,im2);
		       v.ps = ZMM16c4(re3,im3);
		       v.sp = ZMM16c4(re4,im4);
		       return (v);  
	        }


		 __ATTR_REGCALL__
                 __ATTR_ALWAYS_INLINE__
		 __ATTR_ALIGN__(32)
		 static inline
		 JMat4x16c16
		 JMat4x16c16_set_5(const __m512 re1,
		                   const __m512 im1,
				   const __m512 re2,
				   const __m512 im2,
				   const __m512 re3,
				   const __m512 im3,
				   const __m512 re4,
				   const __m512 im4) {

                       JMatrix4x16c16 v;
		       v.pp = ZMM16c4(re1,im1);
		       v.ss = ZMM16c4(re2,im2);
		       v.ps = ZMM16c4(re3,im3);
		       v.sp = ZMM16c4(re4,im4);
		       return (v); 
	       }


	       	 __ATTR_REGCALL__
                 __ATTR_ALWAYS_INLINE__
		 __ATTR_ALIGN__(32)
		 static inline
		 JMat4x16c16
		 JMat4x16c16_set_6(const JMat4x16c16 x) {

                      JMat4x16c16 v;
		      v.pp = x.pp;
		      v.ss = x.ss;
		      v.ps = x.ps;
		      v.sp = x.sp;
		      return (v);
	       }


	        __ATTR_REGCALL__
                __ATTR_ALWAYS_INLINE__
	        __ATTR_ALIGN__(32)
		static inline
		JMat4x16c16
		JMat4x16c16_set_7(const ZMM16c4 pp,
		                  const ZMM16c4 ss,
				  const ZMM16c4 ps,
				  const ZMM16c4 sp) {

                     JMat4x16c16 v;
		     v.pp = pp;
		     v.ss = ss;
		     v.ps = ps;
		     v.sp = sp;
		     return (v);
	      }


	      // Jones Matrix multiplication
	      __ATTR_REGCALL__
              __ATTR_ALWAYS_INLINE__
	      __ATTR_ALIGN__(32)
	      static inline
	      JMat4x16c16
	      JMat4x16c16_mul_JMat4x16c16(const JMat4x16c16 m1,
	                                  const JMat4x16c16 m2) {

                    JMat4x16c16 v;
		    v.pp = m1.sp*m2.ps+m1.pp*m2.pp;
		    v.ss = m1.ss*m2.ss+m1.ps*m2.sp;
		    v.ps = m1.ss*m2.ps+m1.ps*m2.pp;
		    v.sp = m1.sp*m2.ss+m1.pp*m2.sp;
		    return (v);
	    }


	      __ATTR_REGCALL__
              __ATTR_ALWAYS_INLINE__
	      __ATTR_ALIGN__(32)
	      static inline
	      void
	      JMat4x16c16_mul_JMat4x16c16(JMat4x16c16 &m1,
	                                  const JMat4x16c16 m2) {

                    m1.pp = m1.sp*m2.ps+m1.pp*m2.pp;
		    m1.ss = m1.ss*m2.ss+m1.ps*m2.sp;
		    m1.ps = m1.ss*m2.ps+m1.ps*m2.pp;
		    m1.sp = m1.sp*m2.ss+m1.pp*m2.sp;
	   }


	   // Jones Matrix multiplied by constant
	     __ATTR_REGCALL__
             __ATTR_ALWAYS_INLINE__
	     __ATTR_ALIGN__(32)
	     static inline
	     JMat4x16c16
	     JMat4x16c16_mul_ZMM16c4(const JMat4x16c16 m
	                             const ZMM16c4 x) {

                   JMat4x16c16 v;
		   v.pp = m.pp*x;
		   v.ss = m.ss*x;
		   v.ps = m.ps*x;
		   v.sp = m.sp*x;
		   return (v);
	   }


	    __ATTR_REGCALL__
            __ATTR_ALWAYS_INLINE__
	    __ATTR_ALIGN__(32)
	    static inline
            void
	    JMat4x16c16_mul_JMat4x16c16(JMat4x6c16 &m,
	                                const ZMM16c x) {

                   m.pp = m.pp*x;
		   m.ss = m.ss*x;
		   m.ps = m.ps*x;
		   m.sp = m.sp*x;
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JVec2x16c16
	  JVec2x16c16_mul_JMat4x16c16(const JVec2x16c16 x,
	                              const JMat4x16c16 y) {

               JVec2x16c16 v;
	       v.p = y.ss*x.p+y.sp*x.s;
	       v.s = y.sp*x.p+y.pp*x.s;
	       return (v);
	  }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JMat4x16c16
	  JMat4x16c16_div_ZMM16c4( const JMat4x16c16 m
	                           const ZMM16c4 x) {

               JMat4x16c16 v;
	       v.pp = m.pp/x;
	       v.ss = m.ss/x;
	       v.ps = m.ps/x;
	       v.sp = m.sp/x;
	       return (v);
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  void
	  JMat4x16c16_div_ZMM16c4(JMat4x16c16 &m,
	                          const ZMM16c4 x) {

               m.pp = m.pp/x;
	       m.ss = m.ss/x;
	       m.ps = m.ps/x;
	       m.sp = m.sp/x;
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JMat4x16c16
	  JMat4x16c16_add_JMat4x16c16(const JMat4x16c16 x,
	                              const JMat4x16c16 y) {

               JMat4x16c16 v;
	       v.pp = x.pp+y.pp;
	       v.ss = x.ss+y.ss;
	       v.ps = x.ps+y.ps;
	       v.sp = x.sp+y.sp;
	       return (v);
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  void
	  JMat4x16c16_add_JMat4x16c16(JMat4x16c16 &x,
	                              const JMat4x16c16 y) {

               x.pp += y.pp;
	       x.ss += y.ss;
	       x.ps += y.ps;
	       x.sp += y.sp;
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JMat4x16c16
	  JMat4x16c16_sub_JMat4x16c16(const JMat4x16c16 x,
	                              const JMat4x16c16 y) {

               JMat4x16c16 v;
	       v.pp = x.pp-y.pp;
	       v.ss = x.ss-y.ss;
	       v.ps = x.ps-y.ps;
	       v.sp = x.sp-y.sp;
	       return (v);
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  void
	  JMat4x16c16_sub_JMat4x16c16(JMat4x16c16 &x,
	                              const JMat4x16c16 y) {

                x.pp -= y.pp;
		x.ss -= y.ss;
		x.ps -= y.ps;
		x.sp -= y.sp;
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JMat4x16c16
	  JMat4x16c16_negate(const JMat4x16c16 x) {

                JMat4x16c16 v;
		v.pp = -x.pp;
		v.ss = -x.ss;
		v.ps = -x.ps;
		v.sp = -x.sp;
		return (v);
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JMat4x16c16
	  JMat4x16c16_negate(JMat4x16c16 &x) {

               x.pp = -x.pp;
	       x.ss = -x.ss;
	       x.ps = -x.ps;
	       x.sp = -x.sp;
	  }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JMat4x16c16
	  JMat4x16c16_rotator(const __m512 theta) {

              const __m512 vc  = _mm512_cos_ps(theta);
	      const __m512 vs  = _mm512_sin_ps(theta);
	      const __m512 nvs = _mm512_sub_ps(_mm512_setzero_ps(),vs);
	      JMat4x16c16 v;
	      //v.pp.re = vc;
	      //v.ss.re = vc;
	      //v.ps.re = vs;
	      //v.sp.re = vns;
	      v.pp = ZMM16c4(vc);
	      v.ss = ZMM16c4(vc);
	      v.ps = ZMM16c4(vs);
	      v.sp = ZMM16c4(nvs);
	      return (v);
	 }


	  __ATTR_REGCALL__
          __ATTR_ALWAYS_INLINE__
	  __ATTR_ALIGN__(32)
	  static inline
	  JMat4x16c16
	  JMat4x16c16_linear_retarder(const __m512 phi,
	                              const __m512 ang) {

	     const __m512 _n1   = _mm512_set1_ps(-1.0F);
	     const __m512 _0    = _mm512_setzero_ps();
	     const __m512 _0_5  = _mm512_set1_ps(0.5F);
	     const __m512 _1    = _mm512_set1_ps(1.0F);
             const ZMM16c4 j    = ZMM16c4(_0,_n1);
	     const __m512 h_phi = _mm512_mul_ps(phi,_0_5);
	                     
	     const ZMM16c4 phasor = cexp(j*h_phi);
	     JMat4x16c16 v;
	     v.pp = phasor;
	     v.ss = _1/phasor;
	     v.ps = ZMM16c4();
	     v.sp = v.ps;
	     return (v);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 JMat4x16c16
	 JMat4x16c16_circular_retarder(const __m512 phi) {

	     const __m512 _0_5 = _mm512_set1_ps(0.5F);
	     const __m512 h_phi = _mm512_mul_ps(_0_5,phi);
             JMat4x16c16 v;
	     v.pp = ZMM16c4(_mm512_cos_ps(h_phi));
	     v.ss = v.pp;
	     v.ps = ZMM16c4(_mm512_sin_ps(h_phi));
	     v.sp = ZMM16c4(_mm512_sub_ps(_mm512_setzero_ps(),
	                                  v.ps.re));
	     return (v);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 JMat4x16c16
	 JMat4x16c16_circular_polarizer(const __m512 attenuation) {

             const __m512 _1   = _mm512_set1_ps(1.0F);
	     const __m512 _0_5 = _mm512_set1_ps(0.5F);
	     const __m512 _0   = _mm512_setzero_ps();
	     const __m512 t0   = _mm512_div_ps(_mm512_sub_ps(attenuation,_1),
	                                     _mm512_add_ps(attenuation,_1));
	     const __m512 e    = _mm512_sqrt_ps(t0);
	     const __m512 t2   = _mm512_mul_ps(_mm512_add_ps(_1,e),_0_5);
	     const __m512 t1   = _mm512_mul_ps(_mm512_sub_ps(_1,e),_0_5);
	     JMat4x16c16 v;
	     v.pp = ZMM16c4(t2);
	     v.ss = ZMM16c4(t2);
	     v.ps = ZMM16c4(_0,t1);
	     v.sp = ZMM16c4(_0,_mm512_sub_ps(_0,t1));
	     return (v);
	                    
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 JMat4x16c16
	 JMat4x16c16_eigenvalues(const JVec2x16c16 a,
	                         const JVec2x16c16 b,
				 const ZMM16c4 ca,
				 const ZMM16c4 cb) {

               JMat4x16c16 v;
	      //const ZMM16C4 t0 = a.p*b.s*cb-a.s*b.p*ca;
	      const ZMM16c4 t0 = a.p*b.s;
	      const ZMM16c4 t1 = a.s*b.p;
	      
	      v.pp = t0*cb-t1*ca;
              v.ss = t0*ca-t1*cb;
	      v.ps = a.p*b.p*(cb-ca);
	      const ZMM16c4 det = a.p*b.s-a.s*b.p;
	      v.sp = a.s*b.s*(ca-cb);
	      return (v/det);
        }


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 JMat4x16c16
	 JMat4x16c16_transpose(const JMat4x16c16 x) {

             JMat4x16c16 v;
	     v.pp = x.pp;
	     v.ss = x.ss;
	     v.ps = x.sp;
	     v.sp = x.ps;
	     return (v);
	}


         __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 JMat4x16c16
	 JMat4x16c16_hermitian(const JMat4x16c16 x) {

              JMat4x16c16 v;
	      v.pp = conj(x.pp);
	      v.ss = conj(x.ss);
	      v.ps = conj(x.sp);
	      v.sp = conj(x.ps);
	      return (v);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 JMat4x16c16
	 JMat4x16c16_conjugate(const JMat4x16c16 x) {

              JMat4x16c16 v;
	      v.pp = conj(x.pp);
	      v.ss = conj(x.ss);
	      v.ps = conj(x.ps);
	      v.sp = conj(x.sp);
	      return (v);
	}

         __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 SVec4x16v16
	 SVec4x16v16_set_1() {

             SVec4x16v16 sv;
	     sv.s0 = AVX512Vec16();
	     sv.s1 = AVX512Vec16();
	     sv.s2 = AVX512Vec16();
	     sv.s3 = AVX512Vec16();
	     return (sv);
	}


	


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 SVec4x16v16
	 SVec4x16v16_set_2(const float s0,
	                   const float s1,
			   const float s2,
			   const float s3) {

             SVec4x16v16 sv;
	     sv.s0 = AVX512Vec16(s0);
	     sv.s1 = AVX512Vec16(s1);
	     sv.s2 = AVX512Vec16(s2);
	     sv.s3 = AVX512Vec16(s3);
	     return (sv);
	}
	

	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 SVec4x16v16_set_2(SVec4x16v16 &sv
	                   const float s0,
	                   const float s1,
			   const float s2,
			   const float s3) {

             
	     sv.s0 = AVX512Vec16(s0);
	     sv.s1 = AVX512Vec16(s1);
	     sv.s2 = AVX512Vec16(s2);
	     sv.s3 = AVX512Vec16(s3);
	   
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 SVec4x16v16
	 SVec4x16v16_set_3(
	                   const float * __restrict __ATTR_ALIGN__(64) s0,
	                   const float * __restrict __ATTR_ALIGN__(64) s1,
			   const float * __restrict __ATTR_ALIGN__(64) s2,
			   const float * __restrict __ATTR_ALIGN__(64) s3) {

             SVec4x16v16 sv;
	     sv.s0 = AVX512Vec16(s0);
	     sv.s1 = AVX512Vec16(s1);
	     sv.s2 = AVX512Vec16(s2);
	     sv.s3 = AVX512Vec16(s3);
	     return (sv);
       }


         __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         void
	 SVec4x16v16_set_3(SVec4x16v16 &sv
	                   const float * __restrict __ATTR_ALIGN__(64) s0,
	                   const float * __restrict __ATTR_ALIGN__(64) s1,
			   const float * __restrict __ATTR_ALIGN__(64) s2,
			   const float * __restrict __ATTR_ALIGN__(64) s3) {

             
	     sv.s0 = AVX512Vec16(s0);
	     sv.s1 = AVX512Vec16(s1);
	     sv.s2 = AVX512Vec16(s2);
	     sv.s3 = AVX512Vec16(s3);
	   
       }


       	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 SVec4x16v16
	 SVec4x16v16_set_4(const __m512 s0,
	                   const __m512 s1,
			   const __m512 s2,
			   const __m512 s3) {

             SVec4x16v16 sv;
	     sv.s0 = AVX512Vec16(s0);
	     sv.s1 = AVX512Vec16(s1);
	     sv.s2 = AVX512Vec16(s2);
	     sv.s3 = AVX512Vec16(s3);
	     return (sv);
      }


      	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 SVec4x16v16_set_4(SVec4x16v16 &sv,
	                   const __m512 s0,
	                   const __m512 s1,
			   const __m512 s2,
			   const __m512 s3) {

             
	     sv.s0 = AVX512Vec16(s0);
	     sv.s1 = AVX512Vec16(s1);
	     sv.s2 = AVX512Vec16(s2);
	     sv.s3 = AVX512Vec16(s3);
	     
      }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_set_5(const SVec4x16v16 sv1) {

            SVec4x16v16 sv;
	    sv.s0 = sv1.s0;
	    sv.s1 = sv1.s1;
	    sv.s2 = sv1.s2;
	    sv.s3 = sv1.s3;
	    return (sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_set_5( SVec4x16v16 &sv
	                   const SVec4x16v16 sv1) {

            
	    sv.s0 = sv1.s0;
	    sv.s1 = sv1.s1;
	    sv.s2 = sv1.s2;
	    sv.s3 = sv1.s3;
	    
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_from_JVec2x16c16(const JVec2x16c16 jv) {

           SVec4x16c16 sv;
	   const ZMM16c4 t0 = cnorm(jv.p);
	   const ZMM16c4 t1 = cnorm(jv.s)
	   const ZMM16c4 t2 = conjugate(jv.p)*jv.s;
	   const __m512 _2  = _mm512_set1_ps(2.0F);
	   sv.s0 = t0+t1;
	   sv.s1 = t0-t1;
	   sv.s2 = _mm512_mul_ps(_2,t2.m_re);
	   sv.s3 = _mm512_mul_ps(_2,t2.m_im);
	   return (sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_from_JVec2x16c16(SVec4x16v16 &sv
	                             const JVec2x16c16 jv) {

           
	   const ZMM16c4 t0 = cnorm(jv.p);
	   const ZMM16c4 t1 = cnorm(jv.s)
	   const ZMM16c4 t2 = conjugate(jv.p)*jv.s;
	   const __m512 _2  = _mm512_set1_ps(2.0F);
	   sv.s0 = t0+t1;
	   sv.s1 = t0-t1;
	   sv.s2 = _mm512_mul_ps(_2,t2.m_re);
	   sv.s3 = _mm512_mul_ps(_2,t2.m_im);
	   
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_normalized_Q(const SVec4x16v16 jv) {

	      return (jv.s1/jv.s0);
        }


	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_normalized_U(const SVec4x16v16 jv) {

	      return (jv.s2/jv.s0);
        }


	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_normalized_V(const SVec4x16v16 jv) {

	      return (jv.s3/jv.s0);
        }


	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_add_SVec4x16v16(const SVec4x16v16 sv1,
	                            const SVec4x16v16 sv2) {

            SVec4x16v16 sv;
	    sv.s0 = sv1.s0+sv2.s0;
	    sv.s1 = sv1.s1+sv2.s1;
	    sv.s2 = sv1.s2+sv2.s2;
	    sv.s3 = sv1.s3+sv2.s3;
	    return (sv);
       }


       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_add_SVec4x16v16(SVec4x16v16 &sv,
	                            const SVec4x16v16 sv1,
	                            const SVec4x16v16 sv2) {

           
	    sv.s0 = sv1.s0+sv2.s0;
	    sv.s1 = sv1.s1+sv2.s1;
	    sv.s2 = sv1.s2+sv2.s2;
	    sv.s3 = sv1.s3+sv2.s3;
	   
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_add_SVec4x16v16(SVec4x16v16 &sv,
	                            const SVec4x16v16 sv1){
	                            

           
	    sv.s0 += sv1.s0;
	    sv.s1 += sv1.s1;
	    sv.s2 += sv1.s2;
	    sv.s3 += sv1.s3;
	   
       }


       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_sub_SVec4x16v16(const SVec4x16v16 sv1,
	                            const SVec4x16v16 sv2) {

            SVec4x16v16 sv;
	    sv.s0 = sv1.s0-sv2.s0;
	    sv.s1 = sv1.s1-sv2.s1;
	    sv.s2 = sv1.s2-sv2.s2;
	    sv.s3 = sv1.s3-sv2.s3;
	    return (sv);
       }


       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_sub_SVec4x16v16(SVec4x16v16 &sv,
	                            const SVec4x16v16 sv1,
	                            const SVec4x16v16 sv2) {

           
	    sv.s0 = sv1.s0-sv2.s0;
	    sv.s1 = sv1.s1-sv2.s1;
	    sv.s2 = sv1.s2-sv2.s2;
	    sv.s3 = sv1.s3-sv2.s3;
	   
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_sub_SVec4x16v16(SVec4x16v16 &sv,
	                            const SVec4x16v16 sv1){
	                            

           
	    sv.s0 -= sv1.s0;
	    sv.s1 -= sv1.s1;
	    sv.s2 -= sv1.s2;
	    sv.s3 -= sv1.s3;
	   
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_negate(const SVec4x16v16 sv1){
	                           
            SVec4x16v16 sv;
	    const AVX512Vec16 _0 = AVX512Vec16(0.0F);
	    sv.s0 = _0-sv1.s0;
	    sv.s1 = _0-sv1.s1;
	    sv.s2 = _0-sv1.s2;
	    sv.s3 = _0-sv1.s3;
	    return (sv);
       }




       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
	AVX512Vec16
        SVec4x16v16_inner_prod(const SVec4x16v16 sv1,
	                       const SVec4x16v16 sv2) {

	   AVX512Vec16 t0,t1
           t0 = sv1.s0*sv2.s0+sv1.s1*sv2.s1;
	   t1 = sv1.s2*sv2.s2+sv1.s3*sv2.s3;
	   return (t0+t1);
      }


      	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_mul_AVX512Vec16( const SVec4x16v16 sv1
	                             const AVX512Vec16 v) {

          SVec4x16v16 sv;
	  sv.s0 = sv1.s0*v;
	  sv.s1 = sv1.s1*v;
	  sv.s2 = sv1.s2*v;
	  sv.s3 = sv1.s3*v;
	  return (sv);
      }


      	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_mul_AVX512Vec16( SVec4x16v16 &sv,
	                             const SVec4x16v16 sv1
	                             const AVX512Vec16 v) {

         
	  sv.s0 = sv1.s0*v;
	  sv.s1 = sv1.s1*v;
	  sv.s2 = sv1.s2*v;
	  sv.s3 = sv1.s3*v;
	 
      }


       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_div_AVX512Vec16( const SVec4x16v16 sv1
	                             const AVX512Vec16 v) {

          SVec4x16v16 sv;
	  sv.s0 = sv1.s0/v;
	  sv.s1 = sv1.s1/v;
	  sv.s2 = sv1.s2/v;
	  sv.s3 = sv1.s3/v;
	  return (sv);
      }


      	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_div_AVX512Vec16( SVec4x16v16 &sv,
	                             const SVec4x16v16 sv1
	                             const AVX512Vec16 v) {

         
	  sv.s0 = sv1.s0/v;
	  sv.s1 = sv1.s1/v;
	  sv.s2 = sv1.s2/v;
	  sv.s3 = sv1.s3/v;
	 
      }


       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_eta_angle(const SVec4x16v16 sv1) {

           AVX512Vec16 v;
	   const AVX512Vec16 half = AVX512Vec16(0.5F);
	   v = atan2(sv1.s2,sv1.s1)*halfl
	   return (v);
      }


       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_deg_linear_pol(const SVec4x16v16 sv1) {

             AVX512Vec16 v;
	     const AVX512Vec16 t = sv1.s2*sv1.s2+sv1.s1*sv1.s1;
	     v = sqrt(t)/sv1.s0;
	     return (v);
      }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_deg_polarization(const SVec4x16v16 sv1) {

             const AVX512Vec16 t0 = sv1.s2*sv1.s2;
	     const AVX512Vec16 t1 = sv1.s1*sv1.s1;
	     const AVX512Vec16 t2 = sv1.s3*sv1.s3;
	     return ((t0+t1+t2)/sv1.s0);
      }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_deg_circular_pol(const SVec4x16v16 sv1) {

            return (sv1.s3/sv1.s0);
      }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_ellipticity(const SVec4x16v16 sv1) {

            const AVX512Vec16 t0 = sv1.s1*sv1.s1+sv1.s2*sv1.s2;
	    return (sv1.s3/(sqrt(t0)));
      }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        AVX512Vec16
	SVec4x16v16_eccentricity(const SVec4x16v16 sv1) {

            const AVX512Vec16 e = SVec4x16v16_ellipticity(sv1);
	    const AVX512Vec16 one = AVX512Vec16(1.0F);
	    return (sqrt(one-e*e));
      }


      	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_zero() {

           SVec4x16v16 sv;
	   sv.s0 = AVX512Vec16();
	   sv.s1 = AVX512Vec16();
	   sv.s2 = AVX512Vec16();
	   sv.s3 = AVX512Vec16();
	   return(sv);
       }


       	__ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_unit_unpolarized() {

           SVec4x16v16 sv;
	   sv.s0 = AVX512Vec16(1.0F);
	   sv.s1 = AVX512Vec16();
	   sv.s2 = AVX512Vec16();
	   sv.s3 = AVX512Vec16();
	   return(sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_unit_s_polarized() {

           SVec4x16v16 sv;
	   sv.s0 = AVX512Vec16(1.0F);
	   sv.s1 = AVX512Vec16(1.0F);
	   sv.s2 = AVX512Vec16();
	   sv.s3 = AVX512Vec16();
	   return(sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_unit_p_polarized() {

           SVec4x16v16 sv;
	   sv.s0 = AVX512Vec16(1.0F);
	   sv.s1 = AVX512Vec16(-1.0F);
	   sv.s2 = AVX512Vec16();
	   sv.s3 = AVX512Vec16();
	   return(sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_lin_polarized_eta(const AVX512Vec16 eta) {

            SVec4x16v16 sv;
	    const AVX512Vec16 _2eta = eta+eta;
	    sv.s0 = AVX512Vec16(1.0F);
	    sv.s1 = cos(_2eta);
	    sv.s2 = sin(_2eta);
	    sv.s3 = AVX512Vec16(0.0F);
	    return (sv);
       }



        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_right_cir_polarized() {

           SVec4x16v16 sv;
	   sv.s0 = AVX512Vec16(1.0F);
	   sv.s1 = AVX512Vec16();
	   sv.s2 = AVX512Vec16();
	   sv.s3 = AVX512Vec16(-1.0F);
	   return(sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_left_cir_polarized() {

           SVec4x16v16 sv;
	   sv.s0 = AVX512Vec16(1.0F);
	   sv.s1 = AVX512Vec16();
	   sv.s2 = AVX512Vec16();
	   sv.s3 = AVX512Vec16(1.0F);
	   return(sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        SVec4x16v16
	SVec4x16v16_unit_generalized(const AVX512Vec16 eta) {

            SVec4x16v16 sv;
	    const AVX512Vec16 _2eta = eta+eta;
	    sv.s0 = AVX512Vec16(1.0F);
	    sv.s1 = cos(_2eta);
	    sv.s2 = sin(_2eta);
	    sv.s3 = AVX512Vec16(0.0F);
	    return (sv);
       }


        __ATTR_REGCALL__
        __ATTR_ALWAYS_INLINE__
	__ATTR_ALIGN__(32)
	static inline
        void
	SVec4x16v16_unit_generalized( SVec4x16v16 &sv,
	                              const AVX512Vec16 eta) {

           
	    const AVX512Vec16 _2eta = eta+eta;
	    sv.s0 = AVX512Vec16(1.0F);
	    sv.s1 = cos(_2eta);
	    sv.s2 = sin(_2eta);
	    sv.s3 = AVX512Vec16(0.0F);
	   
       }




      

       

	  /*
                           Mueller Matrix based on 16 16-tuple SIMD real types.
                       */
           /*            typedef struct __ATTR_ALIGN__(64)  MMat16x16v16 {

                                  AVX512Vec16 m0;
				  AVX512Vec16 m1;
				  AVX512Vec16 m2;
				  AVX512Vec16 m3;
				  AVX512Vec16 m4;
				  AVX512Vec16 m5;
				  AVX512Vec16 m6;
				  AVX512Vec16 m7;
				  AVX512Vec16 m8;
				  AVX512Vec16 m9;
				  AVX512Vec16 m10;
				  AVX512Vec16 m11;
				  AVX512Vec16 m12;
				  AVX512Vec16 m13;
				  AVX512Vec16 m14;
				  AVX512Vec16 m15;
		       }MMat16x16v16;
              */


         //using Vec16 = AVX512Vec16;
	 
	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         MMat16x16v16
	 MMat16x16v16_set_1() {

             MMat16x16v16 m;
	     m.m0  = AVX512Vec16();
	     m.m1  = AVX512Vec16();
	     m.m2  = AVX512Vec16();
	     m.m3  = AVX512Vec16();
	     m.m4  = AVX512Vec16();
	     m.m5  = AVX512Vec16();
	     m.m6  = AVX512Vec16();
	     m.m7  = AVX512Vec16();
	     m.m8  = AVX512Vec16();
	     m.m9  = AVX512Vec16();
	     m.m10 = AVX512Vec16();
	     m.m11 = AVX512Vec16();
	     m.m12 = AVX512Vec16();
	     m.m13 = AVX512Vec16();
	     m.m14 = AVX512Vec16();
	     m.m15 = AVX512Vec16();
	     return (m);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 MMat16x16v16_set_1(MMat16x16v16 &m) {

             m.m0  = AVX512Vec16();
	     m.m1  = AVX512Vec16();
	     m.m2  = AVX512Vec16();
	     m.m3  = AVX512Vec16();
	     m.m4  = AVX512Vec16();
	     m.m5  = AVX512Vec16();
	     m.m6  = AVX512Vec16();
	     m.m7  = AVX512Vec16();
	     m.m8  = AVX512Vec16();
	     m.m9  = AVX512Vec16();
	     m.m10 = AVX512Vec16();
	     m.m11 = AVX512Vec16();
	     m.m12 = AVX512Vec16();
	     m.m13 = AVX512Vec16();
	     m.m14 = AVX512Vec16();
	     m.m15 = AVX512Vec16();
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         MMat16x16v16
	 MMat16x16v16_set_2(const AVX512Vec16 * __restrict m) {

             MMat16x16v16 mat;
	     mat.m0  = m[0];
	     mat.m1  = m[1];
	     mat.m2  = m[2];
	     mat.m3  = m[3];
	     mat.m4  = m[4];
	     mat.m5  = m[5];
	     mat.m6  = m[6];
	     mat.m7  = m[7];
	     mat.m8  = m[8];
	     mat.m9  = m[9];
	     mat.m10 = m[10];
	     mat.m11 = m[11];
	     mat.m12 = m[12];
	     mat.m13 = m[13];
	     mat.m14 = m[14];
	     mat.m15 = m[15];
	     return (mat);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         MMat16x16v16
	 MMat16x16v16_set_2(const AVX512Vec16 * __restrict m,
	                    MMat16x16v16 &mat) {

             
	     mat.m0  = m[0];
	     mat.m1  = m[1];
	     mat.m2  = m[2];
	     mat.m3  = m[3];
	     mat.m4  = m[4];
	     mat.m5  = m[5];
	     mat.m6  = m[6];
	     mat.m7  = m[7];
	     mat.m8  = m[8];
	     mat.m9  = m[9];
	     mat.m10 = m[10];
	     mat.m11 = m[11];
	     mat.m12 = m[12];
	     mat.m13 = m[13];
	     mat.m14 = m[14];
	     mat.m15 = m[15];
	     
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         MMat16x16v16
	 MMat16x16v16_set_3(const float * __restrict __ATTR_ALIGN__(64) m) {

             MMat16x16v16 mat;
	     mat.m0  = AVX512Vec16(&m[0*16]);
	     mat.m1  = AVX512Vec16(&m[1*16]);
	     mat.m2  = AVX512Vec16(&m[2*16]);
	     mat.m3  = AVX512Vec16(&m[3*16]);
	     mat.m4  = AVX512Vec16(&m[4*16]);
	     mat.m5  = AVX512Vec16(&m[5*16]);
	     mat.m6  = AVX512Vec16(&m[6*16]);
	     mat.m7  = AVX512Vec16(&m[7*16]);
	     mat.m8  = AVX512Vec16(&m[8*16]);
	     mat.m9  = AVX512Vec16(&m[9*16]);
	     mat.m10 = AVX512Vec16(&m[10*16]);
	     mat.m11 = AVX512Vec16(&m[11*16]);
	     mat.m12 = AVX512Vec16(&m[12*16]);
	     mat.m13 = AVX512Vec16(&m[13*16]);
	     mat.m14 = AVX512Vec16(&m[14*16]);
	     mat.m15 = AVX512Vec16(&m[15*16]);
	     return (mat);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         MMat16x16v16
	 MMat16x16v16_set_3( MMat16x16v16 &mat
	                     const float * __restrict __ATTR_ALIGN__(64) m) {

             
	     mat.m0  = AVX512Vec16(&m[0*16]);
	     mat.m1  = AVX512Vec16(&m[1*16]);
	     mat.m2  = AVX512Vec16(&m[2*16]);
	     mat.m3  = AVX512Vec16(&m[3*16]);
	     mat.m4  = AVX512Vec16(&m[4*16]);
	     mat.m5  = AVX512Vec16(&m[5*16]);
	     mat.m6  = AVX512Vec16(&m[6*16]);
	     mat.m7  = AVX512Vec16(&m[7*16]);
	     mat.m8  = AVX512Vec16(&m[8*16]);
	     mat.m9  = AVX512Vec16(&m[9*16]);
	     mat.m10 = AVX512Vec16(&m[10*16]);
	     mat.m11 = AVX512Vec16(&m[11*16]);
	     mat.m12 = AVX512Vec16(&m[12*16]);
	     mat.m13 = AVX512Vec16(&m[13*16]);
	     mat.m14 = AVX512Vec16(&m[14*16]);
	     mat.m15 = AVX512Vec16(&m[15*16]);
	     
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         MMat16x16v16
	 MMat16x16v16_set_4(const MMat16x16v16 m) {

              MMat16x16v16 mat;
	      mat.m0  = m.m0;
	      mat.m1  = m.m1;
	      mat.m2  = m.m2;
	      mat.m3  = m.m3;
	      mat.m4  = m.m4;
	      mat.m5  = m.m5;
	      mat.m6  = m.m6;
	      mat.m7  = m.m7;
	      mat.m8  = m.m8;
	      mat.m9  = m.m9;
	      mat.m10 = m.m10;
	      mat.m11 = m.m11;
	      mat.m12 = m.m12;
	      mat.m13 = m.m13;
	      mat.m14 = m.m14;
	      mat.m15 = m.m15;
	      return (mat);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
         MMat16x16v16
	 MMat16x16v16_set_4( MMat16x16v16 &mat
	                     const MMat16x16v16 m) {

              
	      mat.m0  = m.m0;
	      mat.m1  = m.m1;
	      mat.m2  = m.m2;
	      mat.m3  = m.m3;
	      mat.m4  = m.m4;
	      mat.m5  = m.m5;
	      mat.m6  = m.m6;
	      mat.m7  = m.m7;
	      mat.m8  = m.m8;
	      mat.m9  = m.m9;
	      mat.m10 = m.m10;
	      mat.m11 = m.m11;
	      mat.m12 = m.m12;
	      mat.m13 = m.m13;
	      mat.m14 = m.m14;
	      mat.m15 = m.m15;
	      
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 AVX512Vec16
	 MMat16x16v16_max_transmision(const MMat16x16v16 mat) {

              AVX512Vec16 mt;
	      const AVX512Vec16 t0 = mat.m1*mat.m1;
	      const AVX512Vec16 t1 = mat.m2*mat.m2;
	      const AVX512Vec16 t2 = mat.m3*mat.m3;
	      mt = mat.m0+sqrt(t0+t1+t2);
	      return (mt);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 AVX512Vec16
	 MMat16x16v16_min_transmision(const MMat16x16v16 mat) {

              AVX512Vec16 mt;
	      const AVX512Vec16 t0 = mat.m1*mat.m1;
	      const AVX512Vec16 t1 = mat.m2*mat.m2;
	      const AVX512Vec16 t2 = mat.m3*mat.m3;
	      mt = mat.m0-sqrt(t0+t1+t2);
	      return (mt);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 AVX512Vec16
	 MMat16x16v16_diattenuation(const MMat16x16v16 mat) {

              AVX512Vec16 mt;
	      const AVX512Vec16 t0 = mat.m1*mat.m1;
	      const AVX512Vec16 t1 = mat.m2*mat.m2;
	      const AVX512Vec16 t2 = mat.m3*mat.m3;
	      mt = sqrt(t0+t1+t2)/mat.m0;
	      return (mt);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 AVX512Vec16
	 MMat16x16v16_linear_diattenuation(const MMat16x16v16 mat) {
	
              AVX512Vec16 mt;
	      const AVX512Vec16 t0 = mat.m1*mat.m1;
	      const AVX512Vec16 t1 = mat.m2*mat.m2;
	      mt = sqrt(t0+t1)/mat.m0;
	      return (mt);
       }


       	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 AVX512Vec16
	 MMat16x16v16_polarization_loss(const MMat16x16v16 mat) {

	      const AVX512Vec16 invlog10(0.00000000000000000000000000000001F);
	      const AVX512Vec16 _10(10.0F);
              AVX512Vec16 mt;
	      const AVX512Vec16 t0 = MMat16x16v16_max_transmission(mat)/
	                             MMat16x16v16_min_transmission(mat);
	      
	      mt = _10*t0*invlog10;
	      return (mt);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 AVX512Vec16
	 MMat16x16v16_polarizance(const MMat16x16v16 mat) {

              AVX512Vec16 mt;
	      const AVX512Vec16 t0 = mat.m5*mat.m5;
	      const AVX512Vec16 t1 = mat.m9.mat.m9;
	      const AVX512Vec16 t2 = mat.m13*mat.m13;
	      mt = sqrt(t0+t1+t2)/mat.m0;
	      return (mt);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 AVX512Vec16
	 MMat16x16v16_extinct_ratio(const MMat16x16v16 mat) {

            return (MMat16x16v16_max_transmission(mat)/
	            MMat16x16v16_min_transmission(mat));
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_unit(const AVX512Vec16 transmittance) {

             MMat16x16v16 mat;
	     mat.m0  = transmittance;
	     mat.m5  = mat.m0;
	     mat.m10 = mat.m0;
	     mat.m15 = mat.m0;
	     return (mat);
	}


	 __ATTR_REGCALL__
	 __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_depolarization(const AVX512Vec16 transmittance,
	                             const AVX512Vec16 depolarization) {

             MMat16x16v16 mat;
	     const AVX512Vec16 one(1.0F);
	     mat.m0  = transmittance;
	     mat.m5  = transmittance*(one-depolarization);
	     mat.m10 = mat.m5;
	     mat.m15 = mat.m5;
	     return (mat);
       }


       	 __ATTR_REGCALL__
	 __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_diagonal(const AVX512Vec16 diag1,
	                       const AVX512Vec16 diag2,
			       const AVX512Vec16 diag3,
			       const AVX512Vec16 diag4) {

             MMat16x16v16 mat;
	     mat.m0  = diag1;
	     mat.m5  = diag2;
	     mat.m10 = diag3;
	     mat.m15 = diag4;
	     return (mat);
      }


      	 __ATTR_REGCALL__
	 __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_linear_polarizer(const AVX512Vec16 maxt,
	                               const AVX512Vec16 mint,
				       const AVX512Vec16 ang) {

            MMat16x16v16 mat;
	    const AVX512Vec16 _2ang = ang+ang;
	    const AVX512Vec16 _0_5(0.5F);
	    const AVX512Vec16 t0 = _0_5*(maxt+mint);
	    mat.m0 = t0;
	    const AVX512Vec16 t1 = _0_5*(maxt-mint);
	    const AVX512Vec16 t2 = sqrt(maxt*mint);
	    mat.m15 = t2;
	    const AVX512Vec16 t3  = cos(_2ang);
	    const AVX512Vec16 st3 = sqrt(t3);
	    const AVX512Vec16 t4  = sin(_2ang);
	    const AVX512Vec16 st4 = sqrt(t4);
	    mat.m1  = t1*t3;
	    mat.m2  = t1*t4;
	    mat.m4  = mat.m1;
	    mat.m5  = t0*st3+t2*st4;
	    mat.m6  = (t0-t2)*t3*t4;
	    mat.m8  = mat.m2;
	    mat.m9  = mat.m6;
	    mat.m10 = t2*st3+t0*st4;
	    return (mat);
      }


       


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_mul_AVX512Vec16(const MMat16x16v16 mat,
	                              const AVX512Vec16 v) {

              MMat16x16v16 res;
	      res.m0 = mat.m0*v;
	      res.m1 = mat.m1*v;
	      res.m2 = mat.m2*v;
	      res.m3 = mat.m3*v;
	      res.m4 = mat.m4*v;
	      res.m5 = mat.m5*v;
	      res.m6 = mat.m6*v;
	      res.m7 = mat.m7*v;
	      res.m8 = mat.m8*v;
	      res.m9 = mat.m9*v;
	      res.m10 = mat.m10*v;
	      res.m11 = mat.m11*v;
	      res.m12 = mat.m12*v;
	      res.m13 = mat.m13*v;
	      res.m14 = mat.m14*v;
	      res.m15 = mat.m15*v;
	      return (res);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 MMat16x16v16_mul_AVX512Vec16(const MMat16x16v16 mat,
	                              MMat16x16v16 &res,
	                              const AVX512Vec16 v) {

             
	      res.m0 = mat.m0*v;
	      res.m1 = mat.m1*v;
	      res.m2 = mat.m2*v;
	      res.m3 = mat.m3*v;
	      res.m4 = mat.m4*v;
	      res.m5 = mat.m5*v;
	      res.m6 = mat.m6*v;
	      res.m7 = mat.m7*v;
	      res.m8 = mat.m8*v;
	      res.m9 = mat.m9*v;
	      res.m10 = mat.m10*v;
	      res.m11 = mat.m11*v;
	      res.m12 = mat.m12*v;
	      res.m13 = mat.m13*v;
	      res.m14 = mat.m14*v;
	      res.m15 = mat.m15*v;
	      
	}

	
	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_div_AVX512Vec16(const MMat16x16v16 mat,
	                              const AVX512Vec16 v) {

            MMat16x16v16 res;
	    res.m0  = mat.m0/v;
	    res.m1  = mat.m1/v;
	    res.m2  = mat.m2/v;
	    res.m3  = mat.m3/v;
	    res.m4  = mat.m4/v;
	    res.m5  = mat.m5/v;
	    res.m6  = mat.m6/v;
	    res.m7  = mat.m7/v;
	    res.m8  = mat.m8/v;
	    res.m9  = mat.m9/v;
	    res.m10 = mat.m10/v;
	    res.m11 = mat.m11/v;
	    res.m12 = mat.m12/v;
	    res.m13 = mat.m13/v;
	    res.m14 = mat.m14/v;
	    res.m15 = mat.m15/v;
	    return (res);
       }


       	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 MMat16x16v16_div_AVX512Vec16(const MMat16x16v16 mat
	                              MMat16x16v16 &res,
	                              const AVX512Vec16 v) {

            res.m0  = mat.m0/v;
	    res.m1  = mat.m1/v;
	    res.m2  = mat.m2/v;
	    res.m3  = mat.m3/v;
	    res.m4  = mat.m4/v;
	    res.m5  = mat.m5/v;
	    res.m6  = mat.m6/v;
	    res.m7  = mat.m7/v;
	    res.m8  = mat.m8/v;
	    res.m9  = mat.m9/v;
	    res.m10 = mat.m10/v;
	    res.m11 = mat.m11/v;
	    res.m12 = mat.m12/v;
	    res.m13 = mat.m13/v;
	    res.m14 = mat.m14/v;
	    res.m15 = mat.m15/v;
       }


       
         __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_add_MMat16x16v16(const MMat16x16v16 mat1,
	                               const MMat16x16v16 mat2) {

             MMat16x16v16 mat;
	     mat.m0  = mat1.m0+mat2.m0;
	     mat.m1  = mat1.m1+mat2.m1;
	     mat.m2  = mat1.m2+mat2.m2;
	     mat.m3  = mat1.m3+mat2.m3;
	     mat.m4  = mat1.m4+mat2.m4;
	     mat.m5  = mat1.m5+mat2.m5;
	     mat.m6  = mat1.m6+mat2.m6;
	     mat.m7  = mat1.m7+mat2.m7;
	     mat.m8  = mat1.m8+mat2.m8;
	     mat.m9  = mat1.m9+mat2.m9;
	     mat.m10 = mat1.m10+mat2.m10;
	     mat.m11 = mat1.m11+mat2.m11;
	     mat.m12 = mat1.m12+mat2.m12;
	     mat.m13 = mat1.m13+mat2.m13;
	     mat.m14 = mat1.m14+mat2.m14;
	     mat.m15 = mat1.m15+mat2.m15;
	     return (mat);
        }


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_add_AVX512Vec16(const MMat16x16v16 m,
	                              const AVX512Vec16 v) {

             MMat16x16v16 mat;
	     mat.m0  = m.m0+v;
	     mat.m1  = m.m1+v;
	     mat.m2  = m.m2+v;
	     mat.m3  = m.m3+v;
	     mat.m4  = m.m4+v;
	     mat.m5  = m.m5+v;
	     mat.m6  = m.m6+v;
	     mat.m7  = m.m7+v;
	     mat.m8  = m.m8+v;
	     mat.m9  = m.m9+v;
	     mat.m10 = m.m10+v;
	     mat.m11 = m.m11+v;
	     mat.m12 = m.m12+v;
	     mat.m13 = m.m13+v;
	     mat.m14 = m.m14+v;
	     mat.m15 = m.m15+v;
	     return (mat);
        }


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 MMat16x16v16_add_MMat16x16v16(MMat16x16v16 &mat,
	                               const MMat16x16v16 mat1,
	                               const MMat16x16v16 mat2) {
             mat.m0  = mat1.m0+mat2.m0;
	     mat.m1  = mat1.m1+mat2.m1;
	     mat.m2  = mat1.m2+mat2.m2;
	     mat.m3  = mat1.m3+mat2.m3;
	     mat.m4  = mat1.m4+mat2.m4;
	     mat.m5  = mat1.m5+mat2.m5;
	     mat.m6  = mat1.m6+mat2.m6;
	     mat.m7  = mat1.m7+mat2.m7;
	     mat.m8  = mat1.m8+mat2.m8;
	     mat.m9  = mat1.m9+mat2.m9;
	     mat.m10 = mat1.m10+mat2.m10;
	     mat.m11 = mat1.m11+mat2.m11;
	     mat.m12 = mat1.m12+mat2.m12;
	     mat.m13 = mat1.m13+mat2.m13;
	     mat.m14 = mat1.m14+mat2.m14;
	     mat.m15 = mat1.m15+mat2.m15;
	}


	


         __ATTR_REGCALL__
	 __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_sub_MMat16x16v16(const MMat16x16v16 mat1,
	                               const MMat16x16v16 mat2) {

             MMat16x16v16 mat;
	     mat.m0  = mat1.m0-mat2.m0;
	     mat.m1  = mat1.m1-mat2.m1;
	     mat.m2  = mat1.m2-mat2.m2;
	     mat.m3  = mat1.m3-mat2.m3;
	     mat.m4  = mat1.m4-mat2.m4;
	     mat.m5  = mat1.m5-mat2.m5;
	     mat.m6  = mat1.m6-mat2.m6;
	     mat.m7  = mat1.m7-mat2.m7;
	     mat.m8  = mat1.m8-mat2.m8;
	     mat.m9  = mat1.m9-mat2.m9;
	     mat.m10 = mat1.m10-mat2.m10;
	     mat.m11 = mat1.m11-mat2.m11;
	     mat.m12 = mat1.m12-mat2.m12;
	     mat.m13 = mat1.m13-mat2.m13;
	     mat.m14 = mat1.m14-mat2.m14;
	     mat.m15 = mat1.m15-mat2.m15;
	     return (mat);
        }


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 MMat16x16v16_sub_MMat16x16v16(MMat16x16v16 &mat,
	                               const MMat16x16v16 mat1,
	                               const MMat16x16v16 mat2) {
             mat.m0  = mat1.m0-mat2.m0;
	     mat.m1  = mat1.m1-mat2.m1;
	     mat.m2  = mat1.m2-mat2.m2;
	     mat.m3  = mat1.m3-mat2.m3;
	     mat.m4  = mat1.m4-mat2.m4;
	     mat.m5  = mat1.m5-mat2.m5;
	     mat.m6  = mat1.m6-mat2.m6;
	     mat.m7  = mat1.m7-mat2.m7;
	     mat.m8  = mat1.m8-mat2.m8;
	     mat.m9  = mat1.m9-mat2.m9;
	     mat.m10 = mat1.m10-mat2.m10;
	     mat.m11 = mat1.m11-mat2.m11;
	     mat.m12 = mat1.m12-mat2.m12;
	     mat.m13 = mat1.m13-mat2.m13;
	     mat.m14 = mat1.m14-mat2.m14;
	     mat.m15 = mat1.m15-mat2.m15;
	}


	 __ATTR_REGCALL__
	 __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_negate(const MMat16x16v16 m) {

              const AVX512Vec16 zero();
	      MMat16x16v16 mat;
	      mat.m0  = zero-m.m0;
	      mat.m1  = zero-m.m1;
	      mat.m2  = zero-m.m2;
	      mat.m3  = zero-m.m3;
	      mat.m4  = zero-m.m4;
	      mat.m5  = zero-m.m5;
	      mat.m6  = zero-m.m6;
	      mat.m7  = zero-m.m7;
	      mat.m8  = zero-m.m8;
	      mat.m9  = zero-m.m9;
	      mat.m10 = zero-m.m10;
	      mat.m11 = zero-m.m11;
	      mat.m12 = zero-m.m12;
	      mat.m13 = zero-m.m13;
	      mat.m14 = zero-m.m14;
	      mat.m15 = zero-m.m15;
	      return (mat);
	}




	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 MMat16x16v16
	 MMat16x16v16_mul_MMat16x16v16(const MMat16x16v16 mat1,
	                               const MMat16x16v16 mat2) {

              MMat16x16v16 out;
	      AVX512Vec16 c0();
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m0.m_v16,c0.m_v16);
	      out.m0 = c0;
	      AVX512Vec16 c1();
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m1.m_v16,c1.m_v16);
	      out.m1 = c1;
	      AVX512Vec16 c2()
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m2.m_v16,c2.m_v16);
	      out.m2 = c2;
	      AVX512Vec16 c3();
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m3.m_v16,c3.m_v16);
	      out.m3 = c3;
	      AVX512Vec16 c4();
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m4.m_v16,c4.m_v16);
	      out.m4 = c4;
	      AVX512Vec16 c5();
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m5.m_v16,c5.m_v16);
	      out.m5 = c5;
	      AVX512Vec16 c6();
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m6.m_v16,c6.m_v16);
	      out.m6 = c6;
	      AVX512Vec16 c7();
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m7.m_v16,c7.m_v16);
	      out.m7 = c7;
	      AVX512Vec16 c8();
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m8.m_v16,c8.m_v16);
	      out.m = c8;
	      AVX512Vec16 c9();
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m9.m_v16,c9.m_v16);
	      out.m9 = c9;
	      AVX512Vec16 c10();
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m10.m_v16,c10.m_v16);
	      out.m10 = c10;
	      AVX512Vec16 c11();
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m11.m_v16,c11.m_v16);
	      out.m11 = c11;
	      AVX512Vec16 c12();
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m12.m_v16,c12.m_v16);
	      out.m12 = c12;
	      AVX512Vec16 c13();
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m13.m_v16,c13.m_v16);
	      out.m13 = c13;
	      AVX512Vec16 c14();
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m14.m_v16,c14.m_v16);
	      out.m14 = c14;
	      AVX512Vec16 c15();
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m15.m_v16,c15.m_v16);
	      out.m15 = c15;
	      return (out);
	}


	 __ATTR_REGCALL__
         __ATTR_ALWAYS_INLINE__
	 __ATTR_ALIGN__(32)
	 static inline
	 void
	 MMat16x16v16_mul_MMat16x16v16(const MMat16x16v16 mat1,
	                               const MMat16x16v16 mat2,
				       MMat16x16v16 &out) {

              AVX512Vec16 c0();
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m0.m_v16,c0.m_v16);
	      c0.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m0.m_v16,c0.m_v16);
	      out.m0 = c0;
	      AVX512Vec16 c1();
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m1.m_v16,c1.m_v16);
	      c1.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m1.m_v16,c1.m_v16);
	      out.m1 = c1;
	      AVX512Vec16 c2()
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m2.m_v16,c2.m_v16);
	      c2.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m2.m_v16,c2.m_v16);
	      out.m2 = c2;
	      AVX512Vec16 c3();
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m3.m_v16,c3.m_v16);
	      c3.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m3.m_v16,c3.m_v16);
	      out.m3 = c3;
	      AVX512Vec16 c4();
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m4.m_v16,c4.m_v16);
	      c4.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m4.m_v16,c4.m_v16);
	      out.m4 = c4;
	      AVX512Vec16 c5();
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m5.m_v16,c5.m_v16);
	      c5.m_v16 = _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m5.m_v16,c5.m_v16);
	      out.m5 = c5;
	      AVX512Vec16 c6();
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m6.m_v16,c6.m_v16);
	      c6.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m6.m_v16,c6.m_v16);
	      out.m6 = c6;
	      AVX512Vec16 c7();
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m7.m_v16,c7.m_v16);
	      c7.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m7.m_v16,c7.m_v16);
	      out.m7 = c7;
	      AVX512Vec16 c8();
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m8.m_v16,c8.m_v16);
	      c8.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m8.m_v16,c8.m_v16);
	      out.m = c8;
	      AVX512Vec16 c9();
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m9.m_v16,c9.m_v16);
	      c9.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m9.m_v16,c9.m_v16);
	      out.m9 = c9;
	      AVX512Vec16 c10();
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m10.m_v16,c10.m_v16);
	      c10.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m10.m_v16,c10.m_v16);
	      out.m10 = c10;
	      AVX512Vec16 c11();
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m11.m_v16,c11.m_v16);
	      c11.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m11.m_v16,c11.m_v16);
	      out.m11 = c11;
	      AVX512Vec16 c12();
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m12.m_v16,c12.m_v16);
	      c12.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m12.m_v16,c12.m_v16);
	      out.m12 = c12;
	      AVX512Vec16 c13();
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m13.m_v16,c13.m_v16);
	      c13.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m13.m_v16,c13.m_v16);
	      out.m13 = c13;
	      AVX512Vec16 c14();
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m14.m_v16,c14.m_v16);
	      c14.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m14.m_v16,c14.m_v16);
	      out.m14 = c14;
	      AVX512Vec16 c15();
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m0.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m1.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m2.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m3.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m4.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m5.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m6.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m7.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m8.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m9.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m10.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m11.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m12.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m13.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m14.m_v16,mat2.m15.m_v16,c15.m_v16);
	      c15.m_v16 =  _mm512_fmadd_ps(mat1.m15.m_v16,mat2.m15.m_v16,c15.m_v16);
	      out.m15 = c15;
	}




	
	


	


	

	


	

	


	


	


	


	

			   
    } // math



} // gms














#endif /*__GMS_MUELLER_CALCULUS_AVX512_HPP__*/
