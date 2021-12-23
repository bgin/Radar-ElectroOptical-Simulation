

#ifndef __GMS_MUELLER_CALCULUS_AVX512_HPP__
#define __GMS_MUELLER_CALCULUS_AVX512_HPP__



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

			   
    } // math



} // gms














#endif /*__GMS_MUELLER_CALCULUS_AVX512_HPP__*/
