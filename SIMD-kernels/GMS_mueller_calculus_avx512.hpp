

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
   }



}














#endif /*__GMS_MUELLER_CALCULUS_AVX512_HPP__*/
