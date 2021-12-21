

#ifndef __GMS_JONES_VECTOR_AVX512_HPP__
#define __GMS_JONES_VECTOR_AVX512_HPP__ 211220211503



namespace file_info {

 const unsigned int gGMS_JONES_VECTOR_AVX512_MAJOR = 1U;
 const unsigned int gGMS_JONES_VECTOR_AVX512_MINOR = 1U;
 const unsigned int gGMS_JONES_VECTOR_AVX512_MICRO = 0U;
 const unsigned int gGMS_JONES_VECTOR_AVX512_FULLVER =
  1000U*gGMS_JONES_VECTOR_AVX512_MAJOR+100U*gGMS_JONES_VECTOR_AVX512_MINOR+10U*gGMS_JONES_VECTOR_AVX512_MICRO;
 const char * const pgGMS_JONES_VECTOR_AVX512_CREATION_DATE = "21-12-2021 15:03 +00200 (TUE 21 DEC 2021 15:03 GMT+2)";
 const char * const pgGMS_JONES_VECTOR_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_JONES_VECTOR_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_JONES_VECTOR_AVX512_SYNOPSIS      = "AVX512 based Jones-Vector implementation based on SIMD complex 16-tuple type."


}

#include <immintrin.h>
#include <complex>
#include "GMS_avx512c16f32.h"
#include "GMS_config.h"


namespace  gms {

         namespace math {


                     /*
                           SIMD implementation of Jones Vector (calculus) for 16 points.
                       */
                        struct JVec2x16c16 __ATTR_ALIGN__(64) {


                                  ZMM16c4 j0;
				  ZMM16c4 j1;

				  JVec2x16c16() {}

				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  JVec2x16c16(const ZMM16c4 s,
				              const ZMM16c4 p) {

                                      j0 = s;
				      j1 = p;
				}

				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  JVec2x16c16(const JVec2x16c16 &x) {

				      j0 = x.j0;
				      j1 = x.j1;
				 }

				 __ATTR_ALWAYS_INLINE__
				 __ATTR_ALIGN__(32)
		                 JVec2x16c16(const std::complex<float> s,
				             const std::complex<float> p) {

                                       j0 = ZMM16c4{s};
				       j1 = ZMM16c4{p};
				 }

				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  JVec2x16c16(const __m512 re1,
				              const __m512 im1,
					      const __m512 re2,
					      const __m512 im2) {

                                       j0 = ZMM16c4{re1,im1};
				       j1 = ZMM16c4{re2,im2};
				 }


				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  JVec2x16c16(const float * __restrict __ATTR_ALIGN__(64) re1,
				              const float * __restrict __ATTR_ALIGN__(64) im1,
					      const float * __restrict __ATTR_ALIGN__(64) re2,
					      const float * __restrict __ATTR_ALIGN__(64) im2) {

                                        j0 = ZMM16c4{re1,im1};
					j1 = ZMM16c4{re2,im2};
				 }


				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  JVec2x16c16 &
				  operator=(const JVec2x16c16 &x) {
				      if(this==&x) {return (*this);}
                                      j0 = x.j0;
				      j1 = x.j1;
				      return (*this);
				 }

				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  JVec2x16c16
				  operator*(const ZMM16c16 x) {

                                      return JVec2x16c16{j0*x,j1*x};
				}

				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  friend JVec2x16c16
				  operator*(const JVec2x16c16 x,
				            const ZMM16c4 y) {

				      return(x*y);
				}

				  __ATTR_ALWAYS_INLINE__
				  __ATTR_ALIGN__(32)
		                  __ATTR_VECTORCALL__
				  JVec2x16c16 &
				  operator*=(const ZMM16c4 x) {

                                     *this = x*(*this);
				}
		      };


		      

	 
       } // math


} // gms























#endif /* __GMS_JONES_VECTOR_AVX512_HPP__*/
