
#ifndef __GMS_AVXC8F32_H__
#define __GMS_AVXC8F32_H__



namespace file_info {

     const unsigned int gGMS_AVXC8F32_MAJOR = 1U;
     const unsigned int gGMS_AVXC8F32_MINOR = 0U;
     const unsigned int gGMS_AVXC8F32_MICRO = 0U;
     const unsigned int gGMS_AVXC8F32_FULLVER =
           1000U*gGMS_AVXC8F32_MAJOR+100U*gGMS_AVXC8F32_MINOR+10U*gGMS_AVXC8F32_MICRO;
     const char * const pgGMS_AVXC8F32_CREATE_DATE = "30-12-2019 12:07 +00200 (MON 30 DEC 2019 GMT+2)";
     const char * const pgGMS_AVXC8F32_BUILD_DATE  = __DATE__ " " __TIME__;
     const char * const pgGMS_AVXC8F32_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
     const char * const pgGMS_AVXC8F32_SYNOPSIS    = "AVX complex number class decomposed into real and imaginary parts stored as 8-tuple."
}

#include <cstdint>
#include <immintrin.h>
#include <complex>
#include "GMS_config.h"

namespace gms {

       namespace math {

#if !defined(AVXC8F32_SETPS)
     #define AVXC8F32_SETPS(x)  _mm256_set_ps(1.0f,1.0f,1.0f,1.0f,   \
                                              1.0f,1.0f,1.0f,(x));
#endif


                  struct AVXc8f32_t {

                         __m256 m_re;
			 __m256 m_im;
			 
                         // Set to 0.0f
 __ATTR_HOT__ __ATTR_ALIGN__(16) AVXc8f32_t() {
                                 m_re = _mm256_setzero_ps();
			         m_im = _mm256_setzero_ps();
			 }

__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         AVXc8f32(const float * __restrict __ATTR_ALIGN__(32) re,
                                                  const float * __restrict __ATTR_ALIGN__(32) im) {
                                       #if defined __GNUC__ && !defined __INTEL_COMPILER
				                   re = (const float*)__builtin_assume_aligned(re,32);
						   im = (const float*)__builtin_assume_aligned(im,32);
				       #elif defined __ICC || defined __INTEL_COMPILER
				                    __assume_aligned(re,32);
						    __assume_aligned(im,32);
				       #endif
				                    m_re = _mm256_load_ps(&re[0]);
						    m_im = _mm256_load_ps(&im[0]);
				       }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          AVXc8f32(const float re,
						   const float im) {

						m_re = _mm256_set1_ps(re);
						m_im = _mm256_set1_ps(im);
				       }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         AVXc8f32(const std::complex<float> c) {

                                                m_re = _mm256_set1_ps(c.real());
						m_im = _mm256_set1_ps(c.imag());
				       }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         AVXc8f32(const float re) {

					        m_re = _mm256_set1_ps(re);
						m_im = _mm256_setzero_ps();
				       }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32(const float re0,
					         const float re1,
						 const float re2,
						 const float re3,
						 const float re4,
						 const float re5,
						 const float re6,
						 const float re7) {

						 m_re = _mm256_setr_ps(re0,re1,re2,re3,
						                       re4,re5,re6,re7);
						 m_im = _mm256_setzero_ps();
					}
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32(const float re0,
					         const float re1,
						 const float re2,
						 const float re3,
						 const float re4,
						 const float re5,
						 const float re6,
						 const float re7,
						 const float im0,
						 const float im1,
						 const float im2,
						 const float im3,
						 const float im4,
						 const float im5,
						 const float im6,
						 const float im7) {

						 m_re = _mm256_setr_ps(re0,re1,re2,re3,
						                       re4,re5,re6,re7);
						 m_im = _mm256_setr_ps(im0,im1,im2,im3,
						                       im4,im5,im6,im7);
					}
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32(const __m256 re,
					         const __m256 im) {

                                                 m_re = re;
						 m_im = im;
					}
__ATTR_HOT__ __ATTR_ALIGN__(16)
                                        AVXc8f32(const AVXc8f32 &x) {

                                                 m_re = x.m_re;
						 m_im = x.m_im;
					}
                                        ~AVXc8f32() {

					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 & load_a(const float * __restrict __ATTR_ALIGN__(32) re,
					                  const float * __restrict __ATTR_ALIGN__(32) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER

                                                 re = (const float*)__builtin_assume_aligned(re,32);
						 im = (const float*)__builtin_assume_aligned(im,32);
#elif defined __ICC || defined __INTEL_COMPILER
                                                 __assume_aligned(re,32);
						 __assume_aligned(im,32);
#endif
                                                 m_re = _mm256_load_ps(&re[0]);
						 m_im = _mm256_load_ps(&im[0]);
						 return (*this);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 & load_u(const float * __restrict re,
					                  const float * __restrict im) {

						  m_re = _mm256_loadu_ps(&re[0]);
						  m_im = _mm256_loadu_ps(&im[0]);
						  return (*this);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        void store_a(float * __restrict __ATTR_ALIGN__(32) re,
					             float * __restrict __ATTR_ALIGN__(32) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                                                  re = (float*)__builtin_assume_aligned(re,32);
						  im = (float*)__builtin_assume_aligned(im,32);
#elif defined __ICC || defined __INTEL_COMPILER
                                                  __assume_aligned(re,32);
						  __assume_aligned(im,32);
#endif
                                                  _mm256_store_ps(&re[0],m_re);
						  _mm256_store_ps(&im[0],m_im);
				        }
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        void store_u(float * __restrict re,
					             float * __restrict im) {

						    _mm256_storeu_ps(&re[0],m_re);
						    _mm256_storeu_ps(&im[0],m_im);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        void store_nt(float * __restrict re,
					              float * __restrict im) {

						    _mm256_stream_ps(&re[0],m_re);
						    _mm256_stream_ps(&im[0],m_im);
						    _mm_sfence();
					}
__ATTR_COLD__ __ATTR_ALIGN__(16)        void extract_1c(float &s0,
                                                        float &s1,
							const int32_t pos) {

						    float re[8] __ATTR_ALIGN__(32) = {};
						    float im[8] __ATTR_ALIGN__(32) = {};
						    store_a(re);
						    s0 = re[pos&7];
						    store_a(im);
						    s1 = im[pos&7];
				        }
__ATTR_COLD__ __ATTR_ALIGN__(16)        AVXc8f32 & insert_1c(const float s0,
                                                             const float s1,
							     const int32_t pos) {

						     float re[8] __ATTR_ALIGN__(32) = {};
						     float im[8] __ATTR_ALIGN__(32) = {};
						     store_a(&re[0]);
						     re[pos&7] = s0;
						     load_a(&re[0]);
						     store_a(&im[0]);
						     im[pos&7] = s1;
						     load_a(&im[0]);
						     return (*this);
				        }
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        void concatenate_a(float * __restrict __ATTR_ALIGN__(32) out) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                                                     out = (float*)__builtin_assume_aligned(out,32);
#elif defined __ICC || defined __INTEL_COMPILER
                                                     __assume_aligned(out,32);
#endif
                                                     store_a(&out[0],&out[4]);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        void concatenate_u(float * __restrict  out) {

                                                     store_u(&out[0],&out[4]);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 &  maskload(float const * __restrict re,
                                                             float const * __restrict im,
							     const __m256i mask1,
							     const __m256i mask2) {

						       m_re = _mm256_maskload_ps(&re[0],mask1);
						       m_im = _mm256_maskload_ps(&im[0],mask2);
						       return (*this);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        void maskstore(float * __restrict re,
					               float * __restrict im,
						       const __m256i mask1,
						       const __m256i mask2) {

						       _mm256_maskstore_ps(&re[0],mask1,m_re);
						       _mm256_maskstore_ps(&im[0],mask2,m_im);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 & maskload(const AVXc8f32 c8,
					                    const __m256i mask1,
							    const __m256i mask2) {

						        m_re = _mm256_maskload_ps(&c8.m_re[0],mask1);
							m_im = _mm256_maskload_ps(&c8.m_im[0],mask2);
							return (*this);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        void maskstore(AVXc8f32 & c8,
					               const __m256i mask1,
						       const __m256i mask2) {

						       _mm256_maskstore_ps((float*)&c8.m_re[0],mask1,m_re);
						       _mm256_maskstore_ps((float*)&c8.m_im[0],mask1,m_im);
					}
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 & blendv(const AVXc8f32 x
					                 const __m256 maskx,
							 const AVXc8f32 y,
							 const __m256 masky) {
							 
                                              m_re = _mm256_blendv_ps(x.m_re,y.m_re,maskx);
					      m_im = _mm256_blendv_ps(x.m_im,y.m_im,masky);
					      return (*this);
					}
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 & blend(const AVXc8f32 x,
					                 const int32_t imm8x,
							 const AVXc8f32 y,
							 const int32_t imm8y) {

					       m_re = _mm256_blendv_ps(x.m_re,y.m_re,imm8x);
					       m_im = _mm256_blendv_ps(x.m_im,y.m_im,imm8y);
					       return (*this);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16)
                                        __m128 re_lo() {
					
                                               return (_mm256_extractf128_ps(m_re,0));
					}
__ATTR_COLD__ __ATTR_ALIGN__(16)        __m128 re_hi() {
                                               return (_mm256_extractf128_ps(m_re,1));
                                        }

	       } __ATTR_ALIGN__(64);
      
    } // math   

}  // gms


#endif /*__GMS_AVXC8F32_H__*/
