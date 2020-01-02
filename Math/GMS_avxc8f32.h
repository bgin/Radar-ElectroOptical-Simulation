
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
                                        __m128 re_lo() const {
					
                                               return (_mm256_extractf128_ps(m_re,0));
					}
__ATTR_COLD__ __ATTR_ALIGN__(16)        __m128 re_hi() {
                                               return (_mm256_extractf128_ps(m_re,1));
                                        }
__ATTR_COLD__ __ATTR_ALIGN__(16)        __m128 im_lo() const {
                                               return (_mm256_extractf128_ps(m_im,0));
                                        }
__ATTR_COLD__ __ATTR_ALIGN__(16)        __m128 im_hi() const {
                                               return (_mm256_extractf128_ps(m_im,1));
                                        }
__ATTR_COLD__ __ATTR_ALIGN__(16) 
                                        AVXc8f32 & permute(const int32_t imm8x,
					                   const int32_t imm8y) {

					       m_re = _mm256_permute_ps(m_re,imm8x);
					       m_im = _mm256_permute_ps(m_im,imm8y);
					       return (*this);
					}
__ATTR_COLD__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 & permute(const AVXc8f32 x,
					                   const AVXc8f32 y,
							   const int32_t imm8x,
							   const int32_t imm8y) {

				               m_re = _mm256_permute_ps(x.m_re,imm8x);
					       m_im = _mm256_permute_ps(x.m_im,imm8y);
					       return (*this);
					}
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                        AVXc8f32 & operator(const AVXc8f32 x) {
                                                   if(__builtin_expect(x == &this,0)) return (*this);
						   m_re = x.m_re;
						   m_im = x.m_im;
						   return (*this);
					}	    
					

	       } __ATTR_ALIGN__(64);

__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
					 static inline
					 AVXc8f32 conj(const AVXc8f32 x) {

					         auto tmp = ~x;
						 return (AVXc8f32{x.m_re,tmp.im});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
                                         AVXc8f32 polar(const __m256 ro,
					                const __m256 theta) {

					          const __m256 re_part =
						        _mm256_mul_ps(rho,_mm256_cos_ps(theta));
						  const __m256 im_part =
						        _mm256_mul_ps(rho,_mm256_sin_ps(theta));
						  return (AVXc8f32{re_part,im_part});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 __m256 carg(const AVXc8f32 x) {

					         return (_mm256_atan2_ps(x.m_re,x.m_im));
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16)
                                         static inline
					 _mm256 carg(const float re,
					             const float im) {

						 const __m256 real = AVXC8F32_SETPS(re);
						 const __m256 imag = AVXC8F32_SETPS(im);
						 return (_mm256_atan2_ps(real,imag));
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 csin(const AVXc8f32 x) {

					         const __m256 re_part =
						       _mm256_mul_ps(_mm256_sin_ps(x.m_re),
						                     _mm256_cosh_ps(x.m_im));
						 const __m256 im_part =
						       _mm256_mul_ps(_mm256_cos_ps(x.m_re),
						                     _mm256_sinh_ps(x.m_im));
						 return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 csin(const std::complex<float> x) {

					         const __m256 real = AVXC8F32_SETPS(x.real());
						 const __m256 imag = AVXC8F32_SETPS(x.imag());
						 const __m256 re_part =
						       _mm256_mul_ps(_mm256_sin_ps(real),
						                     _mm256_cosh_ps(imag));
						 const __m256 im_part =
						       _mm256_mul_ps(_mm256_cos_ps(real),
						                     _mm256_sinh_ps(imag));
						 return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 csin(const float re,
					               const float im) {

						 const __m256 real = AVXC8F32_SETPS(re);
						 const __m256 imag = AVXC8F32_SETPS(im);
						 const __m256 re_part =
						       _mm256_mul_ps(_mm256_sin_ps(real),
						                     _mm256_cosh_ps(imag));
						 const __m256 im_part =
						       _mm256_mul_ps(_mm256_cos_ps(real),
						                     _mm256_sinh_ps(imag));
						 return (AVXc8f32{re_part,im_part});     
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 csinh(const AVXc8f32 x) {

					         const __m256 re_part =
						      _mm256_mul_ps(_mm256_sinh_ps(x.m_re),
						                    _mm256_cos_ps(x.m_im));
						 const __m256 im_part =
						      _mm256_mul_ps(_mm256_cosh_ps(x.m_re),
						                    _mm256_sin_ps(x.m_im));
						 return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 csinh(const std::complex<float> x) {
                                                 
                                                 const __m256 real = AVXC8F32_SETPS(x.real());
						 const __m256 imag = AVXC8F32_SETPS(x.imag());
						 const __m256 re_part =
						       _mm256_mul_ps(_mm256_sinh_ps(real),
						                     _mm256_cos_ps(imag));
						 const __m256 im_part =
						       _mm256_mul_ps(_mm256_cosh_ps(real),
						                     _mm256_sin_ps(imag));
						 return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16)
                                         static inline
					 AVXc8f32 csinh(const float re,
					                const float im) {

						 const __m256 real = AVXC8F32_SETPS(x.real());
						 const __m256 imag = AVXC8F32_SETPS(x.imag());
						 const __m256 re_part =
						       _mm256_mul_ps(_mm256_sinh_ps(real),
						                     _mm256_cos_ps(imag));
						 const __m256 im_part =
						       _mm256_mul_ps(_mm256_cosh_ps(real),
						                     _mm256_sin_ps(imag));
						 return (AVXc8f32{re_part,im_part});	
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 ccos(const AVXc8f32 x) {

					           const __m256 re_part =
                                                         _mm256_mul_ps(_mm256_cos_ps(x.m_re),
	                                                               _mm256_cosh_ps(x.m_im));
                                                   const __m256 im_part =
                                                         _mm256_mul_part(_mm256_sin_ps(x.m_re),
	                                                                 _mm256_sinh_ps(x.m_im));
                                                   return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 ccos(const std::complex<float> x) {

					           const __m256 real = AVXC8F32_SETPS(x.real());
						   const __m256 imag = AVXC8F32_SETPS(x.imag());
						   const __m256 re_part =
						         _mm256_mul_ps(_mm256_cos_ps(real),
							               _mm256_cosh_ps(imag));
						   const __m256 im_part =
						         _mm256_mul_ps(_mm256_sin_ps(real),
							               _mm256_sinh_ps(imag));
						   return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 ccos(const float re,
					               const float im) {


						   const __m256 real = AVXC8F32_SETPS(re);
						   const __m256 imag = AVXC8F32_SETPS(im);
						   const __m256 re_part =
						         _mm256_mul_ps(_mm256_cos_ps(real),
							               _mm256_cosh_ps(imag));
						   const __m256 im_part =
						         _mm256_mul_ps(_mm256_sin_ps(real),
							               _mm256_sinh_ps(imag));
						   return (AVXc8f32{re_part,im_part});    
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 ccosh(const AVXc8f32 x) {
					       
                                                   const __m256 re_part =
						         _mm256_mul_ps(_mm256_cosh_ps(x.m_re),
							               _mm256_cos_ps(x.m_im));
                                                   const __m256 im_part =
						         _mm256_mul_ps(_mm256_sinh_ps(x.m_re),
							               _mm256_sin_ps(x.m_im));
						   return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 ccosh(const std::complex<float> x) {

					            const __m256 real = AVXC8F32_SETPS(x.real());
						    const __m256 imag = AVXC8F32_SETPS(x.imag());
						    const __m256 re_part =
						          _mm256_mul_ps(_mm256_cosh_ps(real),
							                _mm256_cos_ps(imag));
						    const __m256 im_part =
						          _mm256_mul_ps(_mm256_sinh_ps(real),
							                _mm256_sin_ps(imag));
						    return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 ccosh(const float re,
					                const float im) {

					            const __m256 real = AVXC8F32_SETPS(re);
						    const __m256 imag = AVXC8F32_SETPS(im);
						    const __m256 re_part =
						          _mm256_mul_ps(_mm256_cosh_ps(real),
							                _mm256_cos_ps(imag));
						    const __m256 im_part =
						          _mm256_mul_ps(_mm256_sinh_ps(real),
							                _mm256_sin_ps(imag));
						    return (AVXc8f32{re_part,im_part});		

					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 cexp(const AVXc8f32 x) {

					            const __m256 re_part =
						          _mm256_mul_ps(_mm256_exp_ps(x.m_re),
							                _mm256_cos_ps(x.m_im));
						    const __m256 im_part =
						          _mm256_mul_ps(_mm256_exp_ps(x.m_re),
							                _mm256_sin_ps(x.m_im));
						    return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 cexp(const std::complex<float> x) {

                                                    const __m256 real = AVXC8F32_SETPS(x.real());
						    const __m256 imag = AVXC8F32_SETPS(x.imag());
						    const __m256 re_part =
						          _mm256_mul_ps(_mm256_exp_ps(real),
							                _mm256_cos_ps(imag));
						    const __m256 im_part =
						          _mm256_mul_ps(_mm256_exp_ps(real),
							                _mm256_sin_ps(imag));
						    return (AVXc8f32{re_part,im_part});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 cexp(const float re,
					               const float im) {

                                                    const __m256 real = AVXC8F32_SETPS(re);
						    const __m256 imag = AVXC8F32_SETPS(im);
						    const __m256 re_part =
						          _mm256_mul_ps(_mm256_exp_ps(real),
							                _mm256_cos_ps(imag));
						    const __m256 im_part =
						          _mm256_mul_ps(_mm256_exp_ps(real),
							                _mm256_sin_ps(imag));
						    return (AVXc8f32{re_part,im_part}); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 __m256 cabs(const AVXc8f32 x) {

					            const __m256 re_part =
						          _mm256_mul_ps(x.m_re,x.m_re);
						    const __m256 im_part =
						          _mm256_mul_ps(x.m_im,x.m_im);
						    return (_mm256_sqrt_ps(_mm256_add_ps(re_part,im_part)));
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 __m256 cabs(const std::complex<float> x) {

					            const __m256 real = AVXC8F32_SETPS(x.real());
						    const __m256 imag = AVXC8F32_SETPS(x.imag());
						    const __m256 re_part =
						          _mm256_mul_ps(real,real);
						    const __m256 im_part =
						          _mm256_mul_ps(imag,imag);
						    return (_mm256_sqrt_ps(_mm256_add_ps(re_part,im_part)));

					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 __m256 cabs(const float re,
					             const float im) {

                                                    const __m256 real = AVXC8F32_SETPS(re);
						    const __m256 imag = AVXC8F32_SETPS(im);
						    const __m256 re_part =
						          _mm256_mul_ps(real,real);
						    const __m256 im_part =
						          _mm256_mul_ps(imag,imag);
						    return (_mm256_sqrt_ps(_mm256_add_ps(re_part,im_part)));
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 cpow(const AVXc8f32 x,
					               const float n) {

                                                    const __m256 re_part =
                                                            _mm256_mul_ps(x.m_re,x.m_re);
                                                    const __m256 im_part =
                                                            _mm256_mul_ps(x.m_im,x.m_im);
                                                    const __m256 r =
                                                            _mm256_sqrt_ps(_mm256_add_ps(re_part,im_part));
                                                    const __m256 theta =
                                                            _mm256_atan_ps(_mm256_div_ps(x.m_im,x.m_re));
                                                    const __m256 vn = _mm256_set1_ps(n);
                                                    const __m256 pow_term = _mm256_pow_ps(r,vn);
                                                    const __m256 trig_arg = _mm256_mul_ps(vn,theta);
                                                    return (AVXc8f32{_mm256_mul_ps(pow_term,_mm256_cos_ps(trig_arg),
                                                                         _mm256_mul_ps(pow_term,_mm256_sin_ps(trig_arg))}));
			 		 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 clog(const AVXc8f32 x) {

					           const __m256 t1 = cabs(x);
						   const __m256 t2 = carg(x);
						   const __m256 re_part =
						         _mm256_log_ps(t1);
						   return (AVXc8f32{re_part,t2});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 clog(const std::complex<float> x) {

                                                   const __m256 t1 = cabs(c);
						   const __m256 t2 = carg(x.real(),x.imag());
						   const __m256 re_part =
						         _mm256_log_ps(t1);
						   return (AVXc8f32{re_part,t2});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 clog(const float re,
					                const float im) {

                                                   const __m256 t1 = cabs(c);
						   const __m256 t2 = carg(re,im);
						   const __m256 re_part =
						         _mm256_log_ps(t1);
						   return (AVXc8f32{re_part,t2});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 csqrt(const AVXc8f32 x) {
                                                  
                                                    const __m256 t = cabs(x);
                                                    const __m256 re_part =
                                                          _mm256_mul_ps(_mm256_set1_ps(0.5f),
	                                                                _mm256_add_ps(t,x.m_re));
                                                    const __m256 im_part =
                                                          _mm256_mul_ps(_mm256_set1_ps(0.5f),
	                                                                _mm256_sub_ps(t,x.m_re));
                                                    return (AVXc8f32{_mm256_sqrt_ps(re_part),
                                                                        _mm256_sqrt_ps(im_part)});
			             	   }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 csqrt(const std::complex<float> x) {

                                                    const __m256 t = cabs(c);
                                                    const __m256 real = AVXC8F32_SETPS(c.real())
                                                    const __m256 imag = AVXC8F32_SETPS(c.imag())
                                                    const __m256 re_part =
                                                            _mm256_mul_ps(_mm256_set1_ps(0.5f),
	                                                                  _mm256_add_ps(t,real));
                                                    const __m256 im_part =
                                                            _mm256_mul_ps(_mm256_set1_ps(0.5f),
	                                                                  _mm256_sub_ps(t,real));
                                                    return (AVXc8f32{_mm256_sqrt_ps(re_part),
                                                                          _mm256_sqrt_ps(im_part)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 csqrt(const float re,
					                 const float im) {

                                                    const __m256 t = cabs(c);
                                                    const __m256 real = AVXC8F32_SETPS(re)
                                                    const __m256 imag = AVXC8F32_SETPS(im)
                                                    const __m256 re_part =
                                                            _mm256_mul_ps(_mm256_set1_ps(0.5f),
	                                                                  _mm256_add_ps(t,real));
                                                    const __m256 im_part =
                                                            _mm256_mul_ps(_mm256_set1_ps(0.5f),
	                                                                  _mm256_sub_ps(t,real));
                                                    return (AVXc8f32{_mm256_sqrt_ps(re_part),
                                                                          _mm256_sqrt_ps(im_part)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 ctan(const AVXc8f32 x) {

                                                    return (csin(x)/ccos(x));
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 ctan(const std::complex<float> x) {

					            return (csin(x)/ccos(x));
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 ctan(const float re,
					                const float im) {
                                                    
                                                    return (csin(re,im)/ccos(re,im));
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 ctanh(const AVXc8f32 x) {

                                                    return (csinh(x)/ccosh(x));
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 ctanh(const std::complex<float> x) {

                                                    return (csinh(x)/ccosh(x));
					  }
_ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 ctanh(const float re,
					                 const float im) {

                                                    return (csinh(re,im)/ccos(re,im));
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 blend_move(const AVXc8f32 x,
					                      const AVXc8f32 y,
							      const __m256 maskx,
							      const __m256 masky) {

                                                    const __m256 re = _mm256_blendv_ps(x.re,y.re,maskx);
						    const __m256 im = _mm256_blendv_ps(x.im,y.im,masky);
						    return (AVXc8f32{re,im});
					   }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  __m256 abs_real(const __m256 x) {

					            static const __m256 mask = _mm256_set1_ps(0x7FFFFFFF);
						    return (_mm256_and_ps(x,mask));
					  }
__ATTR_HOT__ __ATTR_ALIGN__(32) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 cdiv_smith(const AVXc8f32 x,
					                      const AVXc8f32 y) {

                                                    __m256 ratio,denom,re_part,im_part;
						    __m256 mgte;
						    //
						    mgte = _mm256_setzero_ps();
						    mgte = _mm256_cmp_ps(abs_real(y.m_re),
						                         abs_real(y.m_im),
									 _CMP_GE_OQ);
						    ratio = _mm256_setzero_ps();
						    denom = _mm256_setzero_ps();
						    if(_mm256_testz_ps(mgte,mgte)) {
                                                           ratio = _mm256_div_ps(y.m_im,y.m_re);
	                                                   denom = _mm256_add_ps(y.m_re,
	                                                          _mm256_mul_ps(ratio,y.m_im));
	                                                   re_part = _mm256_div_ps(_mm256_add_ps(x.m_re,
	                                                            _mm256_mul_ps(x.m_im,ratio)),denom);
	                                                   im_part = _mm256_div_ps(_mm256_sub_ps(x.m_im,
	                                                            _mm256_mul_ps(x.m_re,ratio)),denom);
	                                                  return (AVXc8f32{re_part,im_part});
						    } else {
                                                            ratio   = _mm256_div_ps(y.m_re,y.m_im);
	                                                    denom   = _mm256_add_ps(y.m_im,_mm256_mul_ps(ratio,y.m_re));
	                                                    re_part = _mm256_div_ps(_mm256_add_ps(
	                                                              _mm256_mul_ps(x.m_re,ratio)),denom);
	                                                    im_part = _mm256_div_ps(_mm256_sub_ps(
	                                                              _mm256_mul_ps(x.m_im,ratio)),denom);
	                                                    return (AVXc8f32{re_part,im_part});
						    }
						    
					   }

__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+(const AVXc8f32 x,
					                     const AVXc8f32 y) {

					           return (AVXc8f32{_mm256_add_ps(x.m_re,y.m_re),
						                    _mm256_add_ps(x.m_im,y.m_im)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+(const AVXc8f32 x,
					                     const __m256 y) {

					           return(AVXc8f32{_mm256_add_ps(x.m_re,y),
						                   _mm256_add_ps(x.m_im,y)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+(const AVXc8f32 x,
					                     const float s) {

					           return (x+AVXc8f32{s});		     
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+(const __m256 x,
					                     const AVXc8f32 y) {

					           return (AVXc8f32{_mm256_add_ps(x,y.m_re),
						                    _mm256_add_ps(x,y.m_im)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+(const float s,
					                     const AVXc8f32 x) {

                                                   return (AVXc8f32{s}+x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+=(AVXc8f32 x,
					                      const AVXc8f32 y) {

                                                   x = x+y;
						   return (x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+=(AVXc8f32 x,
					                      const __m256 y) {

                                                   x = x+y;
						   return (x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+=(const __m256 x,
					                      AVXc8f32 y) {

						   y = y+x;
						   return (y);

					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+=(AVXc8f32 x,
					                      const float s) {

                                                   x = x+AVXc8f32{s};
						   return (x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator+=(const float s,
					                      AVXc8f32 x) {

                                                   x = AVXc8f32{s}+x;
						   return (x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-(const AVXc8f32 x,
					                     const AVXc8f32 y) {

                                                   return (AVXc8f32{_mm256_sub_ps(x.m_re,y.m_re),
						                    _mm256_sub_ps(x.m_im,y.m_im)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-(const AVXc8f32 x,
					                     const __m256 y) {

                                                   return (AVXc8f32{_mm256_sub_ps(x.m_re,y),
						                    _mm256_sub_ps(x.m_im,y)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-(const __m256 x,
					                     const AVXc8f32 y) {

                                                   return (AVXc8f32{_mm256_sub_ps(y.m_re,x),
						                    _mm256_sub_ps(y.m_im,x)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-(const AVXc8f32 x,
					                     const float s) {

					        return (x-AVXc8f32{s});   		     
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-(const float s,
					                     const AVXc8f32 x) {

                                                return (AVXc8f32{s}-x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-=(AVXc8f32 x,
					                      const AVXc8f32 y) {

                                                x = x-y;
						return (x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-=(AVXc8f32 x,
					                      const __m256 y) {
                                               
                                                x = x-y;
						return (x);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-=(const __m256 x,
					                      AVXc8f32 y) {

                                                y = y-x;
						return (y);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-=(const float s,
					                      AVXc8f32 x) {

                                                x = x-AVXc8f32{s};
						return (s);
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
					  AVXc8f32 operator-=(AVXc8f32 x,
					                      const float s) {
                                           
                                                x = AVXc8f32{s}-x;
						return (x);
					  }
 __ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                          static inline
                                          AVXc8f32 operator*(const AVXc8f32 x,
					                     const AVXc8f32 y) {

                                                 const __m256 ymm0 = _mm256_mul_ps(x.m_re,y.m_re);
                                                 const __m256 ymm1 = _mm256_mul_ps(x.m_im,y.m_im);
                                                 const __m256 ymm2 = _mm256_mul_ps(x.m_im,y.m_re);
                                                 const __m256 ymm3 = _mm256_mul_ps(x.m_re,y.m_im);
                                                 return (AVXc8f32{_mm256_sub_ps(ymm0,ymm1),
                                                                  _mm256_sub_ps(ymm2,ymm3)});
					  }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*(const AVXc8f32 x,
					                    const __m256 y) {
                
                                                 return (AVXc8f32{_mm256_mul_ps(x.m_re,y),
						                  _mm256_mul_ps(x.m_im,y)});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*(const __m256 x,
					                    const AVXc8f32 y) {

                                                 return (AVXc8f32{_mm256_mul_ps(x,y.m_re),
						                  _mm256_mul_ps(x,y.m_im)});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*(const AVXc8f32 x,
					                    const float s) {

						 return (x-AVXc8f32{s});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*(const float s,
					                    AVXv8f32 x) {

                                                 return (AVXc8f32{s}-x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*=(AVXc8f32 x,
					                     const AVXc8f32 y) {

						 x = x*y;
						 return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*=(AVXc8f32 x,
					                     const __m256 y) {

						 x = x*y;
						 return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*=(const __m256 x,
					                     AVXc8f32 y) {

                                                 y = y*x;
						 return (y);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*=(AVXc8f32 x,
					                     const float s) {

                                                 x = x*AVXc8f32{s};
						 return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator*=(const float s,
					                     AVXc8f32 x) {

                                                 x = AVXc8f32{s}*x;
						 return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/(const AVXc8f32 x,
					                    const AVXc8f32 y) {
#if defined USE_SAFE_COMPLEX_DIVISION && (USE_SAFE_COMPLEX_DIVISION) == 1
                                                 return (cdiv_smith(x,y));
#else
                                                 const __m256 ymm0 = _mm256_mul_ps(x.m_re,y.m_re);
                                                 const __m256 ymm1 = _mm256_mul_ps(x.m_im,y.m_im);
                                                 const __m256 ymm2 = _mm256_mul_ps(x.m_im,y.m_re);
                                                 const __m256 ymm3 = _mm256_mul_ps(x.m_re,y.m_im);
                                                 const __m256 den  = _mm256_add_ps(_mm256_mul_ps(y.m_re,y.m_re),
                                                                                   _mm256_mul_ps(y.m_im,y.m_im));
                                                 const __m256 re_part = _mm256_add_ps(ymm0,ymm1);
                                                 const __m256 im_part = _mm256_sub_ps(ymm2,ymm3);
                                                 return (AVXc8f32{_mm256_div_ps(re_part,den),
                                                                     _mm256_div_ps(im_part,den)});
#endif

					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/(const AVXc8f32 x,
					                    const __m256 y) {

                                                  return (AVXc8f32{_mm256_div_ps(x.m_re,y),
						                   _mm256_div_ps(x.m_im,y)});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/(const __m256 x,
					                    const AVXc8f32 y) {

                                                  return (AVXc8f32{_mm256_div_ps(x,y.m_re),
						                   _mm256_div_ps(x,y.m_im)});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/(const AVXc8f32 x,
					                    const float s) {

                                                  return (x/AVXc8f32{s});
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/(const float s,
					                    const AVXc8f32 x) {

                                                  return (AVXc8f32{s}/x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/=(AVXc8f32 x,
					                     const AVXc8f32 y) {

                                                  x = x/y;
						  return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/=(AVXc8f32 x,
					                     const __m256 y) {

                                                  x = x/y;
						  return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/=(const __m256 x,
					                     AVXc8f32 y) {

                                                  y = y/x;
						  return (y);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/=(const float s,
					                     AVXc8f32 x) {

                                                  x = AVXc8f32{s}/s;
						  return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator/=(AVXc8f32 x,
					                     const float s) {

                                                  x = x/AVXc8f32{s};
						  return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 AVXc8f32 operator~(AVXc8f32 x) {
                                                  x.m_re = _mm256_sub_ps(_mm256_setzero_ps(),x.m_re);
						  return (x);
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator==(const AVXc8f32 x,
					            const AVXc8f32 y) {

						   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,y.m_re,
						                        _CMP_EQ_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,y.m_im,
						                        _CMP_EQ_OQ);
                                                   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim)));
 					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator==(const AVXc8f32 x,
					            const std::complex<float> y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,
						                        _mm256_set1_ps(y.real()),
									_CMP_EQ_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,
						                        _mm256_set1_ps(y.imag()),
									_CMP_EQ_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator==(const std::complex<float> x,
					            const AVXc8f32 y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(_mm256_set1_ps(x.real()),
						                        y.m_re,
									_CMP_EQ_OQ);
						   eqim = _mm256_cmp_ps(_mm256_set1_ps(x.imag()),
						                        y.m_im,
									_CMP_EQ_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator!=(const AVXc8f32 x,
					            const AVXc8f32 y) {

						   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,y.m_re,
						                        _CMP_NEQ_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,y.m_im,
						                        _CMP_NEQ_OQ);
                                                   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim)));
 					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator!=(const AVXc8f32 x,
					            const std::complex<float> y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,
						                        _mm256_set1_ps(y.real()),
									_CMP_NEQ_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,
						                        _mm256_set1_ps(y.imag()),
									_CMP_NEQ_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator!=(const std::complex<float> x,
					            const AVXc8f32 y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(_mm256_set1_ps(x.real()),
						                        y.m_re,
									_CMP_NEQ_OQ);
						   eqim = _mm256_cmp_ps(_mm256_set1_ps(x.imag()),
						                        y.m_im,
									_CMP_NEQ_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator>(const AVXc8f32 x,
					            const AVXc8f32 y) {

						   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,y.m_re,
						                        _CMP_GT_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,y.m_im,
						                        _CMP_GT_OQ);
                                                   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim)));
 					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator>(const AVXc8f32 x,
					            const std::complex<float> y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,
						                        _mm256_set1_ps(y.real()),
									_CMP_GT_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,
						                        _mm256_set1_ps(y.imag()),
									_CMP_GT_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator>(const std::complex<float> x,
					            const AVXc8f32 y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(_mm256_set1_ps(x.real()),
						                        y.m_re,
									_CMP_GT_OQ);
						   eqim = _mm256_cmp_ps(_mm256_set1_ps(x.imag()),
						                        y.m_im,
									_CMP_GT_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator<(const AVXc8f32 x,
					            const AVXc8f32 y) {

						   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,y.m_re,
						                        _CMP_LT_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,y.m_im,
						                        _CMP_LT_OQ);
                                                   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim)));
 					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator<(const AVXc8f32 x,
					            const std::complex<float> y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,
						                        _mm256_set1_ps(y.real()),
									_CMP_LT_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,
						                        _mm256_set1_ps(y.imag()),
									_CMP_LT_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator<(const std::complex<float> x,
					            const AVXc8f32 y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(_mm256_set1_ps(x.real()),
						                        y.m_re,
									_CMP_LT_OQ);
						   eqim = _mm256_cmp_ps(_mm256_set1_ps(x.imag()),
						                        y.m_im,
									_CMP_LT_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator>=(const AVXc8f32 x,
					            const AVXc8f32 y) {

						   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,y.m_re,
						                        _CMP_GE_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,y.m_im,
						                        _CMP_GE_OQ);
                                                   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim)));
 					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator>=(const AVXc8f32 x,
					            const std::complex<float> y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,
						                        _mm256_set1_ps(y.real()),
									_CMP_GE_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,
						                        _mm256_set1_ps(y.imag()),
									_CMP_GE_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator>=(const std::complex<float> x,
					            const AVXc8f32 y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(_mm256_set1_ps(x.real()),
						                        y.m_re,
									_CMP_GE_OQ);
						   eqim = _mm256_cmp_ps(_mm256_set1_ps(x.imag()),
						                        y.m_im,
									_CMP_GE_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator<=(const AVXc8f32 x,
					            const AVXc8f32 y) {

						   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,y.m_re,
						                        _CMP_LE_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,y.m_im,
						                        _CMP_LE_OQ);
                                                   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim)));
 					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator<=(const AVXc8f32 x,
					            const std::complex<float> y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(x.m_re,
						                        _mm256_set1_ps(y.real()),
									_CMP_LE_OQ);
						   eqim = _mm256_cmp_ps(x.m_im,
						                        _mm256_set1_ps(y.imag()),
									_CMP_LE_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
__ATTR_HOT__ __ATTR_ALIGN__(16) __ATTR_VECTORCALL__
                                         static inline
					 std::pair<int32_t,int32_t>
					 operator<=(const std::complex<float> x,
					            const AVXc8f32 y) {

                                                   const __m256 eqre,eqim;
						   eqre = _mm256_setzero_ps();
						   eqim = _mm256_setzero_ps();
						   eqre = _mm256_cmp_ps(_mm256_set1_ps(x.real()),
						                        y.m_re,
									_CMP_LE_OQ);
						   eqim = _mm256_cmp_ps(_mm256_set1_ps(x.imag()),
						                        y.m_im,
									_CMP_LE_OQ);
						   return (std::make_pair(_mm256_testz_ps(eqre,eqre),
						                          _mm256_testz_ps(eqim,eqim))); 
					 }
    } // math   

}  // gms


#endif /*__GMS_AVXC8F32_H__*/
