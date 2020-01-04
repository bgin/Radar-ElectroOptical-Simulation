
#ifndef __GMS_AVXC4F64_H__
#define __GMS_AVXC4F64_H__


namespace file_info {

      const unsigned int gGMS_AVXC4F64_MAJOR = 1;
      const unsigned int gGMS_AVXC4F64_MINOR = 0;
      const unsigned int gGMS_AVXC4F64_MICRO = 0;
      const unsigned int gGMS_AVXC4F64_FULLVER =
            1000U*gGMS_AVXC4F64_MAJOR+100U*gGMS_AVXC4F64_MINOR+10U*gGMS_AVXC4F64;
      const char * const pgGMS_AVXC4F64_CREATE_DATE = "03-01-2020 11:26 +00200 (FRI 03 DEC 2020 GMT+2)";
      const char * const pgGMS_AVXC4F64_BUILD_DATE  = __DATE__ " " __TIME__;
      const char * const pgGMS_AVXC4F64_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
      const char * const pgGMS_AVXC4F64_SYNOPSYS    = "AVX complex number class decomposed into real and imaginary parts stored as a 4-tuple";

}

#include <cstdint>
#include <immintrin.h>
#include <complex>
#include "GMS_config.h"

namespace gms {

        namespace math {

#if !defined(AVXC4F64_SETPS)
    #define AVXC4F64_SETPS(x) _mm256_set_pd(1.0,1.0,1.0,(x));
#endif

               struct AVXc4f64 {

                      __m256d m_re;
		      __m256d m_im;
		      
__ATTR_HOT__
__ATTR_ALIGN__(16)
                 AVXc4f64() {

                      m_re = _mm256_setzero_pd();
		      m_im = _mm256_setzero_pd()
                  }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                  AVXc4f64(const double * __restrict __ATTR_ALIGN__(32) re,
			   const double * __restrict __ATTR_ALIGN__(32) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                       m_re = (const double*)__builtin_assume_aligned(re,32);
		       m_im = (const double*)__builtin_assume_aligned(im,32);
#elif defined __ICC || defined __INTEL_COMPILER
                       __assume_aligned(re,32);
		       __assume_aligned(im,32);
#endif
                       m_re = _mm256_load_pd(&m_re[0]);
		       m_im = _mm256_load_pd(&m_im[0]);
		  }
__ATTR_HOT__
__ATTR_ALIGN__(16)
                 AVXc4f64(const double re,
		          const double im) {

		       m_re = _mm256_set1_pd(re);
		       m_im = _mm256_set1_pd(im);
		 }
__ATTR_HOT__
__ATTR_ALIGN__(16)
                 AVXc4f64(const std::complex<double> c) {

                       m_re = _mm256_set1_pd(c.real());
		       m_im = _mm256_set1_pd(c.imag());
		 }
__ATTR_HOT__
__ATTR_ALIGN__(16)
                 AVXc4f64(const double re) {

		       m_re = _mm256_set1_pd(re);
		       m_im = _mm256_setzero_pd();

		 }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64(const double re0,
		          const double re1,
			  const double re2,
			  const double re3) {

                        m_re = _mm256_setr_pd(re0,re1,re2,re3);
			m_im = _mm256_setzero_pd();
		 }

__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64(const double re0,
		          const double re1,
			  const double re2,
			  const double re3,
			  const double im0,
			  const double im1,
			  const double im2,
			  const double im3) {

		        m_re = _mm256_setr_pd(re0,re1,re2,re3);
			m_im = _mm256_setr_pd(im0,im1,im2,im3);
		}

__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64(const __m256d re,
		          const __m256d im) {

                        m_re = re;
			m_im = im;
		}

__ATTR_HOT__
__ATTR_ALIGN__(16)
                AVXc4f64(const AVXc4f64 &x) {

                        m_re = x.m_re;
			m_im = x.m_im;
		}

		~AVXc4f64() {}

__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64 & load_a(const double * __restrict __ATTR_ALIGN__(32) re,
		                   const double * __restrict __ATTR_ALIGN__(32) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                                   re = (const double*)__builtin_assume_aligned(re,32);
				   im = (const double*)__builtin_assume_aligned(im,32);
#elif defined __ICC || defined __INTEL_COMPILER
                                   __assume_aligned(re,32);
				   __assume_aligned(im,32);
#endif
                                   m_re = _mm256_load_pd(&re[0]);
				   m_im = _mm256_load_pd(&im[0]);
				   return (*this);
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64 & load_u(const double * __restrict re,
		                   const double * __restrict im) {

                                   m_re = _mm256_loadu_pd(&re[0]);
				   m_im = _mm256_loadu_pd(&im[0]);
				   return (*this);
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 void store_a(double * __restrict __ATTR_ALIGN__(32) re,
		              double * __restrict __ATTR_ALIGN__(32) im) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                                 re = (double)__builtin_assume_aligned(re,32);
				 im = (double*)__builtin_assume_aligned(im,32);
#elif defined __ICC || defined __INTEL_COMPILER
                                 __assume_aligned(re,32);
                                 __assume_aligned(im,32);
#endif
                                 _mm256_store_pd(&re[0],m_re);
				 _mm256_store_pd(&im[0],m_im);
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                void store_u(double * __restrict re,
		             double * __restrict im) {

                                 _mm256_storeu_pd(&re[0],m_re);
				 _mm256_storeu_pd(&im[0],m_im);
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 void store_nt(double * __restrict re,
		               double * __restrict im) {

                                 _mm256_stream_pd(&re[0],m_re);
				 _mm256_stream_pd(&im[0],m_im);
				 _mm_sfence();
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
                 void extract_1c(double &s0,
		                 double &s1,
				 const int32_t posx,
				 const int32_t posy) {

                                 double re[4] __ATTR_ALIGN__(32) = {};
				 double im[4] __ATTR_ALIGN__(32) = {};
				 store_a(&re[0],&im[0]);
				 s0 = re[posx&3];
				 s1 = im[posy&3];
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
                AVXc4f64 & insert_1c(const double s0,
		                     const double s1,
				     const int32_t posx,
				     const int32_t posy) {

                                     double mem[8] __ATTR_ALIGN__(32) = {};
				     store_a(&mem[0],&mem[4]);
				     mem[posx&0x7] = s0;
				     mem[posy&0x7] = s1;
				     load_a(&mem[0],&mem[4]);
				     return (*this);
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                void concatenate_a(double * __restrict __ATTR_ALIGN__(32) out) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                                     out  = (double*)__builtin_assume_aligned(out,32);
#elif defined __ICC || defined __INTEL_COMPILER
                                     __assume_aligned(out,32);
#endif
                                     store_a(&out[0],&out[4]);
		 }
		 
__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 void concatenate_u(double * __restrict  out) {

                                    store_u(&out[0],&out[4]);

		 }
		 
__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc8f32 &  maskload(double const * __restrict re,
                                      double  const * __restrict im,
				      const __m256i mask1,
				      const __m256i mask2) {

				     m_re = _mm256_maskload_pd(&re[0],mask1);
				     m_im = _mm256_maskload_pd(&im[0],mask2);
				     return (*this);
		 }
		 
__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 void maskstore(double * __restrict re,
				double * __restrict im,
				const __m256i mask1,
				const __m256i mask2) {

				    _mm256_maskstore_pd(&re[0],mask1,m_re);
				    _mm256_maskstore_pd(&im[0],mask2,m_im);
		 }

__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                AVXc4f64 & maskload(const AVXc4f64 c8,
				    const __m256i mask1,
				    const __m256i mask2) {

				    m_re = _mm256_maskload_pd(&c8.m_re[0],mask1);
				    m_im = _mm256_maskload_pd(&c8.m_im[0],mask2);
				    return (*this);
		 }
		 
__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                void maskstore(AVXc4f64 & c8,
			       const __m256i mask1,
			       const __m256i mask2) {

				   _mm256_maskstore_pd((double*)&c8.m_re[0],mask1,m_re);
				   _mm256_maskstore_pd((double*)&c8.m_im[0],mask1,m_im);
		  }
		  
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64 & blendv(const AVXc4f64 x
				   const __m256d maskx,
				   const AVXc4f64 y,
				   const __m256d masky) {
							 
                                    m_re = _mm256_blendv_pd(x.m_re,y.m_re,maskx);
				    m_im = _mm256_blendv_pd(x.m_im,y.m_im,masky);
				    return (*this);
		  }
		  
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64 & blend(const AVXc4f64 x,
				  const int32_t imm8x,
				  const AVXc4f64 y,
				  const int32_t imm8y) {

				    m_re = _mm256_blendv_pd(x.m_re,y.m_re,imm8x);
				    m_im = _mm256_blendv_pd(x.m_im,y.m_im,imm8y);
				    return (*this);
		   }

__ATTR_COLD__
__ATTR_ALIGN__(16)
               __m128d re_lo() const {
					
                    return (_mm256_extractf128_pd(m_re,0));
		 }
		 
__ATTR_COLD__
__ATTR_ALIGN__(16)
              __m128d re_hi() {
                    return (_mm256_extractf128_pd(m_re,1));

		 }
		 
__ATTR_COLD__
__ATTR_ALIGN__(16)
              __m128d im_lo() const {
                    return (_mm256_extractf128_pd(m_im,0));

		}
		
__ATTR_COLD__
__ATTR_ALIGN__(16)
              __m128d im_hi() const {
                    return (_mm256_extractf128_pd(m_im,1));
                }
		
__ATTR_COLD__
__ATTR_ALIGN__(16) 
               AVXc4f64 & permute(const int32_t imm8x,
				  const int32_t imm8y) {

			m_re = _mm256_permute_pd(m_re,imm8x);
			m_im = _mm256_permute_pd(m_im,imm8y);
			return (*this);
	        }
		
__ATTR_COLD__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
               AVXc4f64 & permute(const AVXc4f64 x,
				  const AVXc4f64 y,
				  const int32_t imm8x,
				  const int32_t imm8y) {

			 m_re = _mm256_permute_pd(x.m_re,imm8x);
			 m_im = _mm256_permute_pd(x.m_im,imm8y);
			 return (*this);
		}

__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 AVXc4f64 & operator(const AVXc4f64 x) {
                           if(__builtin_expect(x == &this,0)) return (*this);
			   m_re = x.m_re;
			   m_im = x.m_im;
			   return (*this);
		}
              
	   }__ATTR_ALIGN__(64);


__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
		 static inline
		 AVXc4f64 conj(const AVXc4f64 x) {

			 auto tmp = ~x;
			 return (AVXc4f64{x.m_re,tmp.im});
		 }
		 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
                 AVXc4f64 polar(const __m256d ro,
				const __m256d theta) {

		         const __m256d re_part =
				      _mm256_mul_pd(rho,_mm256_cos_pd(theta));
			 const __m256d im_part =
				      _mm256_mul_pd(rho,_mm256_sin_pd(theta));
			 return (AVXc4f64{re_part,im_part});
		 }
		 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 __m256d carg(const AVXc4f64 x) {

			 return (_mm256_atan2_pd(x.m_re,x.m_im));
		 }
		 
__ATTR_HOT__
__ATTR_ALIGN__(16)
                 static inline
		 __256d carg(const double re,
			      const double im) {

		         const __m256d real = AVXC4F64_SETPS(re);
			 const __m256d imag = AVXC4F64_SETPS(im);
			 return (_mm256_atan2_pd(real,imag));
		 }
		 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csin(const AVXc4f64 x) {

		       const __m256d re_part =
				   _mm256_mul_pd(_mm256_sin_pd(x.m_re),
					         _mm256_cosh_pd(x.m_im));
		       const __m256d im_part =
				   _mm256_mul_pd(_mm256_cos_pd(x.m_re),
						 _mm256_sinh_pd(x.m_im));
		       return (AVXc4f64{re_part,im_part});
		  }
		  
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csin(const std::complex<double> x) {

			const __m256d real = AVXC4F64_SETPS(x.real());
			const __m256d imag = AVXC4F64_SETPS(x.imag());
			const __m256 re_part =
				     _mm256_mul_pd(_mm256_sin_pd(real),
						  _mm256_cosh_pd(imag));
		        const __m256d im_part =
				     _mm256_mul_pd(_mm256_cos_pd(real),
						  _mm256_sinh_pd(imag));
			return (AVXc4f64{re_part,im_part});
		  }
		  
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csin(const double re,
			       const double im) {

			 const __m256d real = AVXC4F64_SETPS(re);
			 const __m256d imag = AVXC4F64_SETPS(im);
			 const __m256d re_part =
				    _mm256_mul_pd(_mm256_sin_pd(real),
						  _mm256_cosh_pd(imag));
			 const __m256d im_part =
				    _mm256_mul_pd(_mm256_cos_pd(real),
						  _mm256_sinh_pd(imag));
			 return (AVXc4f64{re_part,im_part});     
		  }
		  
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csinh(const AVXc4f64 x) {

			  const __m256d re_part =
				 _mm256_mul_pd(_mm256_sinh_pd(x.m_re),
					       _mm256_cos_pd(x.m_im));
			  const __m256d im_part =
				 _mm256_mul_pd(_mm256_cosh_pd(x.m_re),
					       _mm256_sin_pd(x.m_im));
			  return (AVXc4f64{re_part,im_part});
		  }
		  
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csinh(const std::complex<double> x) {
                                                 
                          const __m256d real = AVXC4F64_SETPS(x.real());
			  const __m256d imag = AVXC4F64_SETPS(x.imag());
			  const __m256d re_part =
					_mm256_mul_pd(_mm256_sinh_pd(real),
						      _mm256_cos_pd(imag));
			  const __m256d im_part =
					_mm256_mul_pd(_mm256_cosh_pd(real),
						      _mm256_sin_pd(imag));
			  return (AVXc4f64{re_part,im_part});
		   }
		   
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
               static inline
	       AVXc4f64 csinh(const double re,
			      const double im) {

		          const __m256d real = AVXC4F64_SETPS(x.real());
			  const __m256d imag = AVXC4F64_SETPS(x.imag());
			  const __m256d re_part =
				       _mm256_mul_pd(_mm256_sinh_pd(real),
						     _mm256_cos_pd(imag));
			  const __m256d im_part =
				       _mm256_mul_pd(_mm256_cosh_pd(real),
						     _mm256_sin_pd(imag));
			  return (AVXc4f64{re_part,im_part});	
		    }

__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 ccos(const AVXc8f32 x) {

			    const __m256d re_part =
                                           _mm256_mul_pd(_mm256_cos_pd(x.m_re),
	                                                 _mm256_cosh_pd(x.m_im));
                             const __m256d im_part =
                                           _mm256_mul_pd(_mm256_sin_pd(x.m_re),
	                                                  _mm256_sinh_pd(x.m_im));
                             return (AVXccf64{re_part,im_part});
		    }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ccos(const std::complex<double> x) {

			      const __m256d real = AVXC4F64_SETPS(x.real());
			      const __m256d imag = AVXC4F64_SETPS(x.imag());
			      const __m256d re_part =
					    _mm256_mul_pd(_mm256_cos_pd(real),
							   _mm256_cosh_pd(imag));
			      const __m256d im_part =
					    _mm256_mul_pd(_mm256_sin_pd(real),
							   _mm256_sinh_pd(imag));
			      return (AVXc4f64{re_part,im_part});
		     }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 ccos(const double re,
			      const double im) {


			      const __m256d real = AVXC4F64_SETPS(re);
			      const __m256d imag = AVXC4F64_SETPS(im);
			      const __m256d re_part =
					    _mm256_mul_pd(_mm256_cos_pd(real),
							  _mm256_cosh_pd(imag));
			      const __m256d im_part =
					    _mm256_mul_pd(_mm256_sin_pd(real),
							  _mm256_sinh_pd(imag));
			      return (AVXc4f64{re_part,im_part});    
		      }

__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ccosh(const AVXc4f64 x) {
					       
                              const __m256d re_part =
				      _mm256_mul_pd(_mm256_cosh_pd(x.m_re),
						    _mm256_cos_pd(x.m_im));
                              const __m256d im_part =
				      _mm256_mul_pd(_mm256_sinh_pd(x.m_re),
						    _mm256_sin_pd(x.m_im));
			      return (AVXc4f64{re_part,im_part});
		      }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ccosh(const std::complex<double> x) {

			  const __m256d real = AVXC4F64_SETPS(x.real());
			  const __m256d imag = AVXC4F64_SETPS(x.imag());
			  const __m256d re_part =
					 _mm256_mul_pd(_mm256_cosh_pd(real),
						       _mm256_cos_pd(imag));
			  const __m256d im_part =
					 _mm256_mul_pd(_mm256_sinh_pd(real),
						       _mm256_sin_pd(imag));
			  return (AVXc4f64{re_part,im_part});
		      }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ccosh(const double re,
				const double im) {

			   const __m256d real = AVXC4F64_SETPS(re);
			   const __m256d imag = AVXC4F64_SETPS(im);
			   const __m256d re_part =
					 _mm256_mul_pd(_mm256_cosh_pd(real),
						       _mm256_cos_pd(imag));
			   const __m256d im_part =
					  _mm256_mul_pd(_mm256_sinh_pd(real),
						        _mm256_sin_pd(imag));
			   return (AVXc4f64{re_part,im_part});		

		       }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 cexp(const AVXc4f64 x) {

			   const __m256d re_part =
					 _mm256_mul_pd(_mm256_exp_pd(x.m_re),
						       _mm256_cos_pd(x.m_im));
			   const __m256d im_part =
					 _mm256_mul_pd(_mm256_exp_pd(x.m_re),
						       _mm256_sin_pd(x.m_im));
			   return (AVXc4f64{re_part,im_part});
		       }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 cexp(const std::complex<double> x) {

                           const __m256d real = AVXC4F64_SETPS(x.real());
			   const __m256d imag = AVXC4F64_SETPS(x.imag());
			   const __m256d re_part =
				         _mm256_mul_pd(_mm256_exp_pd(real),
						       _mm256_cos_pd(imag));
			   const __m256d im_part =
					 _mm256_mul_pd(_mm256_exp_pd(real),
						       _mm256_sin_pd(imag));
						    return (AVXc4f64{re_part,im_part});
			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 cexp(const double re,
			       const double im) {

                             const __m256d real = AVXC4F64_SETPS(re);
			     const __m256d imag = AVXC4F64_SETPS(im);
			     const __m256 re_part =
					      _mm256_mul_pd(_mm256_exp_pd(real),
							    _mm256_cos_pd(imag));
			     const __m256 im_part =
					      _mm256_mul_pd(_mm256_exp_pd(real),
							    _mm256_sin_pd(imag));
			     return (AVXc4f64{re_part,im_part}); 
			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 __m256d cabs(const AVXc4f64 x) {

			      const __m256d re_part =
					      _mm256_mul_pd(x.m_re,x.m_re);
			      const __m256d im_part =
					      _mm256_mul_pd(x.m_im,x.m_im);
			      return (_mm256_sqrt_pd(_mm256_add_pd(re_part,im_part)));
			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 __m256d cabs(const std::complex<double> x) {

			       const __m256d real = AVXC4F64_SETPS(x.real());
			       const __m256d imag = AVXC4F64_SETPS(x.imag());
			       const __m256d re_part =
						     _mm256_mul_pd(real,real);
			       const __m256d im_part =
					             _mm256_mul_pd(imag,imag);
			       return (_mm256_sqrt_pd(_mm256_add_pd(re_part,im_part)));

			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 __m256d cabs(const double re,
			      const double im) {

                                const __m256d real = AVXC4F64_SETPS(re);
				const __m256d imag = AVXC4F64_SETPS(im);
				const __m256d re_part =
						_mm256_mul_pd(real,real);
				const __m256d im_part =
						_mm256_mul_pd(imag,imag);
				return (_mm256_sqrt_pd(_mm256_add_pd(re_part,im_part)));
			 }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 cpow(const AVXc4f64 x,
			       const double n) {

                                const __m256d re_part =
                                                        _mm256_mul_pd(x.m_re,x.m_re);
                                const __m256d im_part =
                                                        _mm256_mul_pd(x.m_im,x.m_im);
                                const __m256d r =
                                                        _mm256_sqrt_pd(_mm256_add_pd(re_part,im_part));
                                const __m256d theta =
                                                         _mm256_atan_pd(_mm256_div_pd(x.m_im,x.m_re));
                                const __m256d vn = _mm256_set1_pd(n);
                                const __m256d pow_term = _mm256_pow_pd(r,vn);
                                const __m256d trig_arg = _mm256_mul_pd(vn,theta);
                                return (AVXc4f64{_mm256_mul_pd(pow_term,_mm256_cos_pd(trig_arg),
                                                _mm256_mul_pd(pow_term,_mm256_sin_pd(trig_arg))}));
			 }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 clog(const AVXc4f64 x) {

			    const __m256d t1 = cabs(x);
			    const __m256d t2 = carg(x);
			    const __m256d re_part =
						_mm256_log_pd(t1);
			    return (AVXc4f64{re_part,t2});
			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 clog(const std::complex<double> x) {

                           const __m256d t1 = cabs(c);
			   const __m256d t2 = carg(x.real(),x.imag());
			   const __m256d re_part =
						  _mm256_log_pd(t1);
			   return (AVXc4f64{re_part,t2});
			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 clog(const double re,
			       const double im) {

                            const __m256d t1 = cabs(c);
			    const __m256d t2 = carg(re,im);
			    const __m256d re_part =
						  _mm256_log_pd(t1);
			    return (AVXc4f64{re_part,t2});
			 }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csqrt(const AVXc4f64 x) {
                                                  
                             const __m256d t = cabs(x);
                             const __m256d re_part =
                                                _mm256_mul_pd(_mm256_set1_pd(0.5),
	                                                                _mm256_add_pd(t,x.m_re));
                              const __m256d im_part =
                                                _mm256_mul_pd(_mm256_set1_pd(0.5),
	                                                                _mm256_sub_pd(t,x.m_re));
                              return (AVXc4f64{_mm256_sqrt_pd(re_part),
                                               _mm256_sqrt_pd(im_part)});
			  }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csqrt(const std::complex<double> x) {

                                const __m256d t = cabs(c);
                                const __m256d real = AVXC4F64_SETPS(c.real())
                                const __m256d imag = AVXC4F64_SETPS(c.imag())
                                const __m256d re_part =
                                                  _mm256_mul_pd(_mm256_set1_pd(0.5),
	                                                                  _mm256_add_pd(t,real));
                                const __m256d im_part =
                                                  _mm256_mul_pd(_mm256_set1_pd(0.5),
	                                                                  _mm256_sub_pd(t,real));
                                return (AVXc4f64{_mm256_sqrt_pd(re_part),
                                                 _mm256_sqrt_pd(im_part)});
			  }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 csqrt(const double re,
				const double im) {

                               const __m256d t = cabs(c);
                               const __m256d real = AVXC4F64_SETPS(re)
                               const __m256d imag = AVXC4F64_SETPS(im)
                               const __m256d re_part =
                                                   _mm256_mul_pd(_mm256_set1_pd(0.5),
	                                                         _mm256_add_pd(t,real));
                               const __m256d im_part =
                                                   _mm256_mul_pd(_mm256_set1_pd(0.5),
	                                                          _mm256_sub_pd(t,real));
                               return (AVXc4f64{_mm256_sqrt_pd(re_part),
                                                _mm256_sqrt_pd(im_part)});
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ctan(const AVXc4f64 x) {

                             return (csin(x)/ccos(x));
			   }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ctan(const std::complex<double> x) {

			     return (csin(x)/ccos(x));
			   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ctan(const double re,
			       const double im) {
                                                    
                              return (csin(re,im)/ccos(re,im));
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ctanh(const AVXc4f64 x) {

                               return (csinh(x)/ccosh(x));
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ctanh(const std::complex<double> x) {

                                return (csinh(x)/ccosh(x));
			    }
 
_ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 ctanh(const double re,
				const double im) {

                                 return (csinh(re,im)/ccos(re,im));
			     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 blend_move(const AVXc4f64 x,
				     const AVXc4f64 y,
				     const __m256d maskx,
				     const __m256d masky) {

                        const __m256d re = _mm256_blendv_pd(x.re,y.re,maskx);
		        const __m256d im = _mm256_blendv_pd(x.im,y.im,masky);
			return (AVXc4f64{re,im});
			      }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 __m256d abs_real(const __m256d x) {

			 static const __m256d mask = _mm256_set1_pd(0x7FFFFFFF);
			 return (_mm256_and_pd(x,mask));
			      }
 
__ATTR_HOT__
__ATTR_ALIGN__(32)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 cdiv_smith(const AVXc4f64 x,
				     const AVXc4f64 y) {

                          __m256d ratio,denom,re_part,im_part;
			  __m256d mgte;
						    //
			  mgte = _mm256_setzero_pd();
			  mgte = _mm256_cmp_pd(abs_real(y.m_re),
					       abs_real(y.m_im),
					       _CMP_GE_OQ);
			  ratio = _mm256_setzero_pd();
			  denom = _mm256_setzero_pd();
			  if(_mm256_testz_pd(mgte,mgte)) {
                                ratio = _mm256_div_pd(y.m_im,y.m_re);
	                        denom = _mm256_add_pd(y.m_re,
	                                _mm256_mul_pd(ratio,y.m_im));
	                        re_part = _mm256_div_pd(_mm256_add_pd(x.m_re,
	                                                _mm256_mul_pd(x.m_im,ratio)),denom);
	                        im_part = _mm256_div_pd(_mm256_sub_pd(x.m_im,
	                                                _mm256_mul_pd(x.m_re,ratio)),denom);
	                        return (AVXc4f64{re_part,im_part});
						    } else {
                                ratio   = _mm256_div_pd(y.m_re,y.m_im);
	                        denom   = _mm256_add_pd(y.m_im,_mm256_mul_pd(ratio,y.m_re));
	                        re_part = _mm256_div_pd(_mm256_add_pd(
	                                                _mm256_mul_pd(x.m_re,ratio)),denom);
	                        im_part = _mm256_div_pd(_mm256_sub_pd(
	                                                _mm256_mul_pd(x.m_im,ratio)),denom);
	                        return (AVXc4f64{re_part,im_part});
			     }
						    
		    }

__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
	        AVXc4f64 operator+(const AVXc4f64 x,
				   const AVXc4f64 y) {

			  return (AVXc4f64{_mm256_add_pd(x.m_re,y.m_re),
					   _mm256_add_pd(x.m_im,y.m_im)});
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+(const AVXc4f64 x,
				    const __m256d y) {

			   return(AVXc4f64{_mm256_add_pd(x.m_re,y),
					   _mm256_add_pd(x.m_im,y)});
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+(const AVXc4f64 x,
				    const double s) {

			 return (x+AVXc4f64{s});		     
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+(const __m256d x,
				    const AVXc4f64 y) {

			 return (AVXc4f64{_mm256_add_pd(x,y.m_re),
					  _mm256_add_pd(x,y.m_im)});
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+(const double s,
				    const AVXc4f64 x) {

                          return (AVXc4f64{s}+x);
		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+=(AVXc4f64 x,
				     const AVXc4f64 y) {

                                  x = x+y;
				  return (x);
		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+=(AVXc4f64 x,
				     const __m256d y) {

                                    x = x+y;
				    return (x);
		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+=(const __m256d x,
				     AVXc4f64 y) {

				    y = y+x;
				    return (y);

		      }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+=(AVXc4f64 x,
				    const double s) {

                                     x = x+AVXc4f64{s};
				     return (x);
		       }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator+=(const double s,
				     AVXc4f64 x) {

                                     x = AVXc4f64{s}+x;
				     return (x);
			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator-(const AVXc4f64 x,
				    const AVXc4f64 y) {

                         return (AVXc4f64{_mm256_sub_pd(x.m_re,y.m_re),
					  _mm256_sub_pd(x.m_im,y.m_im)});
			 }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator-(const AVXc4f64 x,
				    const __m256d y) {

                           return (AVXc4f64{_mm256_sub_pd(x.m_re,y),
					    _mm256_sub_pd(x.m_im,y)});
			  }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator-(const __m256d x,
				    const AVXc4f64 y) {

                           return (AVXc4f64{_mm256_sub_pd(y.m_re,x),
					   _mm256_sub_pd(y.m_im,x)});
			   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator-(const AVXc4f64 x,
				    const double s) {

			     return (x-AVXc4f64{s});   		     
			   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator-(const double s,
				    const AVXc4f64 x) {

                            return (AVXc4f64{s}-x);
			   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator-=(AVXc4f64 x,
				   const AVXc4f64 y) {

                             x = x-y;
			     return (x);
		          }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator-=(AVXc4f64 x,
				    const __m256d y) {
                                               
                             x = x-y;
			     return (x);
			  }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator-=(const __m256d x,
				    AVXc4f64 y) {

                              y = y-x;
			      return (y);
			   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator-=(const double s,
				     AVXc4f64 x) {

                             x = x-AVXc4f64{s};
			     return (s);
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator-=(AVXc4f64 x,
				    const double s) {
                                           
                              x = AVXc4f64{s}-x;
			      return (x);
			     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
                 AVXc4f64 operator*(const AVXc4f64 x,
				    const AVXc4f64 y) {

                          const __m256d ymm0 = _mm256_mul_pd(x.m_re,y.m_re);
                          const __m256d ymm1 = _mm256_mul_pd(x.m_im,y.m_im);
                          const __m256d ymm2 = _mm256_mul_pd(x.m_im,y.m_re);
                          const __m256d ymm3 = _mm256_mul_pd(x.m_re,y.m_im);
                          return (AVXc4f64{_mm256_sub_pd(ymm0,ymm1),
                                           _mm256_sub_pd(ymm2,ymm3)});
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator*(const AVXc4f64 x,
				    const __m256d y) {
                
                            return (AVXc4f64{_mm256_mul_pd(x.m_re,y),
					     _mm256_mul_pd(x.m_im,y)});
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator*(const __m256d x,
				    const AVXc4f64 y) {

                              return (AVXc4f64{_mm256_mul_pd(x,y.m_re),
					       _mm256_mul_pd(x,y.m_im)});

                           }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator*(const AVXc4f64 x,
				    const double s) {

				 return (x-AVXc4f64{s});
			   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator*(const double s,
				    AVXv4f64 x) {

                             return (AVXc4f64{s}-x);
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
	        AVXc4f64 operator*=(AVXc4f64 x,
				    const AVXc4f64 y) {

			       x = x*y;
			       return (x);
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator*=(AVXc4f64 x,
				    const __m256d y) {

			       x = x*y;
			       return (x);
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator*=(const __m256d x,
				     AVXc4f64 y) {

                                y = y*x;
			        return (y);
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator*=(AVXc4f64 x,
				     const double s) {

                                 x = x*AVXc4f64{s};
				 return (x);
			     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator*=(const double s,
				     AVXc4f64 x) {

                                   x = AVXc4f64{s}*x;
				   return (x);
			   }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator/(const AVXc4f64 x,
				   const AVXc4f64 y) {
#if defined USE_SAFE_COMPLEX_DIVISION && (USE_SAFE_COMPLEX_DIVISION) == 1
               return (cdiv_smith(x,y));
#else
               const __m256d ymm0 = _mm256_mul_pd(x.m_re,y.m_re);
               const __m256d ymm1 = _mm256_mul_pd(x.m_im,y.m_im);
               const __m256d ymm2 = _mm256_mul_pd(x.m_im,y.m_re);
               const __m256d ymm3 = _mm256_mul_pd(x.m_re,y.m_im);
               const __m256d den  = _mm256_add_pd(_mm256_mul_pd(y.m_re,y.m_re),
                                                 _mm256_mul_pd(y.m_im,y.m_im));
               const __m256d re_part = _mm256_add_pd(ymm0,ymm1);
               const __m256d im_part = _mm256_sub_pd(ymm2,ymm3);
               return (AVXc4f64{_mm256_div_pd(re_part,den),
                                _mm256_div_pd(im_part,den)});
#endif

			}
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator/(const AVXc4f64 x,
				    const __m256d y) {

                          return (AVXc4f64{_mm256_div_pd(x.m_re,y),
					   _mm256_div_pd(x.m_im,y)});
			 }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator/(const __m256d x,
				    const AVXc4f64 y) {

                           return (AVXc4f64{_mm256_div_pd(x,y.m_re),
					   _mm256_div_pd(x,y.m_im)});
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator/(const AVXc4f64 x,
				    const double s) {

                           return (x/AVXc4f64{s});
			     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator/(const double s,
				     const AVXc4f64 x) {

                           return (AVXc4f64{s}/x);
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator/=(AVXc4f64 x,
				     const AVXc4f64 y) {

                              x = x/y;
			      return (x);
			  }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator/=(AVXc4f64 x,
				     const __m256d y) {

                              x = x/y;
			      return (x);
			  }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator/=(const __m256d x,
				    AVXc4f64 y) {

                              y = y/x;
			      return (y);
			  }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator/=(const double s,
				    AVXc4f64 x) {

                              x = AVXc4f64{s}/s;
			      return (x);
			   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 AVXc4f64 operator/=(AVXc4f64 x,
				     const double s) {

                               x = x/AVXc4f64{s};
			      return (x);
			    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		AVXc4f64 operator~(AVXc4f64 x) {
                       x.m_re = _mm256_sub_pd(_mm256_setzero_pd(),x.m_re);
		       return (x);
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator==(const AVXc4f64 x,
			    const AVXc4f64 y) {

			   const __m256d eqre,eqim;
			   eqre = _mm256_setzero_pd();
			   eqim = _mm256_setzero_pd();
			   eqre = _mm256_cmp_pd(x.m_re,y.m_re,
						      _CMP_EQ_OQ);
			   eqim = _mm256_cmp_pd(x.m_im,y.m_im,
						      _CMP_EQ_OQ);
                           return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						  _mm256_testz_pd(eqim,eqim)));
 		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
	         operator==(const AVXc4f64 x,
			    const std::complex<double> y) {

                             const __m256d eqre,eqim;
			     eqre = _mm256_setzero_pd();
			     eqim = _mm256_setzero_pd();
			     eqre = _mm256_cmp_pd(x.m_re,
						  _mm256_set1_pd(y.real()),
								_CMP_EQ_OQ);
			     eqim = _mm256_cmp_pd(x.m_im,
						  _mm256_set1_pd(y.imag()),
							        _CMP_EQ_OQ);
			     return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						    _mm256_testz_pd(eqim,eqim))); 
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator==(const std::complex<double> x,
			     const AVXc4f64 y) {

                             const __m256d eqre,eqim;
			     eqre = _mm256_setzero_pd();
			     eqim = _mm256_setzero_ps();
			     eqre = _mm256_cmp_pd(_mm256_set1_pd(x.real()),
						                 y.m_re,
								 _CMP_EQ_OQ);
			     eqim = _mm256_cmp_pd(_mm256_set1_pd(x.imag()),
						                 y.m_im,
								 _CMP_EQ_OQ);
			     return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						    _mm256_testz_pd(eqim,eqim))); 
		   }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator!=(const AVXc4f64 x,
			    const AVXc4f64 y) {

			  const __m256d eqre,eqim;
			  eqre = _mm256_setzero_pd();
			  eqim = _mm256_setzero_pd();
			  eqre = _mm256_cmp_pd(x.m_re,y.m_re,
						  _CMP_NEQ_OQ);
			  eqim = _mm256_cmp_pd(x.m_im,y.m_im,
						  _CMP_NEQ_OQ);
                          return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						 _mm256_testz_pd(eqim,eqim)));
 		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator!=(const AVXc4f64 x,
			    const std::complex<double> y) {

                            const __m256d eqre,eqim;
			    eqre = _mm256_setzero_pd();
			    eqim = _mm256_setzero_ps();
			    eqre = _mm256_cmp_pd(x.m_re,
						 _mm256_set1_ps(y.real()),
							     _CMP_NEQ_OQ);
			    eqim = _mm256_cmp_pd(x.m_im,
						 _mm256_set1_pd(y.imag()),
							      _CMP_NEQ_OQ);
			    return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						   _mm256_testz_pd(eqim,eqim))); 
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator!=(const std::complex<double> x,
			    const AVXc4f64 y) {

                            const __m256d eqre,eqim;
			    eqre = _mm256_setzero_pd();
			    eqim = _mm256_setzero_pd();
			    eqre = _mm256_cmp_pd(_mm256_set1_pd(x.real()),
						             y.m_re,
							     _CMP_NEQ_OQ);
			    eqim = _mm256_cmp_pd(_mm256_set1_pd(x.imag()),
						             y.m_im,
							    _CMP_NEQ_OQ);
			    return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						   _mm256_testz_pd(eqim,eqim))); 
		    }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator>(const AVXc4f64 x,
			   const AVXc4f64 y) {

			   const __m256d eqre,eqim;
			   eqre = _mm256_setzero_pd();
			   eqim = _mm256_setzero_pd();
			   eqre = _mm256_cmp_pd(x.m_re,y.m_re,
						       _CMP_GT_OQ);
			   eqim = _mm256_cmp_pd(x.m_im,y.m_im,
						       _CMP_GT_OQ);
                           return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						  _mm256_testz_pd(eqim,eqim)));
 		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator>(const AVXc4f64 x,
			  const std::complex<double> y) {

                          const __m256d eqre,eqim;
			  eqre = _mm256_setzero_pd();
			  eqim = _mm256_setzero_pd();
			  eqre = _mm256_cmp_pd(x.m_re,
					      _mm256_set1_pd(y.real()),
							    _CMP_GT_OQ);
			  eqim = _mm256_cmp_pd(x.m_im,
					      _mm256_set1_pd(y.imag()),
							    _CMP_GT_OQ);
						   return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						                          _mm256_testz_pd(eqim,eqim))); 
		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator>(const std::complex<double> x,
			   const AVXc4f64 y) {

                           const __m256d eqre,eqim;
			   eqre = _mm256_setzero_pd();
			   eqim = _mm256_setzero_pd();
			   eqre = _mm256_cmp_pd(_mm256_set1_pd(x.real()),
						                 y.m_re,
								 _CMP_GT_OQ);
			   eqim = _mm256_cmp_pd(_mm256_set1_pd(x.imag()),
						                 y.m_im,
								 _CMP_GT_OQ);
			   return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						 _mm256_testz_pd(eqim,eqim))); 
		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator<(const AVXc4f64 x,
			   const AVXc4f64 y) {

			       const __m256d eqre,eqim;
			       eqre = _mm256_setzero_pd();
			       eqim = _mm256_setzero_pd();
			       eqre = _mm256_cmp_pd(x.m_re,y.m_re,
						            _CMP_LT_OQ);
			       eqim = _mm256_cmp_pd(x.m_im,y.m_im,
						            _CMP_LT_OQ);
                               return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						      _mm256_testz_pd(eqim,eqim)));
 		     }
 
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator<(const AVXc4f64 x,
			  const std::complex<double> y) {

                             const __m256d eqre,eqim;
			     eqre = _mm256_setzero_pd();
			     eqim = _mm256_setzero_pd();
			     eqre = _mm256_cmp_pd(x.m_re,
						 _mm256_set1_pd(y.real()),
								_CMP_LT_OQ);
			     eqim = _mm256_cmp_pd(x.m_im,
						 _mm256_set1_pd(y.imag()),
								 CMP_LT_OQ);
			     return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						    _mm256_testz_pd(eqim,eqim))); 
		    }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator<(const std::complex<double> x,
			  const AVXc4f64 y) {

                               const __m256d eqre,eqim;
			       eqre = _mm256_setzero_pd();
			       eqim = _mm256_setzero_pd();
			       eqre = _mm256_cmp_pd(_mm256_set1_pd(x.real()),
						                  y.m_re,
								 _CMP_LT_OQ);
			       eqim = _mm256_cmp_pd(_mm256_set1_ps(x.imag()),
						                  y.m_im,
								 _CMP_LT_OQ);
			       return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						     _mm256_testz_pd(eqim,eqim))); 
		    }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator>=(const AVXc4f64 x,
			    const AVXc4f64 y) {

			    const __m256 eqre,eqim;
			    eqre = _mm256_setzero_pd();
			    eqim = _mm256_setzero_pd();
			    eqre = _mm256_cmp_pd(x.m_re,y.m_re,
						         _CMP_GE_OQ);
			    eqim = _mm256_cmp_pd(x.m_im,y.m_im,
						         _CMP_GE_OQ);
                            return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						  _mm256_testz_pd(eqim,eqim)));
 			  }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator>=(const AVXc4f64 x,
			   const std::complex<double> y) {

                           const __m256d eqre,eqim;
			   eqre = _mm256_setzero_pd();
			   eqim = _mm256_setzero_pd();
			   eqre = _mm256_cmp_pd(x.m_re,
						_mm256_set1_pd(y.real()),
							      _CMP_GE_OQ);
			   eqim = _mm256_cmp_pd(x.m_im,
						 _mm256_set1_pd(y.imag()),
							       _CMP_GE_OQ);
			   return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						 _mm256_testz_pd(eqim,eqim))); 
			   }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator>=(const std::complex<double> x,
			    const AVXc4f64 y) {

                             const __m256d eqre,eqim;
			     eqre = _mm256_setzero_pd();
			     eqim = _mm256_setzero_pd();
			     eqre = _mm256_cmp_pd(_mm256_set1_pd(x.real()),
						                 y.m_re,
								 _CMP_GE_OQ);
			     eqim = _mm256_cmp_pd(_mm256_set1_pd(x.imag()),
						                 y.m_im,
								 _CMP_GE_OQ);
			     return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						    _mm256_testz_pd(eqim,eqim))); 
					 }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                static inline
		std::pair<int32_t,int32_t>
	        operator<=(const AVXc4f64 x,
			   const AVXc4f64 y) {

			   const __m256d eqre,eqim;
			   eqre = _mm256_setzero_pd();
			   eqim = _mm256_setzero_pd();
			   eqre = _mm256_cmp_pd(x.m_re,y.m_re,
						      _CMP_LE_OQ);
			   eqim = _mm256_cmp_pd(x.m_im,y.m_im,
						      _CMP_LE_OQ);
                           return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						  _mm256_testz_pd(eqim,eqim)));
 		    }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator<=(const AVXc4f64 x,
			    const std::complex<double> y) {

                             const __m256d eqre,eqim;
			     eqre = _mm256_setzero_pd();
			     eqim = _mm256_setzero_pd();
			     eqre = _mm256_cmp_pd(x.m_re,
						  _mm256_set1_pd(y.real()),
							        _CMP_LE_OQ);
			     eqim = _mm256_cmp_pd(x.m_im,
						  _mm256_set1_pd(y.imag()),
								 _CMP_LE_OQ);
			    return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						  _mm256_testz_pd(eqim,eqim))); 
		       }
__ATTR_HOT__
__ATTR_ALIGN__(16)
__ATTR_VECTORCALL__
                 static inline
		 std::pair<int32_t,int32_t>
		 operator<=(const std::complex<double> x,
			    const AVXc4f64 y) {

                              const __m256d eqre,eqim;
			      eqre = _mm256_setzero_pd();
			      eqim = _mm256_setzero_pd();
			      eqre = _mm256_cmp_pd(_mm256_set1_pd(x.real()),
						                   y.m_re,
								   _CMP_LE_OQ);
			      eqim = _mm256_cmp_pd(_mm256_set1_pd(x.imag()),
						                  y.m_im,
								  _CMP_LE_OQ);
			      return (std::make_pair(_mm256_testz_pd(eqre,eqre),
						     _mm256_testz_pd(eqim,eqim))); 
		       }
	   

      } // math



} // gms


#endif /*__GMS_AVXC4F64_H__*/
