
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
	   

      } // math



} // gms


#endif /*__GMS_AVXC4F64_H__*/
