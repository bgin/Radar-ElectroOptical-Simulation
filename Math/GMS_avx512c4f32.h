

#ifndef __GMS_AVX512C4F32_H__
#define __GMS_AVX512C4F32_H__


namespace file_info {
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

    const unsigned int gGMS_AVX512C4F32_MAJOR = gms::common::gVersionInfo.m_VersionMajor;
    const unsigned int gGMS_AVX512C4F32_MINOR = gms::common::gVersionInfo.m_VersionMinor;
    const unsigned int gGMS_AVX512C4F32_MICRO = gms::common::gVersionInfo.m_VersionMicro;
    const unsigned int gGMS_AVX512C4F32_FULLVER =
          1000U*gGMS_AVX512C4F32_MAJOR+100U*gGMS_AVX512C4F32_MINOR+10U*gGMS_AVX512C4F32_MICRO;
    const char * const pgGMS_AVX512C4F32_CREATE_DATE = "19-11-2019 10:30 +00200 (TUE 19 NOV 2019 GMT+2)";
    const char * const pgGMS_AVX512C4F32_BUILD_DATE  = __DATE__ " " __TIME__;
    const char * const pgGMS_AVX512C4F32_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_AVX512C4F32_SYNOPSIS = "AVX512 complex number class decomposed into real and imaginary parts stored as 16-tuple.";
}

#include <cstdint>
#include <iostream>
#include <immintrin.h>
#include <complex> // for comparison operators
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#endif

namespace gms {
     namespace math {

#if !defined (AVX512C4F32_SETPS)
    #define AVX512C4F32_SET_PS(x)  _mm512_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,  \
                                                 1.0f,1.0f,1.0f,1.0f,1.0f,   \
						 1.0f,1.0f,1.0f,1.0f,1.0f,    \
						 (x));
#endif

#if !defined (MMASK16_RETTYPE)
    #define MMASK16_RETTYPE std::pair<__mmask16,__mmask16>
#endif

#if defined _WIN64
        __declspec(align(64)) struct AVX512c4Payload {
#elif defined __linux
                              struct AVX512c4Payload {
#endif
              float re0,re1,re2,re3,re4,re5,re6,re7,
	            re8,re9,re10,re11,re12,re13,re14,re15;
	      float im0,im1,im2,im3,im4,im5,im6,im7,
	            im8,im9,im10,im11,im12,im13,im14,im15;
         };
#if defined __linux
          __attribute__((aligned(64)))
#endif

#if defined _WIN64
        __declspec(align(64)) struct AVX512c4f32 {
#elif defined __linux
                              struct AVX512c4f32 {
#endif
               __m512 m_re;
	       __m512 m_im;

	       static const AVX512c4f32 C4F32_ZERO;
               // Sets components to 0.f
	       AVX512c4f32();

	       AVX512c4f32(const AVX512c4Payload);
               //  Arrays: Re and Im must be aligned on 64-bytes boundary, otherwise GP will be signalled.
	       AVX512c4f32(const float * __restrict,
	                   const float * __restrict);
	       // From single complex number
	       AVX512c4f32(const float,
	                   const float);
	       // From std::complex
	       AVX512c4f32(const std::complex<float>);
	       // From scalar
	       AVX512c4f32(const float);
	       // Real parts only
	       AVX512c4f32(const float,
	                   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float,
			   const float );
              // Main workhorse
	      AVX512c4f32(const __m512,
	                  const __m512);
	      AVX512c4f32(const AVX512c4f32 &);
	      ~AVX512c4f32() noexcept(true) = default;
              
	      // Load aligned
              AVX512c4f32 & load_a(const float * __restrict,
	                           const float * __restrict);
	      // Load unaligned
	      AVX512c4f32 & load_u(const float * __restrict,
	                           const float * __restrict);
	      // Store aligned
	      void store_a(float * __restrict,
	                   float * __restrict) const;
	      // Store unaligned
	      void store_u(float * __restrict,
	                   float * __restrict) const;
	      // Stream-store		   
	      void stream_nt(float * __restrict,
	                     float * __restrict) const;
	      //

	      float extract_1f32(const int32_t);			
				
	      std::pair<float,float>
	      extract_2f32(const int32_t,
	                   const int32_t);
			   
	      AVX512c4f32 & insert_1f32(const int32_t,
	                                const float);
					
	      AVX512c4f32 & insert_2f32(const int32_t,
	                                const int32_t,
					const float,
					const float);
	      // Length of 32 floats
	      void concatenate_a(float * __restrict) const;
	      // Length of 32 floats
	      void concatenate_u(float * __restrict) const;
	      
	      AVX512c4f32 & partial_loadu(const float const * __restrict,
	                                  const int32_t,
					  const float const * __restrict,
					  const int32_t);

	      AVX512c4f32 & partial_loada(const float const * __restrict,
	                                  const int32_t,
					  const float const * __restrict,
					  const int32_t);

	      void partial_storeu(float * __restrict,
	                          const int32_t,
				  float * __restrict,
				  const int32_t);

	      void partial_storea(float * __restrict,
	                          const int32_t,
				  float * __restrict,
				  const int32_t);

	      AVX512c4f32 & expand(const AVX512c4f32,
	                           const __mmask16);

	      AVX512c4f32 & expand_loadu(const AVX512c4f32,
	                                const __mmask16,
					const double * __restrict,
					const double * __restrict);

	      AVX512c4f32 & permute(const __mmask16,
	                            const int32_t);

	      __m256 re_low2() const
	      #if defined _WIN64
                  __declspec(const);
	      #elif defined __linux
	          __attribute__((const));
	      #endif
	      __m256 re_hi2()  const
	      #if defined _WIN64
                  __declspec(const);
	      #elif defined __linux
	          __attribute__((const));
	      #endif

	      __m256 im_low2() const
	        #if defined _WIN64
                  __declspec(const);
	      #elif defined __linux
	          __attribute__((const));
	      #endif

	      __m256 im_hi2()  const
	      #if defined _WIN64
                  __declspec(const);
	      #elif defined __linux
	          __attribute__((const));
	      #endif

	      AVX512c4f32 & operator=(const AVX512c4f32);
	      				
	};
#if defined __linux
           __attribute__((aligned(64)))
#endif

            static inline AVX512c4f32
	              conj(const AVX512c4f32);

	    static inline AVX512c4f32
	              polar(const __m512,
		            const __m512);

	    static inline __m512
	              carg(const AVX512c4f32);

	    static inline __m512
	              carg(const float,
		           const float);

	    static inline AVX512c4f32
	              csin(const AVX512c4f32);
		      
	    static inline AVX512c4f32
	              csin(const std::complex<float>);

	    static inline AVX512c4f32
	              csin(const float,
		           const float);

	    static inline AVX512c4f32
	              csinh(const AVX512c4f32);

	    static inline AVX512c4f32
	              csinh(const std::complex<float>);

	    static inline AVX512c4f32
	              csinh(const float,
		            const float);

	    static inline AVX512c4f32
	              ccos(const AVX512c4f32);

	    static inline AVX512c4f32
	              ccos(const std::complex<float>);

	    static inline AVX512c4f32
	              ccos(const float,
		           const float);

	    static inline AVX512c4f32
	              ccosh(const AVX512c4f32);

	    static inline AVX512c4f32
	              ccosh(const std::complex<float>);

	    static inline AVX512c4f32
	              ccosh(const float,
		            const float);

	    static inline AVX512c4f32
	              cexp(const AVX512c4f32);

	    static inline AVX512c4f32
	              cexp(const std::complex<float>);

	    static inline AVX512c4f32
	              cexp(const float,
		           const float);

	    static inline __m512
	              cabs(const AVX512c4f32);

	    static inline __m512
	              cabs(const std::complex<float>);

	    static inline __m512
	              cabs(const float,
		           const float);

	    static inline AVX512c4f32
	              cpow(const AVX512c4f32,
		           const float);

	    static inline AVX512c4f32
	              clog(const AVX512c4f32);

	    static inline AVX512c4f32
	              clog(const std::complex<float>);

	    static inline AVX512c4f32
	              clog(const float,
		           const float);

	    static inline AVX512c4f32
	              csqrt(const AVX512c4f32);

	    static inline AVX512c4f32
	              csqrt(const std::complex<float>);

	    static inline AVX512c4f32
	              csqrt(const float,
		            const float);

	    static inline AVX512c4f32
	              ctan(const AVX512c4f32);

	    static inline AVX512c4f32
	              ctan(const std::complex<float>);

	    static inline AVX512c4f32
	              ctan(const float,
		           const float);

	    static inline AVX512c4f32
	              ctanh(const AVX512c4f32);

	    static inline AVX512c4f32
	              ctanh(const std::complex<float>);

	    static inline AVX512c4f32
	              ctanh(const float,
		            const float);

	    static inline AVX512c4f32
	              select(const AVX512c4f32,
		             const AVX512c4f32,
			     const __mmask16);

	    static inline AVX512c4f32
	              cdiv_smith(const AVX512c4f32,
		                 const AVX512c4f32);

	    static inline AVX512c4f32
	    operator+(const AVX512c4f32,
	              const AVX512c4f32);

	    static inline AVX512c4f32
	    operator+(const AVX512c4f32,
	              const __m512);

	    static inline AVX512c4f32
	    operator+(const __m512,
	              const AVX512c4f32);

	    static inline AVX512c4f32
	    operator+(const AVX512c4f32,
	              const float);

	    static inline AVX512c4f32
	    operator+(const float,
	              const AVX512c4f32);

	    static inline AVX512c4f32
            operator+=(AVX512c4f32,
	               const AVX512c4f32);

	    static inline AVX512c4f32
	    operator+=(AVX512c4f32,
	               const __m512);

	    static inline AVX512c4f32
	    operator+=(const __m512,
	               AVX512c4f32);

	    static inline AVX512c4f32
	    operator+=(AVX512c4f32,
	               const float);

	    static inline AVX512c4f32
	    operator+=(const float,
	               AVX512c4f32);

	    static inline AVX512c4f32
	    operator-(const AVX512c4f32,
	              const AVX512c4f32);

	    static inline AVX512c4f32
	    operator-(const AVX512c4f32,
	              const __m512);

	    static inline AVX512c4f32
	    operator-(const __m512,
	              const AVX512c4f32);

	    static inline AVX512c4f32
	    operator-(const AVX512c4f32,
	              const float);

	    static inline AVX512c4f32
	    operator-(const float,
	              const AVX512c4f32);

	    static inline AVX512c4f32
	    operator-(AVX512c4f32);

	    static inline AVX512c4f32
	    operator-=(AVX512c4f32,
	               const AVX512c4f32);

	    static inline AVX512c4f32
	    operator-=(AVX512c4f32,
	               const __m512);

	    static inline AVX512c4f32
	    operator-=(const __m512,
	               AVX512c4f32);

	    static inline AVX512c4f32
	    operator-=(AVX512c4f32,
	               const float);

	    static inline AVX512c4f32
	    operator-=(const float,
	               AVX512c4f32);

	    static inline AVX512c4f32
	    operator*(const AVX512c4f32,
	              const AVX512c4f32);

	    static inline AVX512c4f32
	    operator*(const AVX512c4f32,
	              const __m512);

	    static inline AVX512c4f32
	    operator*(const __m512,
	              const AVX512c4f32);

            static inline AVX512c4f32
            operator*(const AVX512c4f32,
	              const float);

	    static inline AVX512c4f32
	    operator*(const float,
	              const AVX512c4f32);

	    static inline AVX512c4f32
	    operator*=(AVX512c4f32,
		       const AV512c4f32);

	    static inline AVX512c4f32
	    operator*=(AVX512c4f32,
		       const __m512);

	    static inline AVX512c4f32
	    operator*=(const __m512,
		       AVX512c4f32);

	    static inline AVX512c4f32
	    operator*=(AVX512c4f32,
		       const float);

	    static inline AVX512c4f32
	    operator*=(const float,
		       AVX512c4f32);

	    static inline AVX512c4f32
	    operator/(const AVX512c4f32,
		      const AVX512c4f32);

	    static inline AVX512c4f32
	    operator/(const AVX512c4f32,
		      const __m512);

	    static inline AVX512c4f32
	    operator/(const __m512,
		      const AVX512c4f32);

	    static inline AVX512c4f32
	    operator/(const AVX512c4f32,
		      const float);

	    static inline AVX512c4f32
	    operator/(const float,
		      const AVX512c4f32);

	    static inline AVX512c4f32
	    operator/=(AVX512c4f32,
		       const AVX512c4f32);

	    static inline AVX512c4f32
	    operator/=(AVX512c4f32,
		       const __m512);

	    static inline AVX512c4f32
	    operator/=(const __m512d,
	               AVX512c4f32);

	    static inline AVX512c4f32
	    operator/=(AVX512c4f32,
		       const float);

	    static inline AVX512c4f32
	    operator/=(const float,
		       AVX512c4f32);

	    static inline AVX512c4f32
	    operator~(AVX512c4f32);

		
	    static inline
	    MMASK16_RETTYPE
	    operator==(const AVX512c4f32,
		       const AVX512c4f32);

		

	    static inline
	    MMASK16_RETTYPE
	    operator==(const AVX512c4f32,
		       const std::complex<float>);

	    static inline
	    MMASK16_RETTYPE
	    operator==(const std::complex<float>,
		       const AVX512c4f32)
		           

	    static inline 
	    MMASK16_RETTYPE
	    operator!=(const AVX512c4f32,
		           const AVX512c4f32);

	    static inline
	    MMASK16_RETTYPE
	    operator!=(const AVX512c4f32,
		       const std::complex<float>);

	    static inline
	    MMASK16_RETTYPE
	    operator!=(const std::complex<float>,
		       const AVX512c4f32);

	    static inline 
	    MMASK16_RETTYPE
	    operator>(const AVX512c4f32,
		      const AVX512c14f32);

	    static inline
            MMASK16_RETTYPE
	    operator>(const AVX512c4f32,
		      const std::complex<float>);

	    static inline	  
	    MMASK16_RETTYPE
	    operator>(const std::complex<float>,
		      const AVX512c4f32);
			  
	    static inline 
	    MMASK16_RETTYPE
	    operator<(const AVX512c4f32,
		      const AVX512c4f32);

	    static inline	  
            MMASK16_RETTYPE
	    operator<(const AVX512c4f32,
		      std::complex<float>);

            static inline
	    MMASK16_RETTYPE
	    operator<(const std::complex<float>,
		      const AVX512c4f32);
			  
            static inline 
	    MMASK16_RETTYPE
	    operator>=(const AVX512c4f32,
		           const AVX512c4f32);

	    static inline
	    MMASK16_RETTYPE
	    operator>=(const AVX512c4f32,
		       const std::complex<float>);

	    static inline
	    MMASK16_RETTYPE
	    operator>=(const std::complex<float>,
		       const AVX512c4f32);

	    static inline 
	    MMASK16_RETTYPE
	    operator<=(const AVX512c4f32,
		       const AVX512c4f32);

	    static inline
	    MMASK16_RETTYPE
	    operator<=(const AVX512c4f32,
	    const std::complex<float>);

	    static inline
	    MMASK16_RETTYPE
	    operator<=(const std::complex<float>,
		       const AVX512c4f32);

			   
#include "GMS_avx512c4f32.inl"


   } // math
}  // gms








#endif /*__GMS_AVX512C4F32_H__*/
