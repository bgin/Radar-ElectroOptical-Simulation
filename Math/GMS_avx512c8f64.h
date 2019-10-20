
#ifndef __GMS_AVX512C8F64_H__
#define __GMS_AVX512C8F64_H__

namespace file_info {
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

	const unsigned int gGMS_AVX512C8F64_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_AVX512C8F64_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gLAM_AVX512C8F64_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_AVX512C8F64_FULLVER = 
		1000U*gGMS_AVX512C8F64_MAJOR+100U*gGMS_AVX3C8F64_MINOR+10U*gGMS_AVX3C8F64_MICRO;

	const char * const pgGMS_AVX512C8F64_CREATE_DATE = "29-09-2018 11:13 +00200 (SAT 29 SEP 2018 GMT+2)";

	const char * const pgGMS_AVX512C8F64_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_AVX512C8F64_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVX512C8F64_SYNOPSIS = "AVX512 complex number class decomposed into real and imaginary parts stored as 8-tuple.";
}

#include <cstdint>
#include <iostream>
#include <immintrin.h>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#endif

namespace gms {
	namespace math {


#if defined _WIN64
		__declspec(align(64)) struct AVX512c8Payload {
#elif defined __linux
		__attribute__((align(64))) struct AVX512c8Payload {
#endif
			double re0,re1,re2,re3,re4,re5,re6,re7;
			double im0,im1,im2,im3,im4,im5,im6,im7;
		};

		//static AVX3C8f64 CZERO{ _mm512_setzero_pd(),
		//						_mm512_setzero_pd() };
#if defined _WIN64
		__declspec(align(64)) struct AVX512c8f64 {
#elif defined __linux
		__attribute__((align(64))) struct AVX512c8f64 {
#endif
				__m512d m_re;
				
				__m512d m_im;

				static const AVX512c8f64 CZERO;
				
				// Default Ctor (sets components to 0.0)
				AVX512c8f64();
					
				

				AVX512c8f64(const AVX512c8Payload x);
					
				
				// Arrays: Re and Im must be aligned on 64-bytes boundary, otherwise GP will be signalled.
				AVX512c8f64(const double * __restrict ,
					    const double * __restrict);
					
				
				// From single complex number
				AVX512c8f64(const double ,
					    const double ); 
					
				
				// From scalar
				explicit AVX512c8f64(const double ); 
					
				

				// Real parts only
				AVX512c8f64(const double ,
				            const double ,
				            const double ,
					    const double,
				            const double ,
					    const double ,
					    const double ,
					    const double );
					
				

				AVX512c8f64(const __m512d,
					    const __m512d);
				

				AVX512c8f64(const AVX512c8f64 &);
					
				

				~AVX512c8f64() noexcept(true) = default;

				
				// Load-store functions
				// Load aligned
				AVX512c8f64 & load_a(const double * __restrict ,
						     const double * __restrict );
					
				
				// Load unaligned
				AVX512C8f64 & load_u(const double * __restrict ,
						     const double * __restrict );
					
				

				// Store aligned
				void store_a(double * __restrict ,
					     double * __restrict ) const;
					
				

				// Store unaligned
				void store_u(double * __restrict ,
				             double * __restrict ) const; 
					
				

				void stream(double * __restrict ,
					    double * __restrict ) const;
					

				double extract_f64(const int32_t) const;
					

				std::pair<double, double> 
				extract_2f64(const int32_t ,
					     const int32_t );
					
				
				AVX512c8f64 & insert_1f64(const int32_t ,
						          const double );
					

				AVX512c8f64 & insert_2f64(const int32_t ,
						          const int32_t ,
						          const double ,
						          const double );
					

				// Length of 16 doubles
				void concatenate_a(double * __restrict ) const;
					

				// Length of 16 doubles
				void concatenate_u(double * __restrict ) const;

				AVX512c8f64 & partial_loadu(const double const * __restrict,
							    const int32_t,
							    const double const * __restrict,
							    const int32_t );

				AVX512c8f64 & partial_loada(const double const * __restrict,
							    const int32_t,
							    const double const * __restrict,
							    const int32_t);

				void partial_storeu(double * __restrict,
					            const int32_t,
						    double * __restrict,
					            const int32_t);

				void partial_storea(double * __restrict,
						    const int32_t,
						    double * __restrict,
					            const int32_t);
					

				AVX512c8f64 & expand(const AVX512c8f64,
						     const __mmask8); 

				AVX512c8f64 & expand_load(const AVX3C8f64,
						          const __mmask8 ,
						          const double * __restrict ,
						          const double * __restrict );



				AVX512c8f64 & permute(const __mmask8 ,
						      const int32_t);
					

				__m256d re_low2() const; 

				

				__m256d re_hi2() const; 

				__m256d im_low2() const;

				__m256d im_hi2() const;

				

				AVX512c8f64 & operator=(_In_ const AVX3C8f64);

				
		};

		static inline AVX512c8f64
		        conj(const AVX512c8f64);

		static inline AVX512c8f64
		        polar(const __m512d,
			      const __m512d);

		static inline AVX512c8f64
		        carg(const AVX512c8f64);

		static inline AVX512c8f64
		        carg(const double,
			     const double);

		static inline AVX512c8f64 
			csin(const AVX5128f64);

		static inline AVX512c8f64
		        csin(const double,
			     const double);

		static inline AVX512c8f64
		        csinh(const AVX512c8f64);

		static inline AVX512c8f64
		        csinh(const double,
			      const double);

		static inline AVX512c8f64
			ccos(const AVX512c8f64);

		static inline AVX512c8f64
		        ccos(const double,
			     const double);

		static inline AVX512c8f64
		        ccosh(const AVX512c8f64);

		static inline AVX512c8f64
		        ccosh(const double,
			      const double);

		static inline AVX512c8f64
			cexp(const AVX512c8f64);

		static inline AVX512c8f64
		        cexp(const double,
			     const double);

		static inline AVX512c8f64
			cabs(const AVX512c8f64);

		static inline AVX512c8f64
		        cabs(const double,
			     const double);

		static inline AVX512c8f64
		        cpowi(const AVX512c8f64,
			      const double );

		static inline AVX512c8f64
		        clog(const AVX512c8f64);

		static inline AVX512c8f64
		        clog(const double,
			     const double);

		static inline AVX512c8f64
		        csqrt(const AVX512c8f64);

		static inline AVX512c8f64
		        csqrt(const double,
			      const double);

		static inline AVX512c8f64
		        ctan(const AVX512c8f64);

		static inline AVX512c8f64
		        ctan(const double,
			     const double);

		static inline AVX512c8f64
		        ctanh(const AVX512c8f64);

		static inline AVX512c8f64
		        ctanh(const double,
			      const double);
		
		static inline AVX512c8f64 
		select(const AVX512c8f64 ,
		       const AVX512c8f64,
		       const __mmask8);

		static inline AVX512c8f64 
		operator+(const AVX512c8f64,
			  const AVX512c8f64);

		static inline AVX512c8f64
		operator+(const AVX512c8f64,
		          const double );

		static inline AVX512c8f64
		operator+(const double ,
		          const AVX512c8f64); 

		static inline AVX512c8f64
		operator+=(AVX3C8f64,
		           const AVX512c8f64); 

		static inline AVX512c8f64
		operator+=(AVX512c8f64 &,
		          const double);

		static inline AVX512c8f64
		operator+=(const double,
			   AVX512c8f64);

		static inline AVX512c8f64
		operator-(const AVX512c8f64,
			  const AVX512c8f64); 

		static inline AVX512c8f64
		operator-(const AVX512c8f64,
			  const double );

		static inline AVX512c8f64
		operator-(const double ,
			  const AVX512c8f64);

		static inline AVX512c8f64
		operator-(AVX512c8f64);

		static inline AVX512c8f64
		operator-=(AVX512c8f64,
			   const AVX512c8f64);

		static inline AVX512c8f64
		operator-=(AVX512c8f64,
			   const double );

		static inline AVX512c8f64
		operator-=(const double,
		           AVX512c8f64);

		static inline AVX512c8f64
		operator*(const AVX512C8f64,
			  const AVX512c8f64);

		static inline AVX512c8f64
		operator*(const AVX512c8f64,
		          const double);

		static inline AVX512c8f64
		operator*(const double,
		          const AVX512c8f64);

		static inline AVX512c8f64
		operator*=(AVX512c8f64,
			   const AV512c8f64);

	        static inline AVX512c8f64
		operator*=(AVX512c8f64,
			   const double);

		static inline AVX512c8f64
		operator*=(const double,
		          AVX512c8f64);

		static inline AVX512c8f64
		operator/(const AVX512c8f64,
		          const AVX512c8f64);

		static inline AVX512c8f64
		operator/(const AVX512c8f64,
			  const double);

		static inline AVX512c8f64
		operator/(const double,
			  const AVX512c8f64);

		static inline AVX512c8f64
		operator/=(AVX512c8f64,
			   const AVX512c8f64);

	        static inline AVX512c8f64
		operator/=(AVX512c8f64,
			   const double);

		static inline AVX512c8f64
		operator/=(const double,
		           AVX512c8f64);

		static inline AVX512c8f64
		operator~(AVX512c8f64);

		
		static inline
			std::pair<__mmask8, __mmask8>
		operator==(const AVX512c8f64,
			   const AVX512c8f64);

		static inline 
			std::pair<__mmask8, __mmask8>
		operator!=(const AVX512c8f64,
		           const AVX512c8f64);

		static inline 
			std::pair<__mmask8, __mmask8>
		operator>(const AVX512c8f64,
			  const AVX512c8f64);

		static inline 
			std::pair<__mmask8, __mmask8>
		operator<(const AVX512c8f64,
			  const AVX512c8f64);

		static inline 
			std::pair<__mmask8, __mmask8>
		operator>=(const AVX512c8f64,
		           const AVX512c8f64);

		static inline 
			std::pair<__mmask8, __mmask8>
	        operator<=(const AVX512c8f64,
			   const AVX512c8f64);

#include "GMS_avx512c8f64.inl"	

	}
}



#endif /*__GMS_AVX512C8F64_H__*/
