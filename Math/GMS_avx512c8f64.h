
#ifndef __GMS_AVX512C8F64_H__
#define __GMS_AVX512C8F64_H__

namespace file_info {


        const unsigned int gGMS_AVX512C8F64_MAJOR = 1;

	const unsigned int gGMS_AVX512C8F64_MINOR = 1;

	const unsigned int gLAM_AVX512C8F64_MICRO = 0;

	const unsigned int gGMS_AVX512C8F64_FULLVER = 
		1000U*gGMS_AVX512C8F64_MAJOR+100U*gGMS_AVX3C8F64_MINOR+10U*gGMS_AVX3C8F64_MICRO;

	const char * const pgGMS_AVX512C8F64_CREATE_DATE = "29-09-2018 11:13 +00200 (SAT 29 SEP 2018 GMT+2)";

	const char * const pgGMS_AVX512C8F64_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const pgGMS_AVX512C8F64_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVX512C8F64_SYNOPSIS = "AVX512 complex number class decomposed into real and imaginary parts stored as 8-tuple.";
}

#include <cstdint>
#include <iostream>
#include <immintrin.h>
#include <complex> // for comparison operators


#include "GMS_config.h"


namespace gms {
	namespace math {



	        struct 	  __attribute__((aligned(64))) ZMM16c8Payload {

			double re0,re1,re2,re3,re4,re5,re6,re7;
			double im0,im1,im2,im3,im4,im5,im6,im7;
		};

	


		//static AVX3C8f64 CZERO{ _mm512_setzero_pd(),
		//						_mm512_setzero_pd() };


	        struct 	  __attribute__((aligned(64))) ZMM16c8 {

				__m512d m_re;
				
				__m512d m_im;

				static const ZMM16c8 CZERO;
				
				// Default Ctor (sets components to 0.0)
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
		                 __ATTR_VECTORCALL__
				ZMM16c8();
					
				 __ATTR_ALWAYS_INLINE__
                                 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
		                 __ATTR_VECTORCALL__
				 ZMM16c8(const ZMM16c8Payload);
					
				
				// Arrays: Re and Im must be aligned on 64-bytes boundary, otherwise GP will be signalled.
				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
		                 __ATTR_VECTORCALL__
				 ZMM16c8(const double * __restrict ,
					 const double * __restrict);
					
				
				// From single complex number
				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM16c8(const double ,
					 const double ); 
					
				
				// From scalar
				  __ATTR_ALWAYS_INLINE__
				  __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				explicit  ZMM16c8(const double ); 
					
				

				// Real parts only
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__  
				 ZMM16c8(   const double ,
				            const double ,
				            const double ,
					    const double,
				            const double ,
					    const double ,
					    const double ,
					    const double );
					
				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM16c8(const __m512d,
					 const __m512d);

		               
				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM16c8(const  ZMM16c8 &);
					
				
                                 
				~ZMM16c8() noexcept(true) = default;

				
				// Load-store functions
				// Load aligned
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM16c8 & load_a(const double * __restrict ,
						     const double * __restrict );
					
				
				// Load unaligned
				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				 ZMM16c8 & load_u(const double * __restrict ,
						     const double * __restrict );
					
				

				// Store aligned
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				void store_a(double * __restrict ,
					     double * __restrict ) const;
					
				

				// Store unaligned
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				void store_u(double * __restrict ,
				             double * __restrict ) const; 
					
				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 __ATTR_VECTORCALL__ 
				void stream(double * __restrict ,
					    double * __restrict ) const;
					
                                  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				double extract_f64(const int32_t) const;
					
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				std::pair<double, double> 
				extract_2f64(const int32_t ,
					     const int32_t );
					
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM16c8 & insert_1f64(const int32_t ,
						          const double );
					
                                  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM16c8 & insert_2f64(const int32_t ,
						          const int32_t ,
						          const double ,
						          const double );
					

				// Length of 16 doubles
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				void concatenate_a(double * __restrict ) const;
					

				// Length of 16 doubles
				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				void concatenate_u(double * __restrict ) const;

				   __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				 ZMM16c8 & partial_loadu(const double const * __restrict,
							    const int32_t,
							    const double const * __restrict,
							    const int32_t );

				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 			    
				 ZMM16c8 & partial_loada(const double const * __restrict,
							    const int32_t,
							    const double const * __restrict,
							    const int32_t);

				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 			    
				void partial_storeu(double * __restrict,
					            const int32_t,
						    double * __restrict,
					            const int32_t);

				  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 		    
				void partial_storea(double * __restrict,
						    const int32_t,
						    double * __restrict,
					            const int32_t);
					
                                  __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 
				 ZMM16c8 & expand(const  ZMM16c8,
						     const __mmask8); 

				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				  __ATTR_VECTORCALL__ 		     
				 ZMM16c8 & expand_load(const  ZMM16c8,
						          const __mmask8 ,
						          const double * __restrict ,
						          const double * __restrict );


                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				 ZMM16c8 & permute(const __mmask8 ,
						      const int32_t);
					
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d re_low2() const; 

				
                                 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d re_hi2() const; 

				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d im_low2() const;

				 __ATTR_ALWAYS_INLINE__
				 __ATTR_HOT__
		                 __ATTR_ALIGN__(32)
				__m256d im_hi2() const;

				

				 ZMM16c8 & operator=(const  ZMM16c8);

				
		};

	


		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline  ZMM16c8
		        conj(const  ZMM16c8);

	         __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 		
		static inline  ZMM16c8
		        polar(const __m512d,
			      const __m512d);

	        __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 		      
		static inline __m512d
		        carg(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	
		static inline __m512d 
		        carg(const double,
			     const double);

		 __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline  ZMM16c8 
			csin(const  ZMM16c8);

		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 		
		static inline  ZMM16c8
		        csin(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__     
		static inline  ZMM16c8
		        csinh(const  ZMM16c8);

		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline  ZMM16c8
		        csinh(const double,
			      const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	      
		static inline  ZMM16c8
			ccos(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	
		static inline  ZMM16c8
		        ccos(const double,
			     const double);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	     
		static inline  ZMM16c8
		        ccosh(const  ZMM16c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	
		static inline  ZMM16c8
		        ccosh(const double,
			      const double);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	      
		static inline  ZMM16c8
			cexp(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM16c8
		        cexp(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline __m512d
			cabs(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline __m512d
		        cabs(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline  ZMM16c8
		        cpowi(const ZMM16c8,
			      const double );

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		      
		static inline ZMM16c8
		        clog(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM16c8
		        clog(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline  ZMM16c8
		        csqrt(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM16c8
		        csqrt(const double,
			      const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		      
		static inline  ZMM16c8
		        ctan(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM16c8
		        ctan(const double,
			     const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		     
		static inline  ZMM16c8
		        ctanh(const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		
		static inline  ZMM16c8
		        ctanh(const double,
			      const double);
			      
			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	
		static inline  ZMM16c8
		select(const  ZMM16c8 ,
		       const  ZMM16c8,
		       const __mmask8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	       
		static inline  ZMM16c8
		cdiv_smith(const ZMM16c8 ,
			   const  ZMM16c8 );

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		   
		static inline  ZMM16c8
		operator+(const  ZMM16c8,
			  const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		  
		static inline  ZMM16c8
		operator+(const  ZMM16c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+(const __m512d,
			  const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+(const  ZMM16c8,
		          const double );
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+(const double ,
		          const  ZMM16c8); 
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+=( ZMM16c8,
		           const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+=( ZMM16c8,
			   const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+=(const __m512d,
			    ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+=( ZMM16c8 ,
		          const double);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator+=(const double,
			    ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-(const  ZMM16c8,
			  const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-(const  ZMM16c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline ZMM16c8
		operator-(const __m512d,
			  const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-(const  ZMM16c8,
			  const double );
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-(const double ,
			  const ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-( ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-=( ZMM16c8,
			   const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-=( ZMM16c8,
			   const __m512d);
		  
		static inline  ZMM16c8
		operator-=(const __m512d,
			    ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-=( ZMM16c8,
			   const double );
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator-=(const double,
		            ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline ZMM16c8
		operator*(const  ZMM16c8,
			  const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator*(const  ZMM16c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator*(const __m512d,
			  const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator*(const  ZMM16c8,
		          const double);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		  
		static inline  ZMM16c8
		operator*(const double,
		          const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		  
		static inline  ZMM16c8
		operator*=( ZMM16c8,
			   const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		   
		static inline  ZMM16c8
		operator*=( ZMM16c8,
			   const __m512d);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__		   
		static inline  ZMM16c8
		operator*=(const __m512d,
			    ZMM16c8);

		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__			    
	        static inline  ZMM16c8
		operator*=( ZMM16c8,
			   const double);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator*=(const double,
		           ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator/(const  ZMM16c8,
		          const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator/(const ZMM16c8,
			  const __m512d);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator/(const __m512d,
			  const  ZMM16c8);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline  ZMM16c8
		operator/(const  ZMM16c8,
			  const double);
		__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline ZMM16c8
		operator/(const double,
			  const  ZMM16c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	  
		static inline  ZMM16c8
		operator/=( ZMM16c8,
			   const  ZMM16c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline  ZMM16c8
		operator/=( ZMM16c8,
			   const __m512d);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline  ZMM16c8
		operator/=(const __m512d,
			    ZMM16c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	    
	        static inline  ZMM16c8
		operator/=( ZMM16c8,
			   const double);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline  ZMM16c8
		operator/=(const double,
		            ZMM16c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	    
		static inline  ZMM16c8
		operator~( ZMM16c8);

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__
		static inline
			std::pair<__mmask8, __mmask8>
		operator==(const  ZMM16c8,
			   const  ZMM16c8);

		

				__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator==(const  ZMM16c8,
		           std::complex<double>);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator==(const std::complex<double>,
		           const  ZMM16c8)
		           

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline 
			std::pair<__mmask8, __mmask8>
		operator!=(const  ZMM16c8,
		           const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator!=(const  ZMM16c8,
		           const std::complex<double>);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline
		        std::pair<__mmask8,__mmask8>
		operator!=(const std::complex<double>,
		           const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline 
			std::pair<__mmask8, __mmask8>
		operator>(const  ZMM16c8,
			  const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline
                        std::pair<__mmask8,__mmask8>
		operator>(const  ZMM16c8,
		          const std::complex<double>);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline	  
			std::pair<__mmask8,__mmask8>
		operator>(const std::complex<double>,
		          const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline 
			std::pair<__mmask8, __mmask8>
		operator<(const  ZMM16c8,
			  const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	  
		static inline	  
                        std::pair<__mmask8,__mmask8>
		operator<(const  ZMM16c8,
		          std::complex<double>);

			  	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
                static inline
		        std::pair<__mmask8,__mmask8>
		operator<(const std::complex<double>,
		          const  ZMM16c8);

			  	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline 
			std::pair<__mmask8, __mmask8>
		operator>=(const  ZMM16c8,
		           const  ZMM16c8);

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator>=(const  ZMM16c8,
		           const std::complex<double>);

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator>=(const std::complex<double>,
		           const  ZMM16c8);

			__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 	   
		static inline 
			std::pair<__mmask8, __mmask8>
	        operator<=(const  ZMM16c8,
			   const  ZMM16c8);

			   	__ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator<=(const  ZMM16c8,
		           const std::complex<double>);
			   
                __ATTR_ALWAYS_INLINE__
		__ATTR_HOT__
		__ATTR_ALIGN__(32)
	        __ATTR_VECTORCALL__ 
		static inline
		        std::pair<__mmask8,__mmask8>
		operator<=(const std::complex<double>,
		           const  ZMM16c8);

#include "GMS_avx512c8f64.inl"	

	}
}



#endif /*__GMS_AVX512C8F64_H__*/
