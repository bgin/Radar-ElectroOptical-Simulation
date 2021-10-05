
#ifndef __GMS_AVX512COMPLEX_H__
#define __GMS_AVX512COMPLEX_H__ 101020191913

namespace file_info {

	const unsigned int gGMS_AVX512COMPLEX_MAJOR = 1;

	const unsigned int gGMS_AVX512COMPLEX_MINOR = 0;

	const unsigned int gGMS_AVX512COMPLEX_MICRO = 0;

	const unsigned int gGMS_AVX512COMPLEX_FULLVER = 
	    1000U*gGMS_AVX512COMPLEX_MAJOR + 100U*gGMS_AVX512COMPLEX_MINOR + 10U*gGMS_AVX512COMPLEX_MICRO;

	const  char * const  pgGMS_AVX512COMPLEX_CREATE_DATE = "10-10-2019 19:13 +00200 (FRI 10 OCT 2019 GMT+2)";

	const  char * const  pgGMS_AVX512COMPLEX_BUILD_DATE = __DATE__":"__TIME__;

	const  char * const  pgGMS_AVX512COMPLEX_AUTHOR     = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const  char * const  pgGMS_AVX512COMPLEX_DESCRIPT   =  "AVX512(512-bit) manual vectorization of complex-valued arrays.";

}



//
// Warning:
//				Include these files if and only if you have 
//				CPU and/or Accelarator i.e Xeon Phi which supports AVX3 ISA,
//				otherwise remove these files from compilation.
//


#include <complex>
#include <cstdint>
#include <array>


    #include "GMS_config.h"
    #include "GMS_common.h"
    #include "GMS_constants.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_AVX512COMPLEX_NT_STORES)
#define USE_AVX512COMPLEX_NT_STORES 0
#endif

namespace gms {
	namespace math{



	  
	         struct AVX512VCData{



		         __attribute__((align(8))) double* __restrict m_Re;
		         __attribute__((align(8))) double* __restrict m_Im;
		         __attribute__((align(4))) int32_t __restrict m_nsize;


#if (USE_STRUCT_PADDING) == 1
			 PAD_TO_ALIGNED(4,0,4)
#endif			

#if (USE_STRUCT_PADDING) == 1
			 PAD_TO_ALIGNED(8,1,8)
#endif
		};

		 


		 __attribute__((align64))) struct AVX512VComplex1D{



		        AVX512VCData data;

			//
			//	Default Constructor.
			//
			AVX512VComplex1D();

			//
			// explicit one-arg Constructor which performs default
			// object initialization by allocating member arrays and
			// initialize then by NaN value. This has no mathematical
			// sense, but it mimmicks default initialization.
			// 
			explicit AVX512VComplex1D(const int32_t);

			// Build from two pointers.
			AVX512VComplex1D(const double * __restrict,
					 const double * __restrict,
					 const int32_t);

		    
			//
			// Copy-Constructor implements deep copy semantics.
			//
			AVX512VComplex1D(const AVX512VComplex1D &);

			//
			// Move-Constructor implements shallow copy semantics
			//
			AVX512VComplex1D(AVX512VComplex1D &&) noexcept;

			//
			// Class Destructor.
			//
			~AVX512VComplex1D();

			//
			// ***Class member and friend operators***
			//

			//
			// Copy-assignment implements deep-copy semantics
			//
			AVX512VComplex1D & operator=(const AVX512VComplex1D &);

			//
			// Copy-assignment implements shallow-copy semantics
			//
			AVX512VComplex1D & operator=(AVX512VComplex1D &&) noexcept(true);

	};		
			
			

			
			//********* Global Operators *************//
			

			

			


		 std::ostream &
		  operator<<(std::ostream &, 
		             const AVX512VComplex1D &);


		AVX512VComplex1D   
		operator+(const AVX512VComplex1D &,
			  const AVX512VComplex1D &);

		AVX512VComplex1D
		operator+(const AVX512VComplex1D &,
			  const double * __restrict);


		AVX512VComplex1D   
		operator-(const AVX512VComplex1D &x,
			  const AVX512VComplex1D &y);


		AVX512VComplex1D
		operator-(const AVX512VComplex1D &,
			  const double * __restrict);
			

		AVX512VComplex1D
		operator*(const AVX512VComplex1D &,
			  const AVX512VComplex1D &);
			
			
		AVX512VComplex1D
		operator*(const AVX512VComplex1D &,
			  const double * __restrict);
		

		AVX512VComplex1D	
		operator/(const AVX512VComplex1D &,
			  const AVX512VComplex1D &);

		AVX512VComplex1D
		operator/(const AVX512VComplex1D &,
			  const double * __restrict);
		
		template< size_t N>
		std::pair<std::array<__mmask8,N>,
		std::array<__mmask8,N>>
	    operator==(const AVX512VComplex1D &x,
		       const AVX512VComplex1D &y) {
				   using namespace gms::common;
				   using namespace gms::math::constants;
			if (x.data.m_nsize != y.data.m_nsize) {
				return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
			}
			int32_t i;
			int32_t k1 = 0, k2 = 0;
			std::pair<std::array<__mmask8, N>,
				      std::array<__mmask8, N >> ret_val; 
			for (i = 0; i != ROUND_TO_EIGHT(x.data.m_nsize); i += 8) {
				const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i]);
				const __m512d zmm1 = _mm512_load_pd(&y.data.m_Re[i]);
				ret_val.first[++k1] = _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_EQ_OQ);
				const __m512d zmm2 = _mm512_load_pd(&x.data.m_Im[i]);
				const __m512d zmm3 = _mm512_load_pd(&y.data.m_Im[i]);
				ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_EQ_OQ);
			}
			 unsigned char t1 = 0x00, t2 = 0x00;
			 k1 += 1;
			 k2 += 1;
			for (; i != x.data.m_nsize; ++i) {
				if (approximately_equalf64(x.data.m_Re[i],
					              y.data.m_Re[i], DEPS)) {
					t1 |= 1 << i;
				}
				if (approximately_equalf64(x.data.m_Im[i],
					             y.data.m_Im[i], DEPS)) {
					t2 |= 1 << i
				}
			}
			ret_val.first[k1]  = t1;
			ret_val.second[k2] = t2;
			return (ret_val); // Move Ctor should kick in if compiler will be smart enough
						      // or at least fast memcpy of primitive type
	   }

		template< size_t N>
		std::pair<std::array<__mmask8,N>,
		std::array<__mmask8,N>>
	    operator!=(const AVX512VComplex1D &x,
		       const AVX512VComplex1D &y) {
			using namespace gms::common;
			using namespace gms::math::constants;
			if (x.data.m_nsize != y.data.m_nsize) {
				return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
			}
			int32_t i;
			int32_t k1 = 0, k2 = 0;
			std::pair<std::array<__mmask8,N>,std::array<__mmask8,N>> ret_val;
			for (i = 0; i != ROUND_TO_EIGHT(x.data.m_nsize, 8); i += 8) {
				const __m512d zmm0   = _mm512_load_pd(&x.data.m_Re[i]);
				const __m512d zmm1   = _mm512_load_pd(&y.data.m_Re[i]);
				ret_val.first[++k1]  = _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_NEQ_OQ);
				const __m512 zmm2    = _mm512_load_pd(&x.data.m_Im[i]);
				const __m512 zmm3    = _mm512_load_pd(&y.data.m_Im[i]);
				ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_NEQ_OQ);
			}
			unsigned char t1 = 0x00, t2 = 0x00;
			k1 += 1;
			k2 += 1;
			for (; i != x.data.m_nsize; ++i) {
				if (!approximately_equalf64(x.data.m_Re[i],
					               y.data.m_Re[i], DEPS)) {
					 t1 |= 1 << i;
				}
				if (!approximately_equalf64(x.data.m_Im[i],
					y.data.m_Im[i], DEPS)) {
					t2 |= 1 << i;
				}
			}
			ret_val.first[k1] = t1;
			ret_val.second[k2] = t2;
			return (ret_val);
	   }
		
		
		// 
		// 
	    // Normalization of vector AVX512VComplex1D product.
		// All three arguments must have the same size.
		// No error checking is made on the inputs.
		//
		void v512cnormalize_product(AVX512VComplex1D &, 
					    const AVX512VComplex1D &,
					    const AVX512VComplex1D &,
					    const bool) noexcept(true);
		//
		//  Arithmetic mean product of complex-valued vectors
		//  Two AVX512VComplex1D arguments must have the same length.
		//  No error checking is made on the inputs.
		// 
		void v512cmean_product(std::complex<double> &, 
				       const AVX512VComplex1D  &,
				       const AVX512VComplex1D  &) noexcept(true);
							   

	    //
		// Mean complex quotient of AVX512VComplex1D
		// Two AVX512VComplex1D arguments must have the same length.
		// No error checking is made on the inputs.
		// 
		void v512cmean_quotient(std::complex<double> &,
					const AVX512VComplex1D &,
					const AVX512VComplex1D &) noexcept(true);

		//
		// Conjugate product of AVX512VComplex1D 
		// Three AVX512VComplex1D arguments must have the same length.
		// No error checking is made on the inputs.
		//
		void v512cconj_product(AVX512VComplex1D    &,
				       const AVX512VComplex1D &,
				       const AVX512VComplex1D &,
				       const bool  ) noexcept(true);

		//
		// Normalized conjugate product of AVX512VComplex1D 
		// Three AVX512VComplex1D arguments must have the same length.
		// No error checking is made on the inputs.
		//
		void v512cnorm_conjprod(AVX512VComplex1D &,
					const AVX512VComplex1D &,
					const AVX512VComplex1D &,
					const bool) noexcept(true);

		//
		//	Mean conjugate product of AVX512VComplex1D
		//  Two AVX512VComplex1D arguments must have the sane length.
		//	No error checking is made on the inputs.
		//
		void v512cmean_conjprod(std::complex<double> &,
					const AVX512VComplex1D &,
					const AVX512VComplex1D &) noexcept(true);

		//
		//	Artithmetic mean of complex series (AVX512VComplex1D vectors). 
		//  No error checking is made on the inputs.
		//	Bear in mind that AVX3 reduction will probably
		//  incur some penalty on code execution because
		//  of horizontal addition.
		//
		void v512c_arithmean(std::complex<double> &,
				     const AVX512VComplex1D &) noexcept(true);

		//
		//	Normalize complex vector
		//	In this case a temporary argument, which
		//  is an exact copy of argument being normalized
		//  is passed to this procedure.
		//  This is done in order to trigger HW prefetcher.
		//  An exact copy should be constructed by calling
		//  class Move Constructor.
		//	All vectors must have the same length.
		//
		void v512c_normalize(AVX512VComplex1D    &,
				     const AVX512VComplex1D &,
				     const AVX512VComplex1D &,
				     const bool ) noexcept(true);

		//
		//	Complex vector magnitude.
		//	In this case a temporary argument, which
		//  is an exact copy of argument being normalized
		//  is passed to this procedure.
		//  This is done in order to trigger HW prefetcher.
		//  An exact copy should be constructed by calling
		//  class Move Constructor.
		//	All 2 vectors and array (double) must have the same length.
		//
		void v512c_magnitude(double * __restrict,
				     const AVX512VComplex1D &,
				     const AVX512VComplex1D &) noexcept(true);
	}
}



#endif /*__GMS_AVX512COMPLEX_H__*/
