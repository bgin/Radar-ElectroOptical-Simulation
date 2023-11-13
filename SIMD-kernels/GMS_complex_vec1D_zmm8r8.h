
#ifndef __GMS_COMPLEX_VEC_1D_ZMM8R8_H__
#define __GMS_COMPLEX_VEC_1D_ZMM8R8_H__

namespace file_info {

	const unsigned int GMS_COMPLEX_VEC_1D_ZMM8R8_MAJOR = 1U;

	const unsigned int GMS_COMPLEX_VEC_1D_ZMM8R8_MINOR = 1U;

	const unsigned int GMS_COMPLEX_VEC_1D_ZMM8R8_MICRO = 0U;

	const unsigned int GMS_COMPLEX_VEC_1D_ZMM8R8_FULLVER = 
	    1000U*GMS_COMPLEX_VEC_1D_ZMM8R8_MAJOR + 100U*GMS_COMPLEX_VEC_1D_ZMM8R8_MINOR + 10U*GMS_COMPLEX_VEC_1D_ZMM8R8_MICRO;

	const  char * const  GMS_COMPLEX_VEC_1D_ZMM8R8_CREATE_DATE = "10-10-2019 19:13 +00200 (FRI 10 OCT 2019 GMT+2)";

	const  char * const  GMS_COMPLEX_VEC_1D_ZMM8R8_BUILD_DATE = __DATE__ ":" __TIME__;

	const  char * const  GMS_COMPLEX_VEC_1D_ZMM8R8_AUTHOR     = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const  char * const  GMS_COMPLEX_VEC_1D_ZMM8R8_DESCRIPT   =  "AVX512(512-bit) manual vectorization of complex-valued arrays.";

}





#include <immintrin.h>
#include <complex>
#include <cstdint>
#include <array>
#include "GMS_config.h"
#include "GMS_common.h"

//#include "GMS_constants.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_AVX512COMPLEX_NT_STORES)
#define USE_AVX512COMPLEX_NT_STORES 0
#endif

namespace gms {
	namespace math{




		 

		 __ATTR_ALIGN__(64) struct CV1DZMM8r8{



		         double* __restrict m_Re;
		         double* __restrict m_Im;
		         int32_t   m_nsize;
#if (USE_STRUCT_PADDING) == 1
			 PAD_TO(0,44)
#endif
			//
			//	Default Constructor.
			//
			CV1DZMM8r8();

			//
			// explicit one-arg Constructor which performs default
			// object initialization by allocating member arrays and
			// initialize then by NaN value. This has no mathematical
			// sense, but it mimmicks default initialization.
			// 
			explicit CV1DZMM8r8(const int32_t);

			// Build from two pointers.
			CV1DZMM8r8(const double * __restrict,
					 const double * __restrict,
					 const int32_t);

		    
			//
			// Copy-Constructor implements deep copy semantics.
			//
			CV1DZMM8r8(const CV1DZMM8r8 &);

			//
			// Move-Constructor implements shallow copy semantics
			//
			CV1DZMM8r8(CV1DZMM8r8 &&) noexcept;

			//
			// Class Destructor.
			//
			~CV1DZMM8r8();

			//
			// ***Class member and friend operators***
			//

			//
			// Copy-assignment implements deep-copy semantics
			//
			CV1DZMM8r8 & operator=(const CV1DZMM8r8 &) = delete;

			//
			// Copy-assignment implements shallow-copy semantics
			//
			CV1DZMM8r8 & operator=(CV1DZMM8r8 &&) noexcept(true);

	};		
			
			

			
			//********* Global Operators *************//
			

			

			


		 std::ostream &
		  operator<<(std::ostream &, 
		             const CV1DZMM8r8 &);


		CV1DZMM8r8   
		operator+(const CV1DZMM8r8 &,
			  const CV1DZMM8r8 &);

		CV1DZMM8r8
		operator+(const CV1DZMM8r8 &,
			  const double * __restrict);


		CV1DZMM8r8   
		operator-(const CV1DZMM8r8 &x,
			  const CV1DZMM8r8 &y);


		CV1DZMM8r8
		operator-(const CV1DZMM8r8 &,
			  const double * __restrict);
			

		CV1DZMM8r8
		operator*(const CV1DZMM8r8 &,
			  const CV1DZMM8r8 &);
			
			
		CV1DZMM8r8
		operator*(const CV1DZMM8r8 &,
			  const double * __restrict);
		

		CV1DZMM8r8	
		operator/(const CV1DZMM8r8 &,
			  const CV1DZMM8r8 &);

		CV1DZMM8r8
		operator/(const CV1DZMM8r8 &,
			  const double * __restrict);
		
		template< size_t N>
		std::pair<std::array<__mmask8,N>,
		std::array<__mmask8,N>>
	    operator==(const CV1DZMM8r8 &x,
		       const CV1DZMM8r8 &y) {
				  
				   using namespace gms::math::constants;
			if (x.m_nsize != y.m_nsize) {
				return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
			}
			int32_t i;
			int32_t k1 = 0, k2 = 0;
			std::pair<std::array<__mmask8, N>,
				      std::array<__mmask8, N >> ret_val; 
			for (i = 0; i != ROUND_TO_EIGHT(x.m_nsize); i += 8) {
				const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i]);
				const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i]);
				ret_val.first[++k1] = _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_EQ_OQ);
				const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i]);
				const __m512d zmm3 = _mm512_load_pd(&y.m_Im[i]);
				ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_EQ_OQ);
			}
			 unsigned char t1 = 0x00, t2 = 0x00;
			 k1 += 1;
			 k2 += 1;
			for (; i != x.m_nsize; ++i) {
				if (approximately_equalf64(x.m_Re[i],
					              y.m_Re[i], std::numeric_limits<double>::epsilon())) {
					t1 |= 1 << i;
				}
				if (approximately_equalf64(x.m_Im[i],
					             y.m_Im[i], std::numeric_limits<double>::epsilon())) {
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
	    operator!=(const CV1DZMM8r8 &x,
		       const CV1DZMM8r8 &y) {
			using namespace gms::common;
			using namespace gms::math::constants;
			if (x.m_nsize != y.m_nsize) {
				return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
			}
			int32_t i;
			int32_t k1 = 0, k2 = 0;
			std::pair<std::array<__mmask8,N>,std::array<__mmask8,N>> ret_val;
			for (i = 0; i != ROUND_TO_EIGHT(x.m_nsize, 8); i += 8) {
				const __m512d zmm0   = _mm512_load_pd(&x.m_Re[i]);
				const __m512d zmm1   = _mm512_load_pd(&y.m_Re[i]);
				ret_val.first[++k1]  = _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_NEQ_OQ);
				const __m512 zmm2    = _mm512_load_pd(&x.m_Im[i]);
				const __m512 zmm3    = _mm512_load_pd(&y.m_Im[i]);
				ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_NEQ_OQ);
			}
			unsigned char t1 = 0x00, t2 = 0x00;
			k1 += 1;
			k2 += 1;
			for (; i != x.m_nsize; ++i) {
				if (!approximately_equalf64(x.m_Re[i],
					               y.m_Re[i],std::numeric_limits<double>::epsilon())) {
					 t1 |= 1 << i;
				}
				if (!approximately_equalf64(x.m_Im[i],
					y.m_Im[i], std::numeric_limits<double>::epsilon())) {
					t2 |= 1 << i;
				}
			}
			ret_val.first[k1] = t1;
			ret_val.second[k2] = t2;
			return (ret_val);
	   }
		
		
		// 
		// 
	    // Normalization of vector CV1DZMM8r8 product.
		// All three arguments must have the same size.
		// No error checking is made on the inputs.
		//
		void v512cnormalize_product(CV1DZMM8r8 &, 
					    const CV1DZMM8r8 &,
					    const CV1DZMM8r8 &,
					    const bool) noexcept(true);
		//
		//  Arithmetic mean product of complex-valued vectors
		//  Two CV1DZMM8r8 arguments must have the same length.
		//  No error checking is made on the inputs.
		// 
		void v512cmean_product(std::complex<double> &, 
				       const CV1DZMM8r8  &,
				       const CV1DZMM8r8  &) noexcept(true);
							   

	    //
		// Mean complex quotient of CV1DZMM8r8
		// Two CV1DZMM8r8 arguments must have the same length.
		// No error checking is made on the inputs.
		// 
		void v512cmean_quotient(std::complex<double> &,
					const CV1DZMM8r8 &,
					const CV1DZMM8r8 &) noexcept(true);

		//
		// Conjugate product of CV1DZMM8r8 
		// Three CV1DZMM8r8 arguments must have the same length.
		// No error checking is made on the inputs.
		//
		void v512cconj_product(CV1DZMM8r8    &,
				       const CV1DZMM8r8 &,
				       const CV1DZMM8r8 &,
				       const bool  ) noexcept(true);

		//
		// Normalized conjugate product of CV1DZMM8r8 
		// Three CV1DZMM8r8 arguments must have the same length.
		// No error checking is made on the inputs.
		//
		void v512cnorm_conjprod(CV1DZMM8r8 &,
					const CV1DZMM8r8 &,
					const CV1DZMM8r8 &,
					const bool) noexcept(true);

		//
		//	Mean conjugate product of CV1DZMM8r8
		//  Two CV1DZMM8r8 arguments must have the sane length.
		//	No error checking is made on the inputs.
		//
		void v512cmean_conjprod(std::complex<double> &,
					const CV1DZMM8r8 &,
					const CV1DZMM8r8 &) noexcept(true);

		//
		//	Artithmetic mean of complex series (CV1DZMM8r8 vectors). 
		//  No error checking is made on the inputs.
		//	Bear in mind that AVX3 reduction will probably
		//  incur some penalty on code execution because
		//  of horizontal addition.
		//
		void v512c_arithmean(std::complex<double> &,
				     const CV1DZMM8r8 &) noexcept(true);

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
		void v512c_normalize(CV1DZMM8r8    &,
				     const CV1DZMM8r8 &,
				     const CV1DZMM8r8 &,
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
				     const CV1DZMM8r8 &,
				     const CV1DZMM8r8 &) noexcept(true);
	}
}



#endif /*__GMS_complex_vec1D_zmm8r8_H__*/
