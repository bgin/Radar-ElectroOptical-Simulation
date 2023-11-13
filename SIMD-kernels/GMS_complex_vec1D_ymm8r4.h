
#ifndef __GMS_COMPLEX_VEC1D_YMM8R4_H__
#define __GMS_COMPLEX_VEC1D_YMM8R4_H__ 091020191905

namespace file_info{

  
const unsigned int GMS_COMPLEX_VEC1D_YMM8R4_MAJOR = 1U;

const unsigned int GMS_COMPLEX_VEC1D_YMM8R4_MINOR = 1U;

const unsigned int GMS_COMPLEX_VEC1D_YMM8R4_MICRO = 0U;

const unsigned int GMS_COMPLEX_VEC1D_YMM8R4_FULLVER = 
	1000U*GMS_COMPLEX_VEC1D_YMM8R4_MAJOR + 100U*GMS_COMPLEX_VEC1D_YMM8R4_MINOR + 10U*GMS_COMPLEX_VEC1D_YMM8R4_MICRO;

const char * const GMS_COMPLEX_VEC1D_YMM8R4_CREATE_DATE = "09-10-2019 19:05 +00200 (WED 09 OCT 2019 GMT+2)";


const char * const GMS_COMPLEX_VEC1D_YMM8R4_BUILD_DATE =  __DATE__ ":" __TIME__;

const char * const GMS_COMPLEX_VEC1D_YMM8R4_AUTHOR  = "Programmer: Bernard Gingold e-mail: beniekg@gmail.com";

const char * const GMS_COMPLEX_VEC1D_YMM8R4_DESCRIPT = "AVX manual vectorization of complex-valued arrays.";

}





#include <cstdint>
#include "GMS_config.h"
// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_COMPLEX_VEC1D_YMM8R4_NT_STORES)
#define USE_COMPLEX_VEC1D_YMM8R4_NT_STORES 0
#endif

namespace gms {
	namespace math{

		
	  

	 __ATTR_ALIGN__(64) struct CV1DYMM8r4{
		
				
				 float* __restrict  m_Re;
		                 float* __restrict  m_Im;
                                 int32_t            m_nsize;
#if (USE_STRUCT_PADDING) == 1
		   PAD_TO(0,44)
#endif

				/*
					Default Constructor 
				*/
				CV1DYMM8r4();

				/*
					Single argument default explicit Constructor
					Initialization of array members to default values i.e. (NaN)
				*/
				explicit CV1DYMM8r4(const int32_t);

				CV1DYMM8r4(   const float * __restrict,
					      const float * __restrict,
					      const int32_t);


				
				/*
					Copy-Constructor implements deep-copy semantics.
				*/
				CV1DYMM8r4(const CV1DYMM8r4 &);

				
				/*
					Move-Constructor implements shallow-copy semantics.
				*/
				CV1DYMM8r4(CV1DYMM8r4 &&) noexcept(true);

				/*
					Class Destructor
				*/
				~CV1DYMM8r4() noexcept(true);

				/*
				    Class member operators
				*/

				/*
				    Copy-assignment operator
				    (deep-copy semantics)
				*/
				CV1DYMM8r4 & operator=(const CV1DYMM8r4 &) = delete;

				/*
					Move-assignment operator
					(shallow-copy semantics)
				*/
				CV1DYMM8r4 & operator=(CV1DYMM8r4 &&);

			
		};		


		std::ostream &
		operator<<(std::ostream &,
			   const CV1DYMM8r4 &);
		
		CV1DYMM8r4
		operator+(const CV1DYMM8r4 &,
			  const CV1DYMM8r4 &);

		CV1DYMM8r4
		operator+(const CV1DYMM8r4 &,
			  const float * __restrict);

	    CV1DYMM8r4
		operator-(const CV1DYMM8r4 &,
			  const CV1DYMM8r4 &);

		CV1DYMM8r4
		operator-(const CV1DYMM8r4 &,
			  const float * __restrict);
			
		CV1DYMM8r4
		operator*(const CV1DYMM8r4 &,
			  const CV1DYMM8r4 &);

				  
		CV1DYMM8r4
		operator*(const CV1DYMM8r4 &,
			  const float * __restrict);

		CV1DYMM8r4
		operator/(const CV1DYMM8r4 &,
			  const CV1DYMM8r4 &);

		CV1DYMM8r4
		operator/(const CV1DYMM8r4 &,
			  const float * __restrict);

		CV1DYMM8r4
		operator==(const CV1DYMM8r4 &,
			   const CV1DYMM8r4 &);

		CV1DYMM8r4
		operator!=(const CV1DYMM8r4 &,
			   const CV1DYMM8r4 &);

		
		void v256cnormalize_product( CV1DYMM8r4 &, 
					     const CV1DYMM8r4 &, 
					     const CV1DYMM8r4 &,
					     const bool) noexcept(true);
		
		
		void v256cmean_product(std::complex<float> &,
				       const CV1DYMM8r4 &v1,
				       const CV1DYMM8r4 &v2);

	  
		void v256cmean_quotient(std::complex<float> &,
					const CV1DYMM8r4 &,
					const CV1DYMM8r4 &);

	  
		void v256cconj_product(CV1DYMM8r4 &,
				       const CV1DYMM8r4 &,
				       const CV1DYMM8r4 &,
				       const bool)noexcept(true);

		void v256cnorm_conjprod(CV1DYMM8r4 &,
					const CV1DYMM8r4 &,
					const CV1DYMM8r4 &,
					const bool)noexcept(true);

		void v256cmean_conjprod(std::complex<float> &,
				        const CV1DYMM8r4 &,
					const CV1DYMM8r4 &);

		void v256c_arithmean(std::complex<float> &,
				     const CV1DYMM8r4 &);

		void v256c_normalize(CV1DYMM8r4 &,
				     const CV1DYMM8r4 &,
				     const CV1DYMM8r4 &,
				     const bool);

		void v256c_magnitude(CV1DYMM8r4 &,
				     const CV1DYMM8r4 &,
				     const CV1DYMM8r4 &);

	}
}



#endif /*__GMS_COMPLEX_VEC1D_YMM8R4_H__*/
