
#ifndef __GMS_COMPLEX_VEC1D_YMM4R8_H__
#define __GMS_COMPLEX_VEC1D_YMM4R8_H__

namespace file_info{

  
const unsigned int GMS_COMPLEX_VEC1D_YMM4R8_MAJOR = 1U;

const unsigned int GMS_COMPLEX_VEC1D_YMM4R8_MINOR = 1U;

const unsigned int GMS_COMPLEX_VEC1D_YMM4R8_MICRO = 0U;

const unsigned int GMS_COMPLEX_VEC1D_YMM4R8_FULLVER = 
	1000U*GMS_COMPLEX_VEC1D_YMM4R8_MAJOR + 100U*GMS_COMPLEX_VEC1D_YMM4R8_MINOR + 10U*GMS_COMPLEX_VEC1D_YMM4R8_MICRO;

const char * const GMS_COMPLEX_VEC1D_YMM4R8_CREATE_DATE = "08-10-2019 19:05 +00200 (TUE 08 OCT 2019 GMT+2)";


const char * const GMS_COMPLEX_VEC1D_YMM4R8_BUILD_DATE =  __DATE__ ":" __TIME__;

const char * const GMS_COMPLEX_VEC1D_YMM4R8_AUTHOR  = "Programmer: Bernard Gingold e-mail: beniekg@gmail.com";

const char * const GMS_COMPLEX_VEC1D_YMM4R8_DESCRIPT = "AVX manual vectorization of complex-valued arrays.";

}





#include <cstdint>
#include "GMS_config.h"
// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES)
#define USE_COMPLEX_VEC1D_YMM4R8_NT_STORES 0
#endif

namespace gms {
	namespace math{

		
	    __ATTR_ALIGN__(64) struct DataYMM4r8{

 		  double* __restrict                  m_Re;
		  double* __restrict                  m_Im;
                  int32_t                             m_nsize;
#if (USE_STRUCT_PADDING) == 1
		   PAD_TO(0,44)
#endif
		};

	 __ATTR_ALIGN__(64) struct CV1DYMM4r8{
		
				DataYMM4r8 data;

				/*
					Default Constructor 
				*/
				CV1DYMM4r8();

				/*
					Single argument default explicit Constructor
					Initialization of array members to default values i.e. (NaN)
				*/
				explicit CV1DYMM4r8(const int32_t);

				CV1DYMM4r8(   const double * __restrict,
					      const double * __restrict,
					      const int32_t);


				
				/*
					Copy-Constructor implements deep-copy semantics.
				*/
				CV1DYMM4r8(const CV1DYMM4r8 &);

				
				/*
					Move-Constructor implements shallow-copy semantics.
				*/
				CV1DYMM4r8(CV1DYMM4r8 &&) noexcept(true);

				/*
					Class Destructor
				*/
				~CV1DYMM4r8() noexcept(true);

				/*
				    Class member operators
				*/

				/*
				    Copy-assignment operator
				    (deep-copy semantics)
				*/
				CV1DYMM4r8 & operator=(const CV1DYMM4r8 &) = delete;

				/*
					Move-assignment operator
					(shallow-copy semantics)
				*/
				CV1DYMM4r8 & operator=(CV1DYMM4r8 &&);

			
		};		


		std::ostream &
		operator<<(std::ostream &,
			   const CV1DYMM4r8 &);
		
		CV1DYMM4r8
		operator+(const CV1DYMM4r8 &,
			  const CV1DYMM4r8 &);

		CV1DYMM4r8
		operator+(const CV1DYMM4r8 &,
			  const double * __restrict);

	    CV1DYMM4r8
		operator-(const CV1DYMM4r8 &,
			  const CV1DYMM4r8 &);

		CV1DYMM4r8
		operator-(const CV1DYMM4r8 &,
			  const double * __restrict);
			
		CV1DYMM4r8
		operator*(const CV1DYMM4r8 &,
			  const CV1DYMM4r8 &);

				  
		CV1DYMM4r8
		operator*(const CV1DYMM4r8 &,
			  const double * __restrict);

		CV1DYMM4r8
		operator/(const CV1DYMM4r8 &,
			  const CV1DYMM4r8 &);

		CV1DYMM4r8
		operator/(const CV1DYMM4r8 &,
			  const double * __restrict);

		CV1DYMM4r8
		operator==(const CV1DYMM4r8 &,
			   const CV1DYMM4r8 &);

		CV1DYMM4r8
		operator!=(const CV1DYMM4r8 &,
			   const CV1DYMM4r8 &);

		
		void v256cnormalize_product( CV1DYMM4r8 &, 
					     const CV1DYMM4r8 &, 
					     const CV1DYMM4r8 &,
					     const bool) noexcept(true);
		
		
		void v256cmean_product(std::complex<double> &,
				       const CV1DYMM4r8 &v1,
				       const CV1DYMM4r8 &v2);

	  
		void v256cmean_quotient(std::complex<double> &,
					const CV1DYMM4r8 &,
					const CV1DYMM4r8 &);

	  
		void v256cconj_product(CV1DYMM4r8 &,
				       const CV1DYMM4r8 &,
				       const CV1DYMM4r8 &,
				       const bool)noexcept(true);

		void v256cnorm_conjprod(CV1DYMM4r8 &,
					const CV1DYMM4r8 &,
					const CV1DYMM4r8 &,
					const bool)noexcept(true);

		void v256cmean_conjprod(std::complex<double> &,
				        const CV1DYMM4r8 &,
					const CV1DYMM4r8 &);

		void v256c_arithmean(std::complex<double> &,
				     const CV1DYMM4r8 &);

		void v256c_normalize(CV1DYMM4r8 &,
				     const CV1DYMM4r8 &,
				     const CV1DYMM4r8 &,
				     const bool);

		void v256c_magnitude(CV1DYMM4r8 &,
				     const CV1DYMM4r8 &,
				     const CV1DYMM4r8 &);

	}
}



#endif /*__GMS_COMPLEX_VEC1D_YMM4R8_H__*/
