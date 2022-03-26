
#ifndef __GMS_AVXCOMPLEX_H__
#define __GMS_AVXCOMPLEX_H__ 081020191905

namespace file_info{

const unsigned int gGMS_AVXCOMPLEX_MAJOR = 1;

const unsigned int gGMS_AVXCOMPLEX_MINOR = 0;

const unsigned int gGMS_AVXCOMPLEX_MICRO = 1;

const unsigned int gGMS_AVXCOMPLEX_FULLVER = 
	1000U*gGMS_AVXCOMPLEX_MAJOR + 100U*gGMS_AVXCOMPLEX_MINOR + 10U*gGMS_AVXCOMPLEX_MICRO;

const char * const pgGMS_AVXCOMPLEX_CREATE_DATE = "08-10-2019 19:05 +00200 (TUE 08 OCT 2019 GMT+2)";


const char * const pgGMS_AVXCOMPLEX_BUILD_DATE =  __DATE__":"__TIME__;

const char * const pgGMS_AVXCOMPLEX_AUTHOR  = "Programmer: Bernard Gingold e-mail: beniekg@gmail.com";

const char * const pgGMS_AVXCOMPLEX_DESCRIPT = "AVX manual vectorization of complex-valued arrays.";

}






#include <cstdint>

#include "GMS_config.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_AVXCOMPLEX_NT_STORES)
#define USE_AVXCOMPLEX_NT_STORES 0
#endif

namespace gms {
	namespace math{

		
		struct AVXVCData{


	         

		  __attribute__((align_value(8))) double* __restrict                  m_Re;
		  __attribute__((align_value(8))) double* __restrict                  m_Im;
                  int32_t				                               m_nsize;

#if (USE_STRUCT_PADDING) == 1
			PAD_TO_ALIGNED(4,0,4)
#endif
#if (USE_STRUCT_PADDING) == 1
		   PAD_TO_ALIGNED(8,1,8)
#endif
		};

		__attribute__((align(64)))
				struct AVXComplex1D{
		
				 AVXVCData data;

				/*
					Default Constructor 
				*/
				AVXVComplex1D();

				/*
					Single argument default explicit Constructor
					Initialization of array members to default values i.e. (NaN)
				*/
				explicit AVXVComplex1D(const int32_t);

				AVXVComplex1D(const double * __restrict,
					      const double * __restrict,
					      const int32_t);


				
				/*
					Copy-Constructor implements deep-copy semantics.
				*/
				AVXVComplex1D(const AVXVComplex1D &);

				
				/*
					Move-Constructor implements shallow-copy semantics.
				*/
				AVXVComplex1D(AVXVComplex1D &&) noexcept(true);

				/*
					Class Destructor
				*/
				~AVXVComplex1D() noexcept(true);

				/*
				    Class member operators
				*/

				/*
				    Copy-assignment operator
				    (deep-copy semantics)
				*/
				AVXVComplex1D & operator=(const AVXVComplex1D &);

				/*
					Move-assignment operator
					(shallow-copy semantics)
				*/
				AVXVComplex1D & operator=(AVXVComplex1D &&);

			
		};		


		std::ostream &
		operator<<(std::ostream &,
			   const AVXVComplex1D &);
		
		AVXVComplex1D
		operator+(const AVXVComplex1D &,
			  const AVXVComplex1D &);

		AVXVComplex1D
		operator+(const AVXVComplex1D &,
			  const double * __restrict);

	    AVXVComplex1D
		operator-(const AVXVComplex1D &,
			  const AVXVComplex1D &);

		AVXVComplex1D
		operator-(const AVXVComplex1D &,
			  const double * __restrict);
			
		AVXVComplex1D
		operator*(const AVXVComplex1D &,
			  const AVXVComplex1D &);

		AVXVComplex1D
		operator*(const AVXVComplex1D &,
			  const double * __restrict);

		AVXVComplex1D
		operator/(const AVXVComplex1D &,
			  const AVXVComplex1D &);

		AVXVComplex1D
		operator/(const AVXVComplex1D &,
			  const double * __restrict);

		AVXVComplex1D
		operator==(const AVXVComplex1D &,
			   const AVXVComplex1D &);

		AVXVComplex1D
		operator!=(const AVXVComplex1D &,
			   const AVXVComplex1D &);

		
		void v256cnormalize_product( AVXVComplex1D &, 
					     const AVXVComplex1D &, 
					     const AVXVComplex1D &,
					     const bool) noexcept(true);
		
		
		void v256cmean_product(std::complex<double> &,
				       const AVXVComplex1D &v1,
				       const AVXVComplex1D &v2);

	  
		void v256cmean_quotient(std::complex<double> &,
					const AVXVComplex1D &,
					const AVXVComplex1D &);

	  
		void v256cconj_product(AVXVComplex1D &,
				       const AVXVComplex1D &,
				       const AVXVComplex1D &,
				       const bool)noexcept(true);

		void v256cnorm_conjprod(AVXVComplex1D &,
					const AVXVComplex1D &,
					const AVXVComplex1D &,
					const bool)noexcept(true);

		void v256cmean_conjprod(std::complex<double> &,
				        const AVXVComplex1D &,
					const AVXVComplex1D &);

		void v256c_arithmean(std::complex<double> &,
				     const AVXVComplex1D &);

		void v256c_normalize(AVXVComplex1D &,
				     const AVXVComplex1D &,
				     const AVXVComplex1D &,
				     const bool);

		void v256c_magnitude(AVXVComplex1D &,
				     const AVXVComplex1D &,
				     const AVXVComplex1D &);

	}
}



#endif /*__GMS_AVXCOMPLEX_H__*/
