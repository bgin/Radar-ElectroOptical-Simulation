
#ifndef __GMS_AVXCOMPLEX_H__
#define __GMS_AVXCOMPLEX_H__

namespace file_info{
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif
const unsigned int gGMS_AVXCOMPLEX_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_AVXCOMPLEX_MINOR = gms::common::gVersionInfo.m_VersionMinor;

const unsigned int gGMS_AVXCOMPLEX_MICRO = gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_AVXCOMPLEX_FULLVER = 
	1000U*gGMS_AVXCOMPLEX_MAJOR + 100U*gGMS_AVXCOMPLEX_MINOR + 10U*gGMS_AVXCOMPLEX_MICRO;

const char * const pgGMS_AVXCOMPLEX_CREATE_DATE = "08-10-2019 19:05 +00200 (TUE 08 OCT 2019 GMT+2)";

/*
Set this value to latest build date/time
*/
const char * const pgGMS_AVXCOMPLEX_BUILD_DATE =  "00-00-0000 00:00";

const char * const pgGMS_AVXCOMPLEX_AUTHOR  = "Programmer: Bernard Gingold e-mail: beniekg@gmail.com";

const char * const pgGMS_AVXCOMPLEX_DESCRIPT = "AVX manual vectorization of complex-valued arrays.";

}


/*
	Bernard Gingold copyright notice:
	This file is a part of Guided-Missile-Simulation program
	Copyright(C) 2019 Bernard Gingold
	License : GNU General Public License version 3 or later,
	for additional information check file LICENSE.txt in
	project directory.
*/



#include <cstdint>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#endif
// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_AVXCOMPLEX_NT_STORES)
#define USE_AVXCOMPLEX_NT_STORES 0
#endif

namespace gms {
	namespace math{

		
		struct AVXVCData{
#if defined _WIN64			
			// Real array
			_Field_size_(m_nsize) __declspec(align(64)) double* __restrict m_Re;
			// Imaginary array
			_Field_size_(m_nsize) __declspec(align(64)) double* __restrict m_Im;

			int32_t				                               m_nsize;
#elif defined __linux
		  __attribute__((align_value(64))) double* __restrict                  m_Re;
		  __attribute__((align_value(64))) double* __restrict                  m_Im;
#endif

#if (USE_STRUCT_PADDING) == 1
			PAD_TO_ALIGNED(4,0,4)
#endif
#if (USE_STRUCT_PADDING) == 1
		   PAD_TO_ALIGNED(8,1,8)
#endif
		};
#if defined _WIN64
		__declspec(align(64)) 
				struct  AVXVComplex1D{ // Start of this struct at 64-byte offset.
#elif defined __linux
		__attribute__((align(64)))
				struct AVXComplex1D{
#endif			
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
