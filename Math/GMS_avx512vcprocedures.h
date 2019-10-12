
#ifndef __GMS_AVX512VCPROCEDURES_H__
#define __GMS_AVX512VCPROCEDURES_H__

namespace file_info {
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

	const unsigned int gGMS_AVX512PROCEDURES_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_AVX512PROCEDURES_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_AVX512PROCEDURES_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_AVX512PROCEDURES_FULLVER = 
		1000U*gGMS_AVX512COMPLEX_MAJOR+100U*gGMS_AVX512COMPLEX_MINOR+10U*gGMS_AVX512COMPLEX_MICRO;

	const char * const pgGMS_AVX512PROCEDURES_CREATE_DATE = "12-10-2019 18:12 +00200 (SAT  12 OCT  2019 GMT+2)";

	const char * const pgGMS_AVX512PROCEDURES_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_AVX512PROCEDURES_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVX512PROCEDURES_SYNOPSIS = "Global (namespace) operator-like procedures.";
}


namespace gms {
	namespace math {

class AVX512VComplex1D;
			
		

		
		void
		avx512vcomplex_add(AVX512VComplex1D    &,
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool ) noexcept(true);

		
		void
		avx512vcomplex_add(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

		
		void
		avx512vcomplex_sub(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool) noexcept(true);

	
		void
		avx512vcomplex_sub(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

	
		void
		avx512vcomplex_mul(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool) noexcept(true);

	
		void
		avx512vcomplex_mul(AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

		void
		avx512vcomplex_div(AVX512VComplex1D	 &, 
				   const AVX512VComplex1D &,
				   const AVX512VComplex1D &,
				   const bool);
						

		
		void
		avx512vcomplex_div(AVX512VComplex1D    &, 
				   const AVX512VComplex1D &,
				   const double * __restrict,
				   const bool) noexcept(true);

	  //
	  //	*********** Procedures which compute equality/inequality of complex series implemented in LAM_avxvcprocedures.cpp ************
	  //
	}
}



#endif /*__GMS_AVX512VCPROCEDURES_H__*/
