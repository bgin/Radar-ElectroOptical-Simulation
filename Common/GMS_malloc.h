
#ifndef __GMS_MALLOC_H__
#define __GMS_MALLOC_H__

#if !defined (GMS_MALLOC_MAJOR)
#define GMS_MALLOC_MAJOR 1
#endif

#if !defined (GMS_MALLOC_MINOR)
#define GMS_MALLOC_MINOR 0
#endif

#if !defined (GMS_MALLOC_MICRO)
#define GMS_MALLOC_MICRO 0
#endif

#if !defined (GMS_MALLOC_FULLVER)
#define GMS_MALLOC_FULLVER 1000
#endif

#if !defined (GMS_MALLOC_CREATE_DATE)
#define GMS_MALLOC_CREATE_DATE "01-10-2019 19:14 +00200 (TUE 01 OCT 2019 GMT+2)"
#endif

#if !defined (GMS_MALLOC_BUILD_DATE)
#define GMS_MALLOC_BUILD_DATE " "
#endif

#if !defined (GMS_MALLOC_AUTHOR)
#define GMS_MALLOC_AUTHOR "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com"
#endif

#if !defined (GMS_MALLOC_DESCRIPT)
#define GMS_MALLOC_DESCRIPT "Malloc wrappers."
#endif

#include <cstdint>

class AVXVec8;
class AVX512Vec16;

namespace gms {
	namespace common {

#if defined _WIN64

		//
		//	Unaligned malloc wrapper
		//  Returns: double * 
		//  No error handling implemented
		//

		_Check_return_ _Ret_maybenull_	   double * gms_dmallocu(_In_ const std::size_t);

		//
		// Unaligned malloc wrapper
		// Returns: float * 
		// No error handling implemented
		//

		_Check_return_ _Ret_maybenull_      float * gms_fmallocu(_In_ const std::size_t);

		//
		// Unaligned malloc wrapper
		// Returns: int32_t * 
		// No error handling implemented
		//

		_Check_return_ _Ret_maybenull_      int32_t * gms_imallocu(_In_ const std::size_t);
		

		
		//
		//	Aligned malloc wrapper
		//  Returns: double * 
		//  No error handling implemented
		//							    _In_ const std::size_t);

		_Check_return_ _Ret_maybenull_      double * gms_dmalloca(_In_ const std::size_t, _In_ const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: float *
		// No error handling implemented
		//

		_Check_return_ _Ret_maybenull_      float * gms_fmalloca(_In_ const std::size_t, _In_ const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: int32_t * 
		// No error handling implemented
		//

	
		

		_Check_return_ _Ret_maybenull_      int32_t * gms_imalloca(_In_ const std::size_t, _In_ const int32_t);

		//
		// Error handling wrappers
		//

		//
		//	Aligned malloc wrapper
		//  Returns: double *
		//  Error checking and handling (calls std::exit)
		//

		_Check_return_ _Ret_maybenull_      double * gms_edmalloca(_In_ const std::size_t, _In_ int32_t);

		//
		// Aligned malloc wrapper
		// Returns: float *
		// Error checking and handling (calls std::exit)
		//

		_Check_return_ _Ret_maybenull_      float * gms_efmalloca(_In_ const std::size_t, _In_ const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: int32_t *
		// Error checking and handling (calls std::exit)
		//

		_Check_return_ _Ret_maybenull_     int32_t * gms_eimalloca4(_In_ const std::size_t, _In_ const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: int32_t *
		// Error checking and handling (calls std::exit)
		//

		_Check_return_ _Ret_maybenull_       int64_t * gms_eimalloca(_In_ const std::size_t, _In_ const int32_t);


		//
		// Few special functions for padding of possibly unaligned rows of flat (multidimensional) arrays.
		//

		_Check_return_ _Ret_maybenull_       
		 float * gms_efmalloca_padded2D(_In_ const int32_t,
									   _In_ const int32_t,
								       _In_ const int32_t,
									   _Inout_ int32_t &);

		_Check_return_ _Ret_maybenull_
		double * gms_edmalloca_padded2D(_In_ const int32_t,
										_In_ const int32_t,
										_In_ const int32_t,
										_Inout_ int32_t &);

		_Check_return_ _Ret_maybenull_
		int32_t * gms_eimalloca4_padded2D(_In_ const int32_t,
										 _In_ const int32_t,
										  _In_ const int32_t,
										 _Inout_ int32_t &);
#elif defined __linux

			//
		//	Unaligned malloc wrapper
		//  Returns: double * 
		//  No error handling implemented
		//

	        double * gms_dmallocu( const std::size_t);

		//
		// Unaligned malloc wrapper
		// Returns: float * 
		// No error handling implemented
		//

	        float * gms_fmallocu( const std::size_t);

		//
		// Unaligned malloc wrapper
		// Returns: int64_t * 
		// No error handling implemented
		//

	        int32_t * gms_imallocu( const std::size_t);

		

		//
		//	Aligned malloc wrapper
		//  Returns: double * 
		//  No error handling implemented
		//

		double * gms_dmalloca( const std::size_t, const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: float *
		// No error handling implemented
		//

	        float * gms_fmalloca( const std::size_t, const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: int32_t * 
		// No error handling implemented
		//

	        int32_t * gms_imalloca( const std::size_t, const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: AVXVec8 * 
		// No error handling implemented
		//
	        
		AVXVec8 * gms_avxvec8_malloca(const std::size_t,const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: AVX512Vec16 * 
		// No error handling implemented
		//

		AVX512Vec16 * gms_avx512vec16_malloca(const std::size_t,const int32_t);

		//
		// Error handling wrappers
		//

		//
		//	Aligned malloc wrapper
		//  Returns: double *
		//  Error checking and handling (calls std::exit)
		//

	        double * gms_edmalloca(const std::size_t, int32_t);

		//
		// Aligned malloc wrapper
		// Returns: float *
		// Error checking and handling (calls std::exit)
		//

	        float * gms_efmalloca(const std::size_t, const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: int32_t *
		// Error checking and handling (calls std::exit)
		//

	        int32_t * gms_eimalloca4(const std::size_t, const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: int64_t *
		// Error checking and handling (calls std::exit)
		//

	        int64_t * gms_eimalloca(const std::size_t, const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: AVXVec8 * 
		// Error checking and handling (calls std::exit)
		//

		AVXVec8 * gms_avxvec8_emalloca(const std::size_t,
					       const int32_t);

		//
		// Aligned malloc wrapper
		// Returns: AVXVec8 * 
		// Error checking and handling (calls std::exit)
		//
		AVX512Vec16 * gms_avx512vec16_emalloca(const std::size_t,
						       const int32_t);

		//
		// Few special functions for padding of possibly unaligned rows of flat (multidimensional) arrays.
		//

	       
		float * gms_efmalloca_padded2D(const int32_t,
						const int32_t,
						const int32_t,
						int32_t &);

	
		double * gms_edmalloca_padded2D(const int32_t,
						const int32_t,
						const int32_t,
					        int32_t &);

		
		int32_t * gms_eimalloca4_padded2D(const int32_t,
						  const int32_t,
						  const int32_t,
						  int32_t &);

#endif

		

	}
}


#endif /*__GMS_MALLOC_H__*/
