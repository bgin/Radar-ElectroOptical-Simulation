
#ifndef __GMS_AVXVECF64_H__
#define __GMS_AVXVECF64_H__ 141020191209





namespace file_info {


	const unsigned int gGMS_AVXVECF64_MAJOR = 1;

	const unsigned int gGMS_AVXVECF64_MINOR = 0;

	const unsigned int gGMS_AVXVECF64_MICRO = 1;

	const unsigned int gGMS_AVXVECF64_FULLVER = 
	 1000U*gGMS_AVXVECF64_MAJOR+100U*gGMS_AVXVECF64_MINOR+10U*gGMS_AVXVECF64_MICRO;

	const char * const pgGMS_AVXVECF64_CREATE_DATE = "14-10-2019 12:09 +00200 (MON 14 OCT 2019 GMT+2)";

	const char * const pgGMS_AVXVECF64_BUILD_DATE =  __DATE__":"__TIME__;

	const char * const pgGMS_AVXVECF64_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVXVECF64_SYNOPSIS =  "Wrapper class around __m256d data type.";
}

#include <iostream>
#include <iomanip>
#include <cstdint>

#include "GMS_simd_defs.h"


namespace gms {
	namespace math {

		__attribute__((align(64))) struct AVXVecF64 {


		

				//
				// Class Constructors and Destructor.
				//
			 	AVXVecF64();

				AVXVecF64(const double[4]);

				AVXVecF64(const double, 
					  const double,
					  const double,
					  const double);

				AVXVecF64(const double);

				AVXVecF64(const __m256d);

				AVXVecF64(const __m256i);

				AVXVecF64(const AVXVecF64 &);

				AVXVecF64(const __m128d, 
					  const __m128d);

				~AVXVecF64() = default;

				//
				// getters.
				//
				const __m256d get_vf64() const;

			         __m256d get_vf64();
 
				__m128d lo_part() const;

			 	__m128d hi_part() const;

				//
				// Load-store functions
				//

				// Address argument should be aligned on 32-byte boundary
				AVXVecF64 & load_a(const double* __restrict);

				AVXVecF64 & load_u(const double* __restrict);
				// Address argument should be aligned on 32-byte boundary
				void  store_a(double* __restrict) const;

				void  store_u(double* __restrict) const;

				void  stream_store(double* __restrict) const;

				double extract_scalar(const uint32_t) const;

				// Inserts a single double in location 0...3
				// Based on vectorclass insert function (Agner Fog)
				AVXVecF64 const & insert(const uint32_t, 
				                         const double  );

				//
				// Few member and friend operators
				//
				AVXVecF64 & operator=(const AVXVecF64 &);

				

				//
				// Type cast operator
				//
				operator __m256d () const;

				double operator[](const unsigned int) const;

				friend std::ostream & operator<<(std::ostream &,
				                                 const AVXVecF64 &);


			

		                 __m256d m_vf64;
		};

		//
		// global(namespace) static functions
		//

		//
		// Extract __m128d part. Value of second parameter
		// must be 0 or 1 only.
		static inline __m128d extract(AVXVecF64 &,
					      const int32_t);

		

		// branchless conditional selection
		static inline AVXVecF64 select_vec(const AVXVecF64 &, 
						   const AVXVecF64 &,
						   const __m256d);
	    
		//
		//	Arithmetic and mathematical operations
		//

		// SIMD max
		static inline AVXVecF64 max(const AVXVecF64 &, 
					    const AVXVecF64 &);

		// SIMD min
		static inline AVXVecF64 min(const AVXVecF64 &,
					    const AVXVecF64 &);

		// SIMD abs
		static inline AVXVecF64 abs(const AVXVecF64 &);

		// SIMD SQRT
		static inline AVXVecF64 sqrt(const AVXVecF64 &);

		// SIMD squared
		static inline AVXVecF64 sqr(const AVXVecF64 &);

		// SIMD ceil
		static inline AVXVecF64 ceil(const AVXVecF64 &);

		// SIMD floor
		static inline AVXVecF64 floor(const AVXVecF64 &);

		// SIMD round
		// Caller must pass either 0 or 1.
		static inline AVXVecF64 round(const AVXVecF64 &, 
					      const int32_t);
		
		// SVML sin
		static inline AVXVecF64 sin(const AVXVecF64 &);

		// SVML cos
		static inline AVXVecF64 cos(const AVXVecF64 &);

		// SVML sinh
		static inline AVXVecF64 sinh(const AVXVecF64 &);

		// SVML cosh
		static inline AVXVecF64 cosh(const AVXVecF64 &);

		// SVML tan
		static inline AVXVecF64 tan(const AVXVecF64 &);

		// SVML tanh
		static inline AVXVecF64 tanh(const AVXVecF64 &);

		// SVML asin
		static inline AVXVecF64 asin(const AVXVecF64 &);

		// SVML asinh
		static inline AVXVecF64 asinh(const AVXVecF64 &);

		// SVML acos
		static inline AVXVecF64 acos(const AVXVecF64 &);

		// SVML acosh
		static inline AVXVecF64 acosh(const AVXVecF64 &);

		// SVML atan
		static inline AVXVecF64 atan(const AVXVecF64 &);

		// SVML atanh
		static inline AVXVecF64 atanh(const AVXVecF64 &);

		// AVML atan2
		static inline AVXVecF64 atan2(const AVXVecF64 &,const AVXVecF64 &);

		// Unary minus
		static inline AVXVecF64 unary_minus(const AVXVecF64 &);

		// SVML exp
		static inline AVXVecF64 exp(AVXVecF64 &);

		// SVML log10
		static inline AVXVecF64 log10(const AVXVecF64 &);

		// SVML log
		static inline AVXVecF64 log(const AVXVecF64 &);

		// SVML pow
		static inline AVXVecF64 pow(const AVXVecF64 &,
					    const AVXVecF64 &);

		// FMA functions

		// fmadd
		static inline AVXVecF64 fmadd(const AVXVecF64 &, 
					      const AVXVecF64 &,
					      const AVXVecF64 &);

		// fmadsubb
		static inline AVXVecF64 fmadsubb(const AVXVecF64 &,
						 const AVXVecF64 &,
						 const AVXVecF64 &);

		// fmsub
		static inline AVXVecF64 fmsub(const AVXVecF64 &,
					      const AVXVecF64 &,
					      const AVXVecF64 &);

		// fmsubadd
		static inline AVXVecF64 fmsubadd(const AVXVecF64 &,
						 const AVXVecF64 &,
						 const AVXVecF64 &);

		// fnmadd
		static inline AVXVecF64 fnmadd(const AVXVecF64 &,
					       const AVXVecF64 &,
					       const AVXVecF64 &);

	   // fnmsub
		static inline AVXVecF64 fnmsub(const AVXVecF64 &,
					       const AVXVecF64 &,
					       const AVXVecF64 &);

		//
		// static operators
		//

		// C = A+B, vector + vector
		static inline AVXVecF64 operator+(const AVXVecF64 &, 
						  const AVXVecF64 &);

		// C = A+B, vector + scalar
		static inline AVXVecF64 operator+(const AVXVecF64 &, 
						  const double);

		// C = A+B, scalar + vector
		static inline AVXVecF64 operator+(const double, 
						  const AVXVecF64 &);

		// A = A+B, vector + vector (in-place)
		static inline AVXVecF64 operator+=(AVXVecF64 &,
						   const AVXVecF64 &);

		// A = A+1
		static inline AVXVecF64 operator++(AVXVecF64 &);

	    // C = A-B, vector - vector
		static inline AVXVecF64 operator-(const AVXVecF64 &,
						  const AVXVecF64 &);

		// C = A-B, vector - scalar
		static inline AVXVecF64 operator-(const AVXVecF64 &, 
						  const double);

		// C = A-B, scalar - vector
		static inline AVXVecF64 operator-(const double, 
						  const AVXVecF64 &);

		// A = A-B, vector - vector (in-place)
		static inline AVXVecF64 operator-=(AVXVecF64 &, 
						   const AVXVecF64 &);

		// A = A-1.0L
		static inline AVXVecF64 operator--(AVXVecF64 &);

		// C = A*B, vector * vector
		static inline AVXVecF64 operator*(const AVXVecF64 &, 
						  const AVXVecF64 &);

	    // C = A*B, vector * scalar
		static inline AVXVecF64 operator*(const AVXVecF64 &,
						  const double);

		// C = A*B, scalar * vector
		static inline AVXVecF64 operator*(const double,
						  const AVXVecF64 &);

		// A = A*B, vector * vector (in-place)
		static inline AVXVecF64 operator*=(AVXVecF64 &, 
						   const AVXVecF64 &);

		// C = A/B, vector / vector
		static inline AVXVecF64 operator/(const AVXVecF64 &, 
						  const AVXVecF64 &);

		// C = A/B, vector / scalar
		static inline AVXVecF64 operator/(const AVXVecF64 &,
						  const double);

		// C = A/B, scalar / vector
		static inline AVXVecF64 operator/(_In_ const double, 
						  const AVXVecF64 &);

		// A = A/B, vector / vector (in-place)
		static inline AVXVecF64 operator/=(AVXVecF64 &, 
						   const AVXVecF64 &);

		// C = A==B, vector == vector, C is of type __m256d
		static inline __m256d operator==(const AVXVecF64 &, 
						 const AVXVecF64 &);

		static inline __m256d operator==(const AVXVecF64 &,
						 const double );

		static inline __m256d operator==(const double,
						 const AVXVecF64 &);

		// C = A != b, vector != vector, C is of type __m256d
		static inline __m256d operator!=(const AVXVecF64 &, 
						 const AVXVecF64 &);

		// C = A>B, vector > vector, C is of type __m256d
		static inline __m256d operator>(const AVXVecF64 &, 
						const AVXVecF64 &);

		// C = A<B, vector < vector, C is of type __m256d
		static inline __m256d operator<(const AVXVecF64 &, 
						const AVXVecF64 &);

		// C = A>=B, vector >= B, C is of type __m256d
		static inline __m256d operator>=(const AVXVecF64 &, 
						 const AVXVecF64 &);

		// C = A<=B, vector <= vector, C is of type __m256d
		static inline __m256d operator<=(const AVXVecF64 &, 
						 const AVXVecF64 &);

		// C = A&B, vector & vector
		static inline AVXVecF64 operator&(const AVXVecF64 &,
						  const AVXVecF64 &);

		
		// A = A&B, vector & vector (in-place)
		static inline AVXVecF64 operator&=(AVXVecF64 &,
						   const AVXVecF64 &);


		// C = A | B, vector | vector
		static inline AVXVecF64 operator|(const AVXVecF64 &,
						  const AVXVecF64 &);

		
		// A = A | B, vector | vector
		static inline AVXVecF64 operator|=(AVXVecF64 &,
						   const AVXVecF64 &);

		// C = A ^ B, vector ^ vector
		static inline AVXVecF64 operator^(const AVXVecF64 &,
						  const AVXVecF64 &);

		// A = A ^ B, vector ^ vector
		static inline AVXVecF64 operator^=(AVXVecF64 &,
						   const AVXVecF64 &);


#include "GMS_avxvecf64.inl"
	}
}

#endif /*__GMS_AVXVECF64_H__*/
