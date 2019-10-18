
#ifndef __GMS_AVXC2F64_H__
#define __GMS_AVXC2F64_H__



#include <iosfwd>
#include <iomanip>
#include <cstdint>
#if defined _WIN64
    #include "../GMS_simd_defs.h"
#elif defined __linux
    #include "GMS_simd_defs.h"
#endif

namespace file_info {

#if defined _WIN64  
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

     const unsigned int gGMS_AVXC2F64_MAJOR = gms::common::gVersionInfo.m_VersionMajor;
     const unsigned int gGMS_AVXC2F64_MINOR = gms::common::gVersionInfo.m_VersionMinor;
     const unsigned int gGMS_AVXC2F64_MICRO = gms::common::gVersionInfo.m_VersionMicro;
     const unsigned int gGMS_AVXC2F64_FULLVER =
       1000U*gGMS_AVXC2F64_MAJOR+100U*gGMS_AVXC2F64_MINOR+10U*gGMS_AVXC2F64_MICRO;
     const char * const pgGMS_AVXC2F64_CREATE_DATE =  "17-11-2017 11:13 +00200 (FRI 17 NOV 2017 GMT+2)";
     const char * const pgGMS_AVXC2F64_BUILD_DATE  =  "00-00-0000 00:00";
     const char * const pgGMS_AVXC2F64_AUTHOR      =  "Programmer: Bernard Gingold e-mail: beniekg@gmail.com";
     const char * const pgGMS_AVXC2F64_SYNOPSIS    =  "Wrapper struct around  __m256d which represents 2 complex numbers packed.";
}

namespace gms {
	namespace math {

		

	__declspec(align(64)) struct AVXC2f64 {


			//
			// Constructors and Destructor (default)
			//
		

			// Implemented default constructor sets member to 
			// default values (0.0,0.0) , (0.0,0.0)
			AVXC2f64();

			// Constructor from 4 doubles i.e. re,im,re,im
			AVXC2f64( const double,
				  const double,
				  const double,
				  const double);

			// Constructs from single complex number
			AVXC2f64(const double,
				 const double);

			// Construct from single re value
			explicit AVXC2f64(const double);

			// Constructs from only real parts (bool type is used for overload resolution only)
			AVXC2f64(const double,
				 const double,
				 const double,
				 const double,
				 const bool);

			// Copy-Construcor
			AVXC2f64(const AVXC2vf64 &);

			// Converting Constructor from __m256d type.
			AVXC2f64(const __m256d &);

			// Destructor default.
			~AVXC2f64() = default;

			//
			// getters and setters
			//
			const __m256d get_cv64() const;

			__m256d get_cv64();

			__m128d complex_1() const;

			__m128d complex_2() const;

			void set_complex_1(const double,
					   const double);

			void set_complex_2(const double,
					   const double);

			void set_complex_12(const double,
					    const double,
					    const double,
					    const double);

			//
			// Load-store function
			//
			// Address argument should be aligned on 32-byte boundary
			AVXC2f64 & load_a(const double * __restrict);

			// Unaligned memory load
			AVXC2f64 & load_u(const double * __restrict);

			// Address argument should be aligned on 32-byte boundary
			void store_a(double * __restrict ) const;

			// Unaligned memory store
			void store_u(double * __restrict) const;

			// Stream store (bypass caching stores)
			void stream_store(double * __restrict) const;

			// Extract single component of two packed complex numbers
	                double extract_component(const int32_t ) const;

			// Inserts a single double in location 0...3
			// Based on vectorclass insert function (Agner Fog)
			AVXC2f64 const & insert(const int32_t,
					         const double);

			//
			// Member operators
			//
			AVXC2f64 & operator=(const AVXC2f64 &);

			// Type cast operator
			operator __m256d () const;

			// Subscripting
			double operator[](const uint32_t) const;

			// ostream
			friend std::ostream & operator<<(std::ostream &,
						         const AVXC2f64 &);
		

			//
			// Packed 2-complex numbers i.e. re,im,re,im
			//
			 __m256d m_cv64;
		};

		// Complex sine
		static  AVXC2f64 csin(const AVXC2f64);

		static AVXC2f64 csin(const double,
				     const double);

		// Complex cosine
		static  AVXC2f64 ccos(const AVXC2f64);

		static  AVXC2f64 ccos(const double,
				      const double);

		// Complex exp
		static  AVXC2f64 cexp(const AVXC2f64);

		static  AVXC2f64 cexp(const double,
				      const double);

		// Complex ABS
		static inline std::pair<double,double> cabs(const AVXC2f64 );

		static inline double cabs(const double,
				          const double);

		// Complex to integer power (low part only)
		static inline AVXC2f64 cpowi(const AVXC2f64,
					     const int32_t );

		// Complex to integer power (low and high parts) i.e. 2 complexes
		static inline AVXC2f64 cpowi2(const AVXC2f64,
					      const int32_t);

		//
		// global(namespace) static functions
		//

		//
		// Extract __m128d part i.e. first or second complex component. Value of second parameter
		// must be 0 or 1 only.
		static inline __m128d extract(AVXC2f64,
					      const int32_t);

		// branchless conditional selection
		static inline AVXC2f64  select(const AVXC2f64 ,
					       const AVXC2f64 ,
					       const __m256d );



		//
		// Global operators
		//

		// C = A+B, complex + complex
		static inline AVXC2f64  operator+(const AVXC2f64,
					          const AVXC2vf64);

		// C = A+B, complex + scalar
		static inline AVXC2f64  operator+(const AVXC2f64,
					          const double);

		// C = A+B, scalar + complex
		static inline AVXC2f64  operator+(const double,
					           const AVXC2f64);

		// A = A+B, complex + complex (in-place)
		static inline AVXC2f64  operator+=(AVXC2f64,
						   const AVXC2f64);

		// A = A+B, complex + scalar (in-place)
		static inline AVXC2f64  operator+=(AVXC2f64,
					          const double);

		// A = A+B, scalar + complex (in-place)
		static inline AVXC2f64  operator+=(const double,
						   AVXC2f64);

		// C = A-B, complex - complex
		static inline AVXC2f64  operator-(const AVXC2f64,
					          const AVXC2f64);

		// C = A-B, complex - scalar
		static inline AVXC2f64  operator-(const AVXC2vf64,
						   const double);

		// C = A-B, scalar - complex
		static inline AVXC2f64  operator-(const double,
						  const AVXC2f64);

		// A = A-B, complex - complex (in-place)
		static inline AVXC2f64  operator-=(AVXC2f64,
						   const AVXC2f64);

		// A = A-B complex - scalar (in-place)
		static inline AVXC2f64  operator-=(AVXC2f64,
						   const double );

		// A = -B scalar - complex (in-place)
		static inline AVXC2f64  operator-=(const double,
						   AVXC2f64 );

	    // A = -A
		static inline AVXC2f64 operator-(AVXC2f64) ;

		// C = A*B, complex * complex
		static inline AVXC2f64  operator*(const AVXC2f64,
					         const AVXC2f64);

	    // C = A*B, complex * scalar
		static inline AVXC2f64  operator*(const AVXC2f64,
						  const double);

		// C = A*B, scalar * complex
		static inline AVXC2f64  operator*(const double,
					         const AVXC2f64);

		// A = A*B, complex * complex (in-place)
		static inline AVXC2f64  operator*=(AVXC2f64,
					          const AVXC2f64);

		// A = A*B complex * scalar (in-place)
		static inline AVXC2f64 operator*=(AVXC2f64,
					         const double);

		// A = A*B scalar * complex (in-place)
		static inline AVXC2f64 operator*=(const double,
					         AVXC2f64);

		// C = A/B, complex / complex
		static inline AVXC2f64  operator/(const AVXC2f64,
					          const AVXC2f64);

		// C = A/B, complex / scalar
		static inline AVXC2f64  operator/(const AVXC2f64,
					         const double);

		// C = A/B, scalar / complex
		static inline AVXC2f64  operator/(const double,
					         const AVXC2f64);

		// A = A/B, complex / complex (in-place)
		static inline AVXC2f64  operator/=(AVXC2f64,
					          const AVXC2f64);

		// A = A/B complex / scalar (in-place)
		static inline AVXC2f64  operator/=(AVXC2f64,
					           const double);

		// A = A/B scalar / complex (in-place)
		static inline AVXC2f64 operator/=(const double,
					          AVXC2f64);

		// Complex conjugate.
		static inline AVXC2f64  operator~(AVXC2f64);

		// Complex equality ==
		static inline AVXC2f64  operator==(const AVXC2f64,
						   const AVXC2f64) ;

		// Complex inequality
		static inline AVXC2f64  operator!=(const AVXC2f64,
						   const AVXC2f64);

		// Complex comparison >
		static inline AVXC2f64  operator>(const AVXC2f64,
					          const AVXC2f64);

		// Complex comparison >=
		static inline AVXC2f64  operator>=(const AVXC2f64,
					          const AVXC2f64);

		// Complex comparison <
		static inline AVXC2f64  operator<(const AVXC2f64,
						  const AVXC2f64);

		// Complex comparison <=
		static inline AVXC2f64  operator<=(const AVXC2f64,
						   const AVXC2f64);



#include "GMS_avxc2f64.inl"
	}
}


#endif /*__GMS_AVXC2F64_H__*/
