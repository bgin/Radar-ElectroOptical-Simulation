
#ifndef __GMS_AVXVECF32_H__
#define __GMS_AVXVECF32_H__


namespace file_info{


       const unsigned int gGMS_AVXVECF32_MAJOR = 1;
       const unsigned int gGMS_AVXVECF32_MINOR = 0;
       const unsigned int gGMS_AVXVECF32_MICRO = 0;
       const unsigned int gGMS_AVXVECF32_FULLVER =
          1000U*gGMS_AVXVECF32_MAJOR+100U*gGMS_AVXVECF32_MINOR+10U*gGMS_AVXVECF32_MICRO;
       const char * const pgGMS_AVXVECF32_CREATION_DATE = "22-12-2019 12:09 +00200 (SUN 22 DEC 2019 GMT+2)";
       const char * const pgGMS_AVXVECF32_BUILD_DATE    = __DATE__ " " __TIME__;
       const char * const pgGMS_AVXVECF32_AUTHOR        = "Programmer: Bernard Gingold, contact; beniekg@gmail.com";
       const char * const pgGMS_AVXVECF32_SYNOPSIS      = "Wrapper class around __m256 data type.";
}

#include <immintrin.h>
#include <cstdint>
#include "Gms_config.h"


namespace gms {
       namespace math {

                 struct AVXVec8 {


		     __m256 m_v8 __ATTR_ALIGN__(32) ;

                     AVXVec8();
		     
		     AVXVec8(const float* __restrict __ATTR_ALIGN__(32));

		     AVXVec8(const float,
		               const float,
			       const float,
			       const float,
			       const float,
			       const float,
			       const float,
			       const float);

		     AVXVec8(const float);

		     AVXVec8(const __m256);

		     AVXVec8(const __m256i);

		     AVXVec8(const AVXVecF32 &);

		     AVXVec8(const __m128,
		             const __m128);

		    ~AVXVec8() = default;

		    __m128 lo_part() const;

		    __m128 hi_part() const;

		    //
		    // Load-store functions
		    //
		    // Address aligned on 32-byte boundary
                    AVXVec8   & load_a(const float * __restrict __ATTR_ALIGN__(32));

		    AVXVec8   & load_u(const float * __restrict __ATTR_ALIGN__(32));
		    
		    // Address argument should be aligned on 32-byte boundary
		    void store_a(float * __restrict __ATTR_ALIGN__(32)) const;

		    void store_u(float * __restrict __ATTR_ALIGN__(32)) const;

		    void stream_store(float * __restrict __ATTR_ALIGN__(32)) const;

		    float extract_scalar(const uint32_t) const;
		    
                    // Inserts a single double in location 0...7
		    AVXVec8 & insert_scalar(const uint32_t,
		                            const float);

		    AVXVec8 & operator=(const AVXVec8 &);

		    AVXVec8 & operator=(const AVXVec8);

		    operator __m256 () const;

		    float operator[](const uint32_t) const;

	   } __attribute__ ((aligned(64)));

	   	//
		// global(namespace) static functions
		//

		//
		// Extract __m128 part. Value of second parameter
		// must be 0 or 1 only.
                static inline __m128 extract_vec4(AVXVec8,
		                                  const uint32_t)      __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 select_vec8(const AVXVec8,
		                                  const AVXVec8,
						  const __m256)        __ATTR_HOT__  __ATTR_ALIGN__(16);

		//
		//	Arithmetic and mathematical operations
		//
                static inline AVXVec8 max(const AVXVec8,
		                          const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 min(const AVXVec8,
		                          const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 abs(const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);
		
		static inline AVXVec8 sqrt(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 ceil(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 floor(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 round(const AVXVec8,
		                            int)                       __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 addsub(const AVXVec8,
		                             const AVXVec8)            __ATTR_HOT__  __ATTR_ALIGN__(16);
					     
                static inline AVXVec8 dot(const AVXVec8,
		                          const AVXVec8,
					  const int32_t)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 hadd(const AVXVec8,
		                           const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 hsub(const AVXVec8,
		                           const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline float extract_1f(const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 rsqrt(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);
		
	        static inline AVXVec8 rcp(const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 i32gather(float const* __restrict,
		                                __m256i, const int)    __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 mask_i32gather(AVXVec8,
		                                     float const* __restrict,
						     __m256i, __m256,
						     const int)        __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 mask_load(float const* __restrict,
		                                __m256i)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		// Return true if all bits are 1
	        static inline bool    testc(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		// Return true if at least one bit is set (1)
		static inline bool    testz(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 permute(const AVXVec8,
		                              int)                     __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 unary_minus(const AVXVec8)       __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 sin(const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 cos(const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 sinh(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 cosh(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 tan(const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 tanh(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 asin(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 asinh(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 acos(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 acosh(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 atan(const AVXVec8)              __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 atanh(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 atan2(const AVXVec8,
		                            const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

	        static inline AVXVec8 exp(const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 log10(const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 log(const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 pow(const AVXVec8,
		                          const AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

	        static inline AVXVec8 fmadd(const AVXVec8,
		                            const AVXVec8,
					    const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

	        static inline AVXVec8 fmadsubb(const AVXVec8,
		                               const AVXVec8,
					       const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 fmsub(const AVXVec8,
		                            const AVXVec8,
					    const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 fmsubadd(const AVXVec8,
		                               const AVXVec8,
					       const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 fnmadd(const AVXVec8,
		                             const AVXVec8,
					     const AVXVec8)            __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline AVXVec8 fnmsub(const AVXVec8,
		                             const AVXVec8,
					     const AVXVec8)            __ATTR_HOT__  __ATTR_ALIGN__(16);

		//
		// static operators
		//

		// C = A+B, vector + vector
		static inline AVXVec8 operator+(const AVXVec8, 
						const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16); 

		// C = A+B, vector + scalar
		static inline AVXVec8 operator+(const AVXVec8, 
						  const float)         __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A+B, scalar + vector
		static inline AVXVec8 operator+(const float, 
						  const AVXVec8)        __ATTR_HOT__  __ATTR_ALIGN__(16);

		// A = A+B, vector + vector (in-place)
		static inline AVXVec8 operator+=(AVXVec8,
						 const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16);

		// A = A+1
		static inline AVXVec8 operator++(AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

	    // C = A-B, vector - vector
		static inline AVXVec8 operator-(const AVXVec8,
						const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16); 

		// C = A-B, vector - scalar
		static inline AVXVec8 operator-(const AVXVec8, 
						  const float)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A-B, scalar - vector
		static inline AVXVec8 operator-(const float, 
						const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16);  

		// A = A-B, vector - vector (in-place)
		static inline AVXVec8 operator-=(AVXVec8, 
						 const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16);

			// A = A-1.0L
		static inline AVXVec8 operator--(AVXVec8)               __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A*B, vector * vector
		static inline AVXVec8 operator*(const AVXVec8, 
						  const AVXVec8)        __ATTR_HOT__  __ATTR_ALIGN__(16);

	    // C = A*B, vector * scalar
		static inline AVXVec8 operator*(const AVXVec8,
						  const float)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A*B, scalar * vector
		static inline AVXVec8 operator*(const float,
						const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		// A = A*B, vector * vector (in-place)
		static inline AVXVec8 operator*=(AVXVec8, 
						 const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A/B, vector / vector
		static inline AVXVec8 operator/(const AVXVec8, 
						  const AVXVec8)        __ATTR_HOT__  __ATTR_ALIGN__(16); 

		// C = A/B, vector / scalar
		static inline AVXVec8 operator/(const AVXVec8
						  const float)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A/B, scalar / vector
		static inline AVXVec8 operator/(const float, 
						  const AVXVec8)        __ATTR_HOT__  __ATTR_ALIGN__(16);

		// A = A/B, vector / vector (in-place)
		static inline AVXVec8 operator/=(AVXVec8, 
						 const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16);

	        //

		// C = A==B, vector == vector, C is of type __m256d
		static inline __m256 operator==(const AVXVec8, 
						 const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16); 

		static inline __m256 operator==(const AVXVec8
						 const float )           __ATTR_HOT__  __ATTR_ALIGN__(16);

		static inline __m256 operator==(const float,
						 const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A != b, vector != vector, C is of type __m256d
		static inline __m256 operator!=(const AVXVec8, 
						 const AVXVec8)           __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A>B, vector > vector, C is of type __m256d
		static inline __m256 operator>(const AVXVec8, 
						const AVXVec8)             __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A<B, vector < vector, C is of type __m256d
		static inline __m256 operator<(const AVXVec8, 
						const AVXVec8)            __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A>=B, vector >= B, C is of type __m256d
		static inline __m256 operator>=(const AVXVec8, 
						 const AVXVec8)            __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A<=B, vector <= vector, C is of type __m256d
		static inline __m256 operator<=(const AVXVec8, 
						 const AVXVec8)            __ATTR_HOT__  __ATTR_ALIGN__(16);

			// C = A&B, vector & vector
		static inline AVXVec8 operator&(const AVXVec8,
						  const AVXVec8)           __ATTR_HOT__  __ATTR_ALIGN__(16);

		
		// A = A&B, vector & vector (in-place)
		static inline AVXVec8 operator&=(AVXVec8,
						   const AVXVec8)          __ATTR_HOT__  __ATTR_ALIGN__(16);


		// C = A | B, vector | vector
		static inline AVXVec8 operator|(const AVXVec8,
						  const AVXVec8)           __ATTR_HOT__  __ATTR_ALIGN__(16);

		
		// A = A | B, vector | vector
		static inline AVXVec8 operator|=(AVXVec8,
						   const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16);

		// C = A ^ B, vector ^ vector
		static inline AVXVec8 operator^(const AVXVec8,
						  const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16);

		// A = A ^ B, vector ^ vector
		static inline AVXVec8 operator^=(AVXVec8,
						   const AVXVec8)         __ATTR_HOT__  __ATTR_ALIGN__(16);				 




#include "GMS_avxvecf32.inl"
    }
}





#endif  /*__GMS_AVXVECF32_H__*/
