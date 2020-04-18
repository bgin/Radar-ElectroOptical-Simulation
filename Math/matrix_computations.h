
#ifndef __MATRIX_COMPUTATIONS_H__
#define __MATRIX_COMPUTATIONS_H__


namespace file_info {

   const unsigned int gGMS_MATRIX_COMPUTATIONS_MAJOR = 1U;
   const unsigned int gGMS_MATRIX_COMPUTATIONS_MINOR = 0U;
   const unsigned int gGMS_MATRIX_COMPUTATIONS_MICRO = 0U;
   const unsigned int gGMS_MATRIX_COMPUTATIONS_FULLVER =
     1000U*gGMS_MATRIX_COMPUTATIONS_MAJOR+100U*gGMS_MATRIX_COMPUTATIONS_MINOR+
     10U*gGMS_MATRIX_COMPUTATIONS_MICRO;
   const char * const pgGMS_MATRIX_COMPUTATIONS_CREATE_DATE = "17-04-2020 11:10 +00200 (FRI 17 APR 2020 GMT+2)";
   const char * const pgGMS_MATRIX_COMPUTATIONS_BUILD_DATE  = __DATE__ ":" __TIME__;
   const char * const pgGMS_MATRIX_COMPUTATIONS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const pgGMS_MATRIX_COMPUTATIONS_SYNOPSIS    = "Matrix real and complex helper functions";
}

#include <cstdint>
#include <complex>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_indices.h"

namespace gms {

       namespace math {

               /*
                    Matrix exponential computation --  unrolled version 1
                */
		__ATTR_HOT__
		__ATTR_ALIGN__(16)
		static inline 
		void exp4x4m_cmplxr4v1( const std::complex<float> * __restrict __ATTR_ALIGN__(64) L,
		                      const std::complex<float> * __restrict __ATTR_ALIGN__(64) Q,
				      const std::complex<float> * __restrict __ATTR_ALIGN__(64) INVQ,
				      const float z,
				      std::complex<float> * __restrict __ATTR_ALIGN__(64) result) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        L      = (const std::complex<float>*)__builtin_assume_aligned(L,64);
			Q      = (const std::complex<float>*)__builtin_assume_aligned(Q,64);
			INVQ   = (const std::complex<float>*)__builtin_assume_aligned(INVQ,64);
			result = (std::complex<float>*)__builtin_assume_aligned(result,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(L,64);
			__assume_aligned(Q,64);
			__assume_aligned(INVQ,64);
			__assume_aligned(result,64);
#endif
                        __ATTR_ALIGNED__(64) std::complex<float> t0[16] = {};
			t0[0]  = std::exp(L[0]*z);
			t0[5]  = std::exp(L[1]*z);
			t0[10] = std::exp(L[2]*z);
			t0[15] = std::exp(L[3]*z);
	        }

		/*
                        Complex 4x4 matrix multiplication (single precision)
                 */
		 __ATTR_HOT__
		 __ATTR_ALIGN__(16)
		 static inline
		 void mul4x4m_cmplxr4(const std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
		                      const std::complex<float> * __restrict __ATTR_ALIGN__(64) b,
				      std::complex<float> * __restrict __ATTR_ALIGN__(64) c) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                      a   = (const std::complex<float>*)__builtin_assume_aligned(a,64);
		      b   = (const std::complex<float>*)__builtin_assume_aligned(b,64);
		      c   = (std::complex<float>*)__builtin_assume_aligned(c,64);
#elif defined __ICC || defined __INTEL_COMPILER
                      __assume_aligned(a,64);
		      __assume_aligned(b,64);
		      __assume_aligned(c,64);
#endif
                      c = {};
                      for(int32_t i = 0; i != 4; i += 2) {
                          for(int32_t j = 0; j != 4; j += 2) {
			      std::complex<float> s0 = {0.0f,0.0f};
			      std::complex<float> s1 = {0.0f,0.0f};
			      std::complex<float> s2 = {0.0f,0.0f};
			      std::complex<float> s3 = {0.0f,0.0f};
			      for(int32_t k = 0; k != 4; ++k) {
                                  std::complex<float> r1 = a[Ix2D(j,4,k)];
				  std::complex<float> r2 = b[Ix2D(k,4,i)];
				  s0 = s0 + r1 * r2;
				  std::complex<float> r3 = a[Ix2D(j+1,4,k)];
				  s1 = s1 + r3 * r2;
				  std::complex<float> r4 = b[Ix2D(k+1,4,i)];
				  s2 = s2 + r1 * r4;
				  s3 = s3 + r3 * r4;
			      }
			      c[Ix2D(i,4,j)]     = s0;
			      c[Ix2D(i+1,4,j)]   = s1;
			      c[Ix2D(i,4,j+1)]   = s2;
			      c[Ix2D(i+1,4,j+1)] = s3;
			  }
		      }
                 }
             /*
                  Helper function to multiply 3 complex matrices 4x4
              */
	      __ATTR_HOT__
	      __ATTR_ALIGN__(16)
	      static inline
	      void mul4x4m_cmplxr4_helper(const std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
	                                  const std::complex<float> * __restrict __ATTR_ALIGN__(64) b,
					  const std::complex<float> * __restrict __ATTR_ALIGN__(64) c,
					  std::complex<float> * __restrict __ATTR_ALIGN__(64) result) {
                  __ATTR_ALIGN__(64) std::complex<float> d[16];
		  mul4x4m_cmplxr4(a,b,d);
		  mul4x4m_cmplxr4(d,c,result);
	       }
             /*
                 The exponential of complex matrix 4x4 version 2
              */
	    __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	    static inline
	    void exp4x4m_cmplxr4v2(const std::complex<float> * __restrict __ATTR_ALIGN__(64) L,
	                           const std::complex<float> * __restrict __ATTR_ALIGN__(64) Q,
				   const std::complex<float> * __restrict __ATTR_ALIGN__(64) INVQ,
				   const float z,
				  float * __restrict __ATTR_ALIGN__(64) result) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        L      = (const std::complex<float>*)__builtin_assume_aligned(L,64);
			Q      = (const std::complex<float>*)__builtin_assume_aligned(Q,64);
			INVQ   = (const std::complex<float>*)__builtin_assume_aligned(INVQ,64);
			result = (float*)__builtin_assume_aligned(result,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(L,64);
			__assume_aligned(Q,64);
			__assume_aligned(INVQ,64);
			__assume_aligned(result,64);
#endif
#if defined __AVX512F__
			__ATTR_ALIGN__(64) std::complex<float> Diag[16] = {};
                        __ATTR_ALIGN__(64) float Qre[16]    = {};
			__ATTR_ALIGN__(64) float Qim[16]    = {};
			__ATTR_ALIGN__(64) float INVQre[16] = {};
			__ATTR_ALIGN__(64) float INVQim[16] = {};
			__ATTR_ALIGN__(64) float Lre[16]    = {};
			__ATTR_ALIGN__(64) float Lim[16]    = {};
			__ATTR_ALIGN__(64) float FRre[16]   = {};
			__ATTR_ALIGN__(64) float FRim[16]   = {};
			__ATTR_ALIGN__(64) float FIre[16]   = {};
			__ATTR_ALIGN__(64) float FIim[16]   = {};
			__ATTR_ALIGN__(64) float Diff[16]   = {};
			__ATTR_ALIGN__(16) float Sum[16]    = {};
#else
			__ATTR_ALIGN__(32) std::complex<float> Diag[16] = {};
                        __ATTR_ALIGN__(32) float Qre[16]    = {};
			__ATTR_ALIGN__(32) float Qim[16]    = {};
			__ATTR_ALIGN__(32) float INVQre[16] = {};
			__ATTR_ALIGN__(32) float INVQim[16] = {};
			__ATTR_ALIGN__(32) float Lre[16]    = {};
			__ATTR_ALIGN__(32) float Lim[16]    = {};
			__ATTR_ALIGN__(32) float FRre[16]   = {};
			__ATTR_ALIGN__(32) float FRim[16]   = {};
			__ATTR_ALIGN__(32) float FIre[16]   = {};
			__ATTR_ALIGN__(32) float FIim[16]   = {};
			__ATTR_ALIGN__(32) float Diff[16]   = {};
			__ATTR_ALIGN__(32) float Sum[16]    = {};
#endif
			for(int32_t i = 0; i != 16; ++i) {
                            Qre[i]    = Q[i].real();
			    Qim[i]    = Q[i].imag();
			    INVQre[i] = INVQ[i].real();
			    INVQim[i] = INVQ[i].imag();
			}
                        Diag[0] = std::exp(L[0]*z);
			Lre[0]  = Diag[0].real();
			Lim[0]  = Diag[0].imag();
			Diag[5] = std::exp(L[1]*z);
			Lre[5]  = Diag[5].real();
			Lim[5]  = Diag[5].imag();
			Diag[10] = std::exp(L[10]*z);
			Lre[10]  = Diag[10].real();
			Lim[10]  = Diag[10].imag();
			Diag[15] = std::exp(L[15]*z);
			Lre[15]  = Diag[15].real();
			Lim[15]  = Diag[15].imag();
#if defined __AVX512F__
                      
			__m512 zmm0,zmm1;
		      
			zmm0 = _mm512_setr_ps(Lre[15],Lre[10],Lre[5],Lre[0],
			                      Lre[15],Lre[10],Lre[5],Lre[0],
					      Lre[15],Lre[10],Lre[5],Lre[0],
					      Lre[15],Lre[10],Lre[5],Lre[0]);
			zmm1 = _mm512_setr_ps(Lim[15],Lim[10],Lim[5],Lim[0],
			                      Lim[15],Lim[10],Lim[5],Lim[0],
					      Lim[15],Lim[10],Lim[5],Lim[0],
					      Lim[15],Lim[10],Lim[5],Lim[0]);
			const __m512 qre = _mm512_load_ps(&INVQre[0]);
			const __m512 qim = _mm512_load_ps(&INVQim[0]);
			_mm512_store_ps(&FRre[0],_mm512_mul_ps(qre,zmm0));
			_mm512_store_ps(&FRim[0],_mm512_mul_ps(qim,zmm1));
			_mm512_store_ps(&FIre[0],_mm512_mul_ps(qre,zmm0));
			_mm512_store_ps(&FIim[0],_mm512_mul_ps(qim,zmm1));
			_mm512_store_ps(&Diff[0],_mm512_sub_ps(_mm512_load_ps(&FRre[0]),
			                                       _mm512_load_ps(&FIim[0])));
			_mm512_store_ps(&Sum[0], _mm512_add_ps(_mm512_load_ps(&FRim[0]),
			                                       _mm512_load_ps(&FIre[0])));
#else
                       	__m256 ymm0,ymm1;
			ymm0 = _mm256_setr_ps(Lre[15],Lre[10],Lre[5],Lre[0],
					      Lre[15],Lre[10],Lre[5],Lre[0]);
			ymm1 = _mm256_setr_ps(Lim[15],Lim[10],Lim[5],Lim[0],
					      Lim[15],Lim[10],Lim[5],Lim[0]);
		   
			const __m256 qrel = _mm256_load_ps(&QINVre[0]);
			const __m256 qreh = _mm256_load_ps(&QINVre[8]);
			const __m256 qiml = _mm256_load_ps(&QINVim[0]);
			const __m256 qimh = _mm256_load_ps(&QINVim[8]);
			_mm256_store_ps(&FRre[0],_mm256_mul_ps(qrel,ymm0));
			_mm256_store_ps(&FRre[8],_mm256_mul_ps(qreh,ymm0));
			_mm256_store_ps(&FRim[0],_mm256_mul_ps(qiml,ymm1));
			_mm256_store_ps(&FRim[8],_mm256_mul_ps(qimh,ymm1));
			_mm256_store_ps(&FIre[0],_mm256_mul_ps(qrel,ymm0));
			_mm256_store_ps(&FIre[8],_mm256_mul_ps(qreh,ymm0));
			_mm256_store_ps(&FIim[0],_mm256_mul_ps(qiml,ymm1));
			_mm256_store_ps(&FIim[8],_mm256_mul_ps(qimh,ymm1));
			_mm256_store_ps(&Diff[0],_mm256_sub_ps(_mm256_load_ps(&FRre[0]),
							       _mm256_load_ps(&FIim[0])));
			_mm256_store_ps(&Sum[0], _mm256_add_ps(_mm256_load_ps(&FRim[0]),
							       _mm256_load_ps(&FIre[0])));
			_mm256_store_ps(&Diff[8],_mm256_sub_ps(_mm256_load_ps(&FRre[8]),
							       _mm256_load_ps(&FIim[8])));
			_mm256_store_ps(&Sum[8], _mm256_add_ps(_mm256_load_ps(&FRim[8]),
							       _mm256_load_ps(&FIre[8])));
#endif
			// This should be vectorized.
			result = {};
		        for(int32_t i = 0; i != 4; i += 2) {
                           for(int32_t j = 0; j != 4; j += 2) {
                              float s0 = 0.0f;
			      float s1 = 0.0f;
			      float s2 = 0.0f;
			      float s3 = 0.0f;
			      for(int32_t k = 0; k != 4; ++k) {
                                  float r1 = Qre[Ix2D(j,4,k)];
				  float r2 = Diff[Ix2D(k,4,i)];
				  float r3 = Qim[Ix2D(j,4,k)];
				  float r4 = Sum[Ix2D(k,4,i)];
				  s0 = s0 + r1*r2-r3*r4;
				  float r5 = Qre[Ix2D(j+1,4,k)];
				  float r6 = Diff[Ix2D(k,4,i+1)];
				  float r7 = Qim[Ix2D(j+1,4,k)];
				  s1 = s1 + r5*r2-r7*r4;
				  float r8 = Sum[Ix2D(k,4,i+1)];
				  s2 = s2 + r1*r6-r3*r8;
				  s3 = s3 + r5*r6-r7*r8;
			      }
			      result[Ix2D(i,4,j)]     = s0;
			      result[Ix2D(i+1,4,j)]   = s1;
			      result[Ix2D(i,4,j+1)]   = s2;
			      result[Ix2D(i+1,4,j+1)] = s3; 
			   }
			}
		
		     }

		   /*
                        4x4 real matrix extinction
                    */
		   __ATTR_HOT__
		   __ATTR_ALIGN__(16)
		   static inline
		   void extinct_m4x4r4(const std::complex<float> * __restrict __ATTR_ALIGN__(32) M,
		                       float * __restrict __ATTR_ALIGN__(64) K) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                             M = (const std::complex<float>*)__builtin_assume_aligned(M,32);
			     K = (float*)__builtin_assume_aligned(K,64);
#elif defined __ICC || defined __INTEL_COMPILER
                             __assume_aligned(M,32);
			     __assume_aligned(K,64);
#endif
                             K[0] = -2.0f*M[0].real();
			     K[1] = 0.0f;
			     K[2] = -2.0f*M[1].real();
			     K[3] = -2.0f*M[1].imag();
			     K[4] = 0.0f;
			     K[5] = -2.0f*M[3].real();
			     K[6] = -2.0f*M[2].real();
			     K[7] = -2.0f*M[2].imag();
			     K[8] = -M[2].real();
			     K[9] = -M[1].real();
			     K[10] = -(M[0].real()+M[3].real());
			     K[11] = -(M[0].imag()-M[1].imag());
			     K[12] = -2.0f*M[2].real();
			     K[13] = -2.0f*M[1].real();
			     K[14] = -(M[0].imag()-M[1].imag());
			     K[15] = -(M[0].real()+M[3].real());
		    }
	           

   } // math

} // gms








#endif /*__MATRIX_COMPUTATIONS_H__*/
