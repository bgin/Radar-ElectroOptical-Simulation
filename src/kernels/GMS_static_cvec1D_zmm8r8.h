
#ifndef __GMS_STATIC_CVEC1D_ZMM8R8_H__
#define __GMS_STATIC_CVEC1D_ZMM8R8_H__ 131120231523

namespace file_info {


	const unsigned int GMS_STATIC_CVEC1D_ZMM8R8_MAJOR = 1U;

	const unsigned int GMS_STATIC_CVEC1D_ZMM8R8_MINOR = 0U;

	const unsigned int GMS_STATIC_CVEC1D_ZMM8R8_MICRO = 0U;

	const unsigned int GMS_STATIC_CVEC1D_ZMM8R8_FULLVER = 
		1000U*GMS_STATIC_CVEC1D_ZMM8R8_MAJOR+100U*GMS_STATIC_CVEC1D_ZMM8R8_MINOR+10U*GMS_STATIC_CVEC1D_ZMM8R8_MICRO;

	const char * const GMS_STATIC_CVEC1D_ZMM8R8_CREATE_DATE = "13-11-2023 15:23 +00200 (MON 13 NOV 2023 GMT+2)";

	const char * const GMS_STATIC_CVEC1D_ZMM8R8_BUILD_DATE = __DATE__ ":" __TIME__;

	const char * const GMS_STATIC_CVEC1D_ZMM8R8_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const GMS_STATIC_CVEC1D_ZMM8R8_SYNOPSIS = "AVX512 complex vector (1D) stack-allocated storage.";
}

//
// Warning:
//				Include these files if and only if you have 
//				CPU and/or Accelarator i.e Xeon Phi which supports AVX512 ISA,
//				otherwise remove these files from compilation.
//

#include <cstdint>
#include <array>
#include <iostream>
#include <immintrin>
#include "GMS_config.h"
#include "GMS_common.h"
#include "GMS_complex_common_zmm8r8.h"



#if !defined (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) // Streaming stores defined per this struct (default set to 0)
#define USE_STATIC_CVEC1D_ZMM8R8_NT_STORES 0
#endif

namespace gms {
	namespace math {

#if !defined (STATIC_CVEC1D_ZMM8R8_LOAD_ZMM)
#define STATIC_CVEC1D_ZMM8R8_LOAD_ZMM(reg1,reg2,reg3,reg4,v1,v2,idx,off) \
	(reg1) = _mm512_load_pd(&(v1).m_Re[(idx)+(off)]);			   \
	(reg2) = _mm512_load_pd(&(v2).m_Re[(idx)+(off)]);			   \
	(reg3) = _mm512_load_pd(&(v1).m_Im[(idx)+(off)]);			   \
	(reg4) = _mm512_load_pd(&(v2).m_Im[(idx)+(off)]);
#endif



	


		template<int32_t N>
		__ATTR_ALIGN__(64) struct SCVec1DZMM8r8 {

		       __ATTR_ALIGN__(64) double m_Re[(m_nsize == 0) ? 8 : N];
		       __ATTR_ALIGN__(64) double m_Im[(m_nsize == 0) ? 8 : N];
                        int32_t m_nsize = N;
			
			SCVec1DZMM8r8() noexcept(true) {

				m_Re[N];
				m_Im[N];
			}

			SCVec1DZMM8r8(const double (&Re)[N],
				      const double (&Im)[N]
							  ) {
				using namespace gms::common;
				
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1	
				avx512_uncached_memmove(&m_Re[0],Re,N);	
				avx512_uncached_memmove(&m_Im[0],Im,N);
#else
				avx512_cached_memmove(&m_Re[0],Re,N);
				avx512_cached_memmove(&m_Im[0],Im,N);
#endif					
			}

			SCVec1DZMM8r8(const SCVec1DZMM8r8 &x)  {
				using namespace gms::common;
				m_nsize = x.m_nsize;
				
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
				avx512_uncached_memmove(&m_Re[0], &x.m_Re[0], N);
				avx512_uncached_memmove(&m_Im[0], &x.m_Im[0], N);
#else
				avx512_cached_memmove(&m_Re[0], &x.m_Re[0], N);
				avx512_cached_memmove(&m_Im[0], &x.m_Im[0], N);
#endif
			}

		
		    // operator= not needed here, although implemented
			SCVec1DZMM8r8 & 
			operator=(const SCVec1DZMM8r8 &x){
			    using namespace GMS::common;
				if (this == &x || m_nsize != x.m_nsize){
				    return (*this);
				}
				// Destructive copy
				
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
				avx512_uncached_memmove(&m_Re[0],&x.m_Re[0], N);
				avx512_cached_memmove(&m_Im[0],&x.m_Im[0], N);
#else
				avx512_cached_memmove(&m_Re[0], &x.m_Re[0], N);
				avx512_cached_memmove(&m_Im[0], &x.m_Im[0], N);
#endif
				return (*this);
			}

			inline int32_t  size() { return (m_nsize); };

		};

		// Global operators

		template<int32_t N> std::ostream & 
		operator<<(std::ostream &os,
			   const SCVec1DZMM8r8<N> &x) {
			for (int32_t i = 0; i != x.m_nsize; ++i) {
				os << std::fixed << std::showpoint << std::setprecision(15) <<
					std::setw(4)  <<  "Re: " << "{"  << x.m_Re[i] << "}" <<
					std::setw(12) << "Im: "  << "{"  << x.m_Im[i] << "}" << std::endl;
			}
			return (os);
		}

		template<int32_t N> SCVec1DZMM8r8<N>
		inline operator+(const SCVec1DZMM8r8<N> &x,
				 const SCVec1DZMM8r8<N> &y) {
			if (x.size() != y.size()) { return (SCVec1DZMM8r8<N>{}); }
			SCVec1DZMM8r8<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.size(), 8); i += 16) {
				// Linearly growing indices, no need for software prefetch.
				// HW prefetch will kick in after 2 maybe 3 cache misses.
				const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
				const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i + 0]);
				_mm512_stream_pd(&ret_vec.m_Re[i], _mm512_add_pd(zmm0, zmm1));
				const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
				const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 8]);
				_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
				const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
				const __m512d zmm5 = _mm512_load_pd(&y.m_Im[i + 0]);
				_mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(zmm4, zmm5));
				const __m512d zmm6 = _mm512_load_pd(&x.m_Im[i + 8]);
				const __m512d zmm7 = _mm512_load_pd(&y.m_Im[i + 8]);
				_mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(zmm6, zmm7));
			}
			_mm_sfence();
			for (; i != ret_vec.m_nsize; ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
				ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
			}
#else			
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.size(), 8); i += 16) {
				// Linearly growing indices, no need for software prefetch.
				// HW prefetch will kick in after 2 maybe 3 cache misses.
				const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
				const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i + 0]);
				_mm512_store_pd(&ret_vec.m_Re[i], _mm512_add_pd(zmm0, zmm1));
				const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
				const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 8]);
				_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
				const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
				const __m512d zmm5 = _mm512_load_pd(&y.m_Im[i + 0]);
				_mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(zmm4, zmm5));
				const __m512d zmm6 = _mm512_load_pd(&x.m_Im[i + 8]);
				const __m512d zmm7 = _mm512_load_pd(&y.m_Im[i + 8]);
				_mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(zmm6, zmm7));
			}
			for (; i != ret_vec.m_nsize; ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
				ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
			}
#endif				
			return (ret_vec);
				
		}
				


			




			

				


				

				


				


			
			
			
		

		template<int32_t N> SCVec1DZMM8r8<N>
		inline   operator+(const SCVec1DZMM8r8<N> &x,
				  const double  Re[N]) {         // If Re is not equal to x --> udefined behaviour.
				  using namespace gms::common;
				  if (!Is_ptr_aligned64(Re)) { return SCVec1DZMM8r8<N>{}; }
				  SCVec1DZMM8r8<N> ret_vec;
				  int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
					  // Linearly growing indices, no need for software prefetch.
					  // HW prefetch will kick in after 2 maybe 3 cache misses.
					  const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
					  const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					  _mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
					  const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					  _mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
					} 
				  _mm_sfence();
				  for (; i != ret_vec.size(); ++i) {
					  ret_vec.m_Re[i] = x.m_Re[i] + Re[i];
				  }
#else
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
					  // Linearly growing indices, no need for software prefetch.
					  // HW prefetch will kick in after 2 maybe 3 cache misses.

					  const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
					  const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					  _mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
					  const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					  _mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));

				  }
				  for (; i != ret_vec.size(); ++i) {
					  ret_vec.m_Re[i] = x.m_Re[i] + Re[i];
				  }
#endif
				  return (ret_vec);
		 }
				 
				 
		

		template<int32_t N> SCVec1DZMM8r8<N>
		inline operator+(const double  Re[N],
				 const SCVec1DZMM8r8<N> &x) {
				  using namespace gms::common;
				  if (!Is_ptr_aligned64(Re)) { return (SCVec1DZMM8r8<N>{}); }
				  SCVec1DZMM8r8<N> ret_vec;
				  int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					  const __m512d zmm0(_mm512_load_pd(&Re[i + 0]));
					  const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
					  _mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2(_mm512_load_pd(&Re[i + 8]));
					  const __m512d zmm3(_mm512_load_pd(&x.m_Re[i + 8]));
					  _mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
					}
				  _mm_sfence();
				  for (; i != x.size(); ++i) {
					  ret_vec.m_Re[i] = Re[i] + x.m_Re[i];
				  }
#else
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
					  const __m512d zmm0(_mm512_load_pd(&Re[i + 0]));
					  const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
					  _mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2(_mm512_load_pd(&Re[i + 8]));
					  const __m512d zmm3(_mm512_load_pd(&x.m_Re[i + 8]));
					  _mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
					} 
				  for (; i != x.size(); ++i) {
					  ret_vec.m_Re[i] = Re[i] + x.m_Re[i];
				  }
#endif				  
				  return (ret_vec);
	       }
				 
				 
				
		  

		  template<int32_t> SCVec1DZMM8r8<N>
		  inline operator+(const SCVec1DZMM8r8<N> &x,
				   const double c) {
					SCVec1DZMM8r8<N> ret_vec;
					int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0(_mm512_set1_pd(c));
						const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
						_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
						const __m512d zmm2(_mm512_load_pd(&x.m_Re[i+8]));
						_mm512_stream_pd(&ret_vec.m_Re[i+8],_mm512_add_pd(zmm0,zmm2));
						}
					_mm_sfence();
					for (; i != x.size(); ++i) {
						ret_vec.m_Re[i] = x.m_Re[i] + c;
					}
#else
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0(_mm512_set1_pd(c));
						const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
						_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
						const __m512d zmm2(_mm512_load_pd(&x.m_Re[i + 8]));
						_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm0, zmm2));
					}

					for (; i != x.size(); ++i) {
						ret_vec.m_Re[i] = x.m_Re[i] + c;
					}
#endif
					return (ret_vec);
			}		
		   

		   template<int32_t N> SCVec1DZMM8r8<N>
		   inline operator-(const SCVec1DZMM8r8<N> &x,
				    const SCVec1DZMM8r8<N> &y) {
					 using namespace gms::common;
					 if (x.size() != y.size()) { return (SCVec1DZMM8r8<N>{}); }
					 SCVec1DZMM8r8<N> ret_vec;
					 int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
					 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						 const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
						 const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i + 0]);
						 _mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
						 const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
						 const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 8]);
						 _mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
						 const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
						 const __m512d zmm5 = _mm512_load_pd(&y.m_Im[i + 0]);
						 _mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_sub_pd(zmm4, zmm5));
						 const __m512d zmm6 = _mm512_load_pd(&x.m_Im[i + 8]);
						 const __m512d zmm7 = _mm512_load_pd(&y.m_Im[i + 8]);
						 _mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_sub_pd(zmm6, zmm7));
					 }
					 _mm_sfence();
					 for (; i != x.size(); ++i) {
						 ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
						 ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
					 }
#else
						 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

							 const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
							 const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i + 0]);
							 _mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
							 const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
							 const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 8]);
							 _mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
							 const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
							 const __m512d zmm5 = _mm512_load_pd(&y.m_Im[i + 0]);
							 _mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_sub_pd(zmm4, zmm5));
							 const __m512d zmm6 = _mm512_load_pd(&x.m_Im[i + 8]);
							 const __m512d zmm7 = _mm512_load_pd(&y.m_Im[i + 8]);
							 _mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_sub_pd(zmm6, zmm7));
						}
						 for (; i != x.size(); ++i) {
							 ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
							 ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
						 }
#endif
					 return (ret_vec);
			}		 
					
		  

		  template<int32_t N> SCVec1DZMM8r8<N>
		  inline operator-(SCVec1DZMM8r8<N> &x,
				   const double Re[N]) {
					using namespace gms::common;
					if (!Is_ptr_aligned64(Re)) { return (SCVec1DZMM8r8<N>{}); }
					SCVec1DZMM8r8<N> ret_vec;
					int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
						const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
						_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
						const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
						const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
						_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));

					}
					_mm_sfence();
					for (; i != x.size(); ++i) {
						ret_vec.m_Re[i] = x.m_Re[i] - Re[i];
					}
#else
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
						const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
						_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
						const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
						const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
						_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
					}

					for (; i != x.size(); ++i) {
						ret_vec.m_Re[i] = x.m_Re[i] - Re[i];
					}
#endif
					
					return (ret_vec);
			}

			template<int32_t N> SCVec1DZMM8r8<N>
			inline operator-(  const SCVec1DZMM8r8<N> &x,
					   const double c) {
					  SCVec1DZMM8r8<N> ret_vec;
					  int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
					  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						  const __m512d zmm0(_mm512_set1_pd(c));
						  const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
						  _mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm1, zmm0));
						  const __m512d zmm2(_mm512_load_pd(&x.m_Re[i+8]));
						  _mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2,zmm0));
						  }
					  _mm_sfence();
					  for (; i != x.size(); ++i) {
						  ret_vec.m_Re[i] = x.m_Re[i] - c;
					  }
#else
					  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						  const __m512d zmm0(_mm512_set1_pd(c));
						  const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
						  _mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm1, zmm0));
						  const __m512d zmm2(_mm512_load_pd(&x.m_Re[i + 8]));
						  _mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm0));
					  }

					  for (; i != x.size(); ++i) {
						  ret_vec.m_Re[i] = x.m_Re[i] - c;
					  }
#endif
					 
					  return (ret_vec);
		   }

		   template<int32_t N> SCVec1DZMM8r8<N>
		   inline operator*(const SCVec1DZMM8r8<N> &x,
				    const SCVec1DZMM8r8<N> &y) {
			   if (x.size() != y.size()) { return (SCVec1DZMM8r8<N>{}); }
			    __m512d zmm0,zmm1,zmm2,zmm3,zmm4,zmm5,zmm6,zmm7;
			   SCVec1DZMM8r8<N> ret_vec;
			   int32_t i; 
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
			   for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm0, zmm1, zmm2, zmm3, x, y, i, 0)
					   _mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(
					   _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm3)));
				   _mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(
					   _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));

				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm4, zmm5, zmm6, zmm7, x, y, i, 8)
					   _mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(
					   _mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6, zmm7)));
				   _mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(
					   _mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));
			   }    
			   _mm_sfence();
			   for (; i != x.size(); ++i) {
				   ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
				   ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
			   }
#else
			   for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm0, zmm1, zmm2, zmm3, x, y, i, 0)
					   _mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(
					   _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm3)));
				   _mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(
					   _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));

				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm4, zmm5, zmm6, zmm7, x, y, i, 8)
					   _mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(
					   _mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6, zmm7)));
				   _mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(
					   _mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));
				}
			   for (; i != x.size(); ++i) {
				   ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
				   ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
			   }
#endif				 
			   return (ret_vec);
		}	  

			 
			  
		 


			  
			 
			 

			 template<int32_t N> SCVec1DZMM8r8<N>
			 inline operator*(const SCVec1DZMM8r8<N> &x,
					  const double  Re[N]) {
				 using namespace gms::common;
				 if (!Is_ptr_aligned64(Re)) { return (SCVec1DZMM8r8<N>{}); }
				 SCVec1DZMM8r8<N> ret_vec;
				 int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
				 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					 const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
					 const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					 _mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
					 const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
					 const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					 _mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
					 const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
					 _mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
					 const __m512d zmm5 = _mm512_load_pd(&x.m_Im[i + 8]);
					 _mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
					 }
				 _mm_sfence();
				 for (; i != ret_vec.m_nsize; ++i) {
					 ret_vec.m_Re[i] = x.m_Re[i] * Re[i];
					 ret_vec.m_Im[i] = x.m_Im[i] * Re[i];
				 }
#else
				 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					 const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
					 const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					 _mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
					 const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
					 const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					 _mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
					 const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
					 _mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
					 const __m512d zmm5 = _mm512_load_pd(&x.m_Im[i + 8]);
					 _mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
				 }
				 for (; i != ret_vec.m_nsize; ++i) {
					 ret_vec.m_Re[i] = x.m_Re[i] * Re[i];
					 ret_vec.m_Im[i] = x.m_Im[i] * Re[i];
				 }
#endif
				 return (ret_vec);
		 }
				
				 
			

			template<int32_t N> SCVec1DZMM8r8<N>
			inline operator*(const SCVec1DZMM8r8<N> &x,
					 const double c) {
				SCVec1DZMM8r8<N> ret_vec;
				int32_t i;
#if (USE_AVX512COMPLEX_SMALLV_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					const __m512d zmm0(_mm512_set1_pd(c));
					const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
					_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_mul_pd(zmm1, zmm0));
					const __m512d zmm2(_mm512_load_pd(&x.m_Re[i+8]));
					_mm512_stream_pd(&ret_vec.m_Re[i+8], _mm512_mul_pd(zmm2,zmm0));

				}
				_mm_sfence();
				for (; i != x.size(); ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] * c;
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					const __m512d zmm0(_mm512_set1_pd(c));
					const __m512d zmm1(_mm512_load_pd(&x.m_Re[i + 0]));
					_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_mul_pd(zmm1, zmm0));
					const __m512d zmm2(_mm512_load_pd(&x.m_Re[i + 8]));
					_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm0));
				}

				for (; i != x.size(); ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] * c;
				}
#endif
				return (ret_vec);
			}	
		  

		  template<int32_t N> SCVec1DZMM8r8<N>
		  inline operator/(const SCVec1DZMM8r8<N> &x,
				   const SCVec1DZMM8r8<N> &y) {
			  if (x.size() != y.size()) { return (SCVec1DZMM8r8<N>{}); }
			  SCVec1DZMM8r8<N> ret_vec;
			  int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
			  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
				  // Will unrolling 2x not saturate divider unit.
				  // We have two parallel division so at least second
				  // operation will be pipelined at divider level.

				  const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
				  const __m512d zmm1 = _mm512_load_pd(&y.m_Im[i + 0]);
				  const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i + 0]);
				  const __m512d re_term1 = _mm512_add_pd(
					  _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm1));
				  const __m512d re_term2 = _mm512_add_pd(
					  _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm1));
				  const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 0]);
				  const __m512d den_term = _mm512_add_pd(
					  _mm512_mul_pd(zmm3, zmm3), _mm512_mul_pd(zmm1, zmm1));

				  _mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
				  _mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));
				  }
			  _mm_sfence();
			  for (; i != x.size(); ++i) {
				  const double tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
				  const double tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
				  const double den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
				  ret_vec.m_Re[i] = tre / den;
				  ret_vec.m_Im[i] = tim / den;
			  }
#else
			  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
				  // Will unrolling 2x not saturate divider unit.
				  // We have two parallel division so at least second
				  // operation will be pipelined at divider level.

				  const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
				  const __m512d zmm1 = _mm512_load_pd(&y.m_Im[i + 0]);
				  const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i + 0]);
				  const __m512d re_term1 = _mm512_add_pd(
					  _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm1));
				  const __m512d re_term2 = _mm512_add_pd(
					  _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm1));
				  const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 0]);
				  const __m512d den_term = _mm512_add_pd(
					  _mm512_mul_pd(zmm3, zmm3), _mm512_mul_pd(zmm1, zmm1));

				  _mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
				  _mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));
			  }
			  for (; i != x.size(); ++i) {
				  const double tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
				  const double tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
				  const double den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
				  ret_vec.m_Re[i] = tre / den;
				  ret_vec.m_Im[i] = tim / den;
			  }
#endif
			  return (ret_vec);
		}
			 
			 
		 	 
			 
		 template<int32_t N> SCVec1DZMM8r8<N>
		 inline operator/(const SCVec1DZMM8r8<N> &x,
				  const double  Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned64(Re)) { return (SCVec1DZMM8r8<N>{}); }
				SCVec1DZMM8r8<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_ZMM8R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
					// Will unrolling 2x not saturate divider unit.
					// We have two parallel division so at least second
					// operation will be pipelined at divider level.
					const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i]);
					const __m512d zmm1 = _mm512_load_pd(&Re[i]);
					_mm512_stream_pd(&ret_vec.m_Re[i], _mm512_div_pd(zmm0, zmm1));
					const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i]);
					_mm512_stream_pd(&ret_vec.m_Im[i], _mm512_div_pd(zmm2, zmm1));
				}	
				_mm_sfence();
				for (; i != x.size(); ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] / Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] / Re[i];
				}	
#else
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
					// Will unrolling 2x not saturate divider unit.
					// We have two parallel division so at least second
					// operation will be pipelined at divider level.
					const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i]);
					const __m512d zmm1 = _mm512_load_pd(&Re[i]);
					_mm512_store_pd(&ret_vec.m_Re[i], _mm512_div_pd(zmm0, zmm1));
					const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i]);
					_mm512_store_pd(&ret_vec.m_Im[i], _mm512_div_pd(zmm2, zmm1));
				}
				for (; i != x.size(); ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] / Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] / Re[i];
				}
#endif
				return (ret_vec);
			}
				
					
				

				
				
				
		  

		 
		  template<int32_t N>
		  std::pair<std::array<__mmask8,Size>,
		  std::array<__mmask8,Size>>
		  inline   operator==(const SCVec1DZMM8r8<N> &x,
				      const SCVec1DZMM8r8<N> &y) {
			  using namespace gms::common;
			  
			 
			  if (x.size() != y.size()) {
				  return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
			  }
			  int32_t i;
			  int32_t k1 = 0, k2 = 0;
			  std::pair<std::array<__mmask8,N>,
						std::array<__mmask8,N>> ret_val;
			  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
				  const __m512d zmm0   = _mm512_load_pd(&x.m_Re[i]);
				  const __m512d zmm1   = _mm512_load_pd(&y.m_Re[i]);
				  ret_val.first[++k1]  = _mm512_cmp_pd_mask(zmm0, zmm1, _CMP_EQ_OQ);
				  const __m512d zmm2   = _mm512_load_pd(&x.m_Im[i]);
				  const __m512d zmm3   = _mm512_load_pd(&y.m_Im[i]);
				  ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2, zmm3, _CMP_EQ_OQ);
			  }
			  unsigned char t1 = 0x00, t2 = 0x00;
			  k1 += 1;
			  k2 += 1;
			  for (; i != x.size(); ++i) {
				  if (approximately_equalf64(x.m_Re[i],
					  y.m_Re[i],std::numeric_limits<double>::epsilon())) {
					  t1 |= 1 << i;
				  }
				  if (approximately_equalf64(x.m_Im[i],
					  y.m_Im[i],std::numeric_limits<double>::epsilon())) {
					  t2 |= 1 << i
				  }
			  }
			  ret_val.first[k1] = t1;
			  ret_val.second[k2] = t2;
			  return (ret_val); // Move Ctor should kick in if compiler will be smart enough
							    // or at least fast memcpy of primitive type
			}

			template<int32_t N>
			std::pair<std::array<__mmask8,N>,
			std::array<__mmask8,N>>
		    inline operator!=(const SCVec1DZMM8r8<N> &x,
				      const SCVec1DZMM8r8<N> &y) {
				using namespace gms::common;
				
				std::pair<std::array<__mmask8,N>,
				std::array<__mmask8,N>> ret_val;
				if (x.size() != y.size()) {
					return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
				}
				int32_t i;
				int32_t k1 = 0, k2 = 0;
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
					const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i]);
					const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i]);
					ret_val.first[++k1] = _mm512_cmp_pd_mask(zmm0, zmm1, _CMP_NEQ_OQ);
					const __m512 zmm2 = _mm512_load_pd(&x.m_Im[i]);
					const __m512 zmm3 = _mm512_load_pd(&y.m_Im[i]);
					ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2, zmm3, _CMP_NEQ_OQ);
				}
				unsigned char t1 = 0x00, t2 = 0x00;
				k1 += 1;
				k2 += 1;
				for (; i != x.size(); ++i) {
					if (!approximately_equalf64(x.m_Re[i],
						y.m_Re[i],std::numeric_limits<double>::epsilon())) {
						t1 |= 1 << i;
					}
					if (!approximately_equalf64(x.m_Im[i],
						y.m_Im[i],std::numeric_limits<double>::epsilon())) {
						t2 |= 1 << i;
					}
				}
				ret_val.first[k1] = t1;
				ret_val.second[k2] = t2;
				return (ret_val);
			}

			template<int32_t N>
			void cnormalize_product_zmm8r8(SCVec1DZMM8r8<N> &vout,
						     const SCVec1DZMM8r8<N> &v1,
						     const SCVec1DZMM8r8<N> &v2,
						     const bool do_nt_stream) {

				avx512_cnormalize_prod<SCVec1DZMM8r8<N>>(vout,v1,v2,do_nt_stream);
			}

			template<int32_t N>
			void cmean_product_zmm8r8(std::complex<double> &mean,
						const SCVec1DZMM8r8<N> &v1,
						const SCVec1DZMM8r8<N> &v2) {
				
				avx512_cmean_prod<SCVec1DZMM8r8<N>>(mean,v1,v2);
			}

			template<int32_t N>
			void cmean_quotient_zmm8r8(std::complex<double> &mean,
						 const SCVec1DZMM8r8<N> &v1,
						 const SCVec1DZMM8r8<N> &v2) {

				avx512_cmean_quot<SCVec1DZMM8r8<N>>(mean,v1,v2);
			}

			template<int32_t N>
			void cconj_product_zmm8r8(SCVec1DZMM8r8<N> &vout,
						const SCVec1DZMM8r8<N> &v1,
						const SCVec1DZMM8r8<N> &v2,
						const bool do_nt_store) {
				
				avx512_cconj_prod<SCVec1DZMM8r8<N>>(vout,v1,v2,do_nt_store);
			}

			template<int32_t N>
			void cnorm_conjprod_zmm8r8(SCVec1DZMM8r8<N> &vout,
						 const SCVec1DZMM8r8<N> &v1,
						 const SCVec1DZMM8r8<N> &v2,
						 const bool do_nt_store) {
				
				avx512_cnorm_conjprod<SCVec1DZMM8r8<N>>(vout,v1,v2,do_nt_store);
			}

			template<int32_t N>
			void cmean_conjprod_zmm8r8(std::complex<double> &mean,
						 const SCVec1DZMM8r8<N> &v1,
						 const SCVec1DZMM8r8<N> &v2) {
				
				avx512_cmean_conjprod<SCVec1DZMM8r8<N>>(mean,v1,v2);
			}

			template<int32_t N>
			void carith_mean_zmm8r8(std::complex<double> &mean,
					      const SCVec1DZMM8r8<N> &v) {
				
				avx512_arith_mean<SCVec1DZMM8r8<N>>(mean,v);
			}

			template<int32_t N>
			void cnormalize_zmm8r8(SCVec1DZMM8r8<N> &norm,
					      const SCVec1DZMM8r8<N> &v,
					      const SCVec1DZMM8r8<N> &cv,
					      const bool do_nt_store) {
					
				avx512_cnormalize<SCVec1DZMM8r8<N>>(norm,v,cv,do_nt_store);
			}

			template<int32_t N>
			void cmagnitude_zmm8r8(const double * __restrict vmag,
					      const SCVec1DZMM8r8<N> &v,
					      const SCVec1DZMM8r8<N> &cv) {
				
				avx512_cmagnitude<SCVec1DZMM8r8<N>>(vmag,v,cv);
			}

	}
}




#endif /*__GMS_STATIC_CVEC1D_ZMM8R8_H__*/
