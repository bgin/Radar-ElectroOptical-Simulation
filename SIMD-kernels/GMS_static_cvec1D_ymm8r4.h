
#ifndef __GMS_STATIC_CVEC1D_YMM8R4_H__
#define __GMS_STATIC_CVEC1D_YMM8R4_H__ 061020191358

namespace file_info {

	
	const unsigned int GMS_STATIC_CVEC1D_YMM8R4_MAJOR = 1U;

	const unsigned int GMS_STATIC_CVEC1D_YMM8R4_MINOR = 0U;

	const unsigned int GMS_STATIC_CVEC1D_YMM8R4_MICRO = 0U;

	const unsigned int gGMS_STATIC_CVEC1D_YMM8R4_FULLVER = 
		1000U*GMS_STATIC_CVEC1D_YMM8R4_MAJOR+100U*GMS_STATIC_CVEC1D_YMM8R4_MINOR+10U*GMS_STATIC_CVEC1D_YMM8R4_MICRO;

	const char * const GMS_STATIC_CVEC1D_YMM8R4_CREATE_DATE = "06-10-2019 13:58 + 00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const GMS_STATIC_CVEC1D_YMM8R4_BUILD_DATE = __DATE__ ":" __TIME__;

	const char * const GMS_STATIC_CVEC1D_YMM8R4_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const GMS_STATIC_CVEC1D_YMM8R4_SYNOPSIS = "AVX complex vector (1D) stack-allocated storage.";
}

#include <cstdint>
#include <iostream>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_common.h"
#include "GMS_complex_common_ymm8r4.h"

#if !defined (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) // Streaming stores defined per this struct (default set to 0)
#define USE_STATIC_CVEC1D_YMM8R4_NT_STORES 0
#endif

namespace gms {
	namespace math {

#if !defined (STATIC_CVEC1D_YMM8R4_LOAD_YMM)
#define STATIC_CVEC1D_YMM8R4_LOAD_YMM(reg1,reg2,reg3,reg4,v1,v2,idx,off) \
	(reg1) = _mm256_load_ps(&(v1).m_Re[(idx)+(off)]);            \
	(reg2) = _mm256_load_ps(&(v2).m_Re[(idx)+(off)]);            \
	(reg3) = _mm256_load_ps(&(v1).m_Im[(idx)+(off)]);			  \
	(reg4) = _mm256_load_ps(&(v2).m_Im[(idx)+(off)]);
#endif

		
		template<int32_t N>
		struct __ATTR_ALIGN__(32) SVec1DYMM8r4{

                       __ATTR_ALIGN__(32) float m_Re[(m_nsize == 0) ? 4 : N];
		       __ATTR_ALIGN__(32) float m_Im[(m_nsize == 0) ? 4 : N];
                        int32_t m_nsize = N;
			
			SVec1DYMM8r4() noexcept(true) {
				m_Re[N];
				m_Im[N];
			}

			SVec1DYMM8r4(const float Re[N],
				       const float Im[N]) {
				using namespace gms::common;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				avx256_uncached_memmove(&m_Re[0], &Re[0], N);
				avx256_uncached_memmove(&m_Im[0], &Im[0], N);
#else
				avx256_cached_memmove(&m_Re[0], &Re[0], N);
				avx256_cached_memmove(&m_Im[0], &Im[0], N);
#endif
			}

			SVec1DYMM8r4(const SVec1DYMM8r4 &x) {
				using namespace gms::common;
				m_nsize = x.m_nsize;
#if (USE_AVXCOMPLEX_SMALLV_NT_STORES) == 1
				avx256_uncached_memmove(&m_Re[0], &x.m_Re[0], x.m_nsize);
				avx256_uncached_memmove(&m_Im[0], &x.m_Im[0], x.m_nsize);
#else
				avx256_cached_memmove(&m_Re[0], &x.m_Re[0], x.m_nsize);
				avx256_cached_memmove(&m_Im[0], &x.m_Im[0], x.m_nsize);
#endif
			}

			~SVec1DYMM8r4() = default;

			SVec1DYMM8r4 &
			operator=(const SVec1DYMM8r4 &x) {
				using namespace gms::common;
				if (this == &x || m_nsize != x.m_nsize) 
				    { return (*this); }
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				avx256_uncached_memmove(&m_Re[0], &x.m_Re[0],x.m_nsize);
				avx256_uncached_memmove(&m_Im[0], &x.m_Im[0],x.m_nsize);
#else
				avx256_cached_memmove(&m_Re[0], &x.m_Re[0], x.m_nsize);
				avx256_cached_memmove(&m_Im[0], &x.m_Im[0], x.m_nsize);
#endif				
				return (*this);
				}
		};

		template<int32_t N> std::ostream &
		operator<<(std::ostream &os,
			   const SVec1DYMM8r4<N> &x) {
			for (int32_t i = 0; i != x.m_nsize; ++i) {
				os << std::fixed << std::showpoint << std::setprecision(7) <<
					std::setw(4)  << "Re: " << "{" << x.m_Re[i] << "}" <<
					std::setw(12) << "Im: " << "{" << x.m_Im[i] << "}" << std::endl;
			}
			 return (os);
		 }


		 template<int32_t N> SVec1DYMM8r4<N>
		 inline operator+(const SVec1DYMM8r4<N> &x,
				  const SVec1DYMM8r4<N> &y) {
			 if (x.m_nsize != y.m_nsize) 
			     { return (SVec1DYMM8r4<N>{}); }
				 SVec1DYMM8r4<N> ret_vec;
				 int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				 for (i = 0; i != ROUND_TO_EIGHT(x.m_nsize, 8); i += 16) {
					 // Linearly growing indices, no need for software prefetch.
					 // HW prefetch will kick in after 2 maybe 3 cache misses.
					 const __m256 ymm0(_mm256_load_ps(&x.m_Re[i + 0]));
					 const __m256 ymm1(_mm256_load_ps(&y.m_Re[i + 0]));
					 _mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					 const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
					 const __m256 ymm3(_mm256_load_ps(&y.m_Re[i + 8]));
					 _mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm2, ymm3));
					 const __m256 ymm4(_mm256_load_ps(&x.m_Im[i + 0]));
					 const __m256 ymm5(_mm256_load_ps(&y.m_Im[i + 0]));
					 _mm256_stream_ps(&ret_vec.m_Im[i + 0], _mm256_add_ps(ymm4, ymm5));
					 const __m256 ymm6(_mm256_load_ps(&x.m_Im[i + 8]));
					 const __m256 ymm7(_mm256_load_ps(&y.m_Im[i + 8]));
					 _mm256_stream_ps(&ret_vec.m_Im[i + 8], _mm256_add_ps(ymm6, ymm7));
				}	
				 _mm_sfence();
				 for (; i != ret_vec.m_nsize; ++i) {
					 ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
					 ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
				 }
#else					
				 for (i = 0; i != ROUND_TO_EIGHT(x.m_nsize, 8); i += 16) {
					 // Linearly growing indices, no need for software prefetch.
					 // HW prefetch will kick in after 2 maybe 3 cache misses.
					 const __m256 ymm0(_mm256_load_ps(&x.m_Re[i + 0]));
					 const __m256 ymm1(_mm256_load_ps(&y.m_Re[i + 0]));
					 _mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					 const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
					 const __m256 ymm3(_mm256_load_ps(&y.m_Re[i + 8]));
					 _mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm2, ymm3));
					 const __m256 ymm4(_mm256_load_ps(&x.m_Im[i + 0]));
					 const __m256 ymm5(_mm256_load_ps(&y.m_Im[i + 0]));
					 _mm256_store_ps(&ret_vec.m_Im[i + 0], _mm256_add_ps(ymm4, ymm5));
					 const __m256 ymm6(_mm256_load_ps(&x.m_Im[i + 8]));
					 const __m256 ymm7(_mm256_load_ps(&y.m_Im[i + 8]));
					 _mm256_store_ps(&ret_vec.m_Im[i + 8], _mm256_add_ps(ymm6, ymm7));
				 }
				 for (; i != ret_vec.m_nsize; ++i) {
					 ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
					 ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
				 }
#endif					
				 return (ret_vec)
			}		

				 
				



				
					 

					 

					

				 
				 
				 
			 

			 template<int32_t N> SVec1DYMM8r4<N>	
			 inline operator+(const SVec1DYMM8r4<N> &x, 
					  const float __restrict Re[N]) {   // If Re is not equal to x --> udefined behaviour.
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM8r4<N>{}); }
				SVec1DYMM8r4<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&Re[i + 0]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&Re[i + 8]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&Re[i + 0]);
					_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&Re[i + 8]);
					_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm2, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM8r4<N>
		        inline operator+(const float Re[N],
					 const SVec1DYMM8r4<N> &x) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM8r4<N>{}); }
				SVec1DYMM8r4<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0(_mm256_load_ps(&Re[i + 0]));
					const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
					_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					const __m256 ymm2(_mm256_load_ps(&Re[i + 8]));
					const __m256 ymm3(_mm256_load_ps(&x.m_Re[i + 8]));
					_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = Re[i] + x.m_Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0(_mm256_load_ps(&Re[i + 0]));
					const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
					_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					const __m256 ymm2(_mm256_load_ps(&Re[i + 8]));
					const __m256 ymm3(_mm256_load_ps(&x.m_Re[i + 8]));
					_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm2, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = Re[i] + x.m_Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM8r4<N>
			inline operator+(SVec1DYMM8r4<N> &x,
					 const float c) {
				SVec1DYMM8r4<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0(_mm256_set1_ps(c));
					const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
					_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
					_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm0, ymm2));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + c;
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0(_mm256_set1_ps(c));
					const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
					_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
					const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
					_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm0, ymm2));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + c;
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM8r4<N>
			inline operator-( const SVec1DYMM8r4<N> &x,
					  const SVec1DYMM8r4<N> &y) {
				using namespace gms::common;
				if (x.m_nsize != y.m_nsize) { return (SVec1DYMM8r4<N>{}); }
				int32_t i;
				SVec1DYMM8r4<N> ret_vec;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&y.m_Re[i + 0]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&y.m_Re[i + 8]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(ymm2, ymm3));
					const __m256 ymm4 = _mm256_load_ps(&x.m_Im[i + 0]);
					const __m256 ymm5 = _mm256_load_ps(&y.m_Im[i + 0]);
					_mm256_stream_ps(&ret_vec.m_Im[i + 0], _mm256_sub_ps(ymm4, ymm5));
					const __m256 ymm6 = _mm256_load_ps(&x.m_Im[i + 8]);
					const __m256 ymm7 = _mm256_load_ps(&y.m_Im[i + 8]);
					_mm256_stream_ps(&ret_vec.m_Im[i + 8], _mm256_sub_ps(ymm6, ymm7));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&y.m_Re[i + 0]);
					_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&y.m_Re[i + 8]);
					_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(ymm2, ymm3));
					const __m256 ymm4 = _mm256_load_ps(&x.m_Im[i + 0]);
					const __m256 ymm5 = _mm256_load_ps(&y.m_Im[i + 0]);
					_mm256_store_ps(&ret_vec.m_Im[i + 0], _mm256_sub_ps(ymm4, ymm5));
					const __m256 ymm6 = _mm256_load_ps(&x.m_Im[i + 8]);
					const __m256 ymm7 = _mm256_load_ps(&y.m_Im[i + 8]);
					_mm256_store_ps(&ret_vec.m_Im[i + 8], _mm256_sub_ps(ymm6, ymm7));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM8r4<N>
			inline operator-(const SVec1DYMM8r4<N> &x,
					 const float  Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM8r4<N>{}) };
				SVec1DYMM8r4<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&Re[i + 0]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&Re[i + 8]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&Re[i + 0]);
					_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&Re[i + 8]);
					_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(ymm2, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM8r4<N>
			inline operator-(const SVec1DYMM8r4<N> &x,
					 const float c) {
				SVec1DYMM8r4<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0(_mm256_set1_ps(c));
					const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
					_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(ymm1, ymm0));
					const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
					_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(ymm2, ymm0));


				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - c;
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0(_mm256_set1_ps(c));
					const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
					_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(ymm1, ymm0));
					const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
					_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(ymm2, ymm0));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - c;
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM8r4<N>
			inline operator*(const SVec1DYMM8r4<N> &x,
					 const SVec1DYMM8r4<N> &y) {
				if (x.m_nsize != y.m_nsize) { return (SVec1DYMM8r4<N>{}); }
				SVec1DYMM8r4<N> ret_vec;

				__ATTR_ALIGN__(32) struct {
                                        __m256 ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
				}ca;

				ca.ymm0 = _mm256_setzero_ps();
				ca.ymm1 = _mm256_setzero_ps();
				ca.ymm2 = _mm256_setzero_ps();
				ca.ymm3 = _mm256_setzero_ps();
				ca.ymm4 = _mm256_setzero_ps();
				ca.ymm5 = _mm256_setzero_ps();
				ca.ymm6 = _mm256_setzero_ps();
				ca.ymm7 = _mm256_setzero_ps();
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm0, ca.ymm1, ca.ymm2, ca.ymm3, x, y, i, 0)
						_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(
						_mm256_mul_ps(ca.ymm0, ca.ymm1), _mm256_mul_ps(ca.ymm2, ca.ymm3)));
					_mm256_stream_ps(&ret_vec.m_Im[i + 0], _mm256_add_ps(
						_mm256_mul_ps(ca.ymm2, ca.ymm1), _mm256_mul_ps(ca.ymm0, ca.ymm3)));

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm4, ca.ymm5, ca.ymm6, ca.ymm7, x, y, i, 8)
						_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(
						_mm256_mul_ps(ca.ymm4, ca.ymm5), _mm256_mul_ps(ca.ymm6, ca.ymm7)));
					_mm256_stream_ps(&ret_vec.m_Im[i + 8], _mm256_add_ps(
						_mm512_mul_ps(ca.ymm6, ca.ymm5), _mm512_mul_ps(ca.ymm4, ca.ymm7)));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
					ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm0, ca.ymm1, ca.ymm2, ca.ymm3, x, y, i, 0)
						_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_sub_ps(
						_mm256_mul_ps(ca.ymm0, ca.ymm1), _mm256_mul_ps(ca.ymm2, ca.ymm3)));
					_mm256_store_ps(&ret_vec.m_Im[i + 0], _mm256_add_ps(
						_mm256_mul_ps(ca.ymm2, ca.ymm1), _mm256_mul_ps(ca.ymm0, ca.ymm3)));

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm4, ca.ymm5, ca.ymm6, ca.ymm7, x, y, i, 8)
						_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_sub_ps(
						_mm256_mul_ps(ca.ymm4, ca.ymm5), _mm256_mul_ps(ca.ymm6, ca.ymm7)));
					_mm256_store_ps(&ret_vec.m_Im[i + 8], _mm256_add_ps(
						_mm512_mul_ps(ca.ymm6, ca.ymm5), _mm512_mul_ps(ca.ymm4, ca.ymm7)));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
					ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
				}
#endif
				
				return (ret_vec);
		  }

		  template<int32_t N> SVec1DYMM8r4<N>
		  inline operator*(const SVec1DYMM8r4<N> &x,
				   const float   Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM8r4<N>{}); }
				SVec1DYMM8r4<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&Re[i + 0]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_mul_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&Re[i + 8]);
					_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_mul_ps(ymm2, ymm3));
					const __m256 ymm4 = _mm256_load_ps(&x.m_Im[i + 0]);
					_mm256_stream_ps(&ret_vec.m_Im[i + 0], _mm256_mul_ps(ymm4, ymm1));
					const __m256 ymm5 = _mm256_load_ps(&x.m_Im[i + 8]);
					_mm256_stream_ps(&ret_vec.m_Im[i + 8], _mm256_mul_ps(ymm5, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] * Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] * Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

					const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i + 0]);
					const __m256 ymm1 = _mm256_load_ps(&Re[i + 0]);
					_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_mul_ps(ymm0, ymm1));
					const __m256 ymm2 = _mm256_load_ps(&x.m_Re[i + 8]);
					const __m256 ymm3 = _mm256_load_ps(&Re[i + 8]);
					_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_mul_ps(ymm2, ymm3));
					const __m256 ymm4 = _mm256_load_ps(&x.m_Im[i + 0]);
					_mm256_store_ps(&ret_vec.m_Im[i + 0], _mm256_mul_ps(ymm4, ymm1));
					const __m256 ymm5 = _mm256_load_ps(&x.m_Im[i + 8]);
					_mm256_store_ps(&ret_vec.m_Im[i + 8], _mm256_mul_ps(ymm5, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] * Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] * Re[i];
				}
#endif
				
				return (ret_vec);
		  }

		  template<int32_t N> SVec1DYMM8r4<N>
		  inline operator*(const SVec1DYMM8r4<N> &x,
				   const float c) {
			SVec1DYMM8r4<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

				const __m256 ymm0(_mm256_set1_ps(c));
				const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
				_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_mul_ps(ymm1, ymm0));
				const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
				_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_mul_ps(ymm2, ymm0));

			}
			_mm_sfence();
			for (; i != ret_vec.m_nsize; ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] * c;
			}
#else
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {

				const __m256 ymm0(_mm256_set1_ps(c));
				const __m256 ymm1(_mm256_load_ps(&x.m_Re[i + 0]));
				_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_mul_ps(ymm1, ymm0));
				const __m256 ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
				_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_mul_ps(ymm2, ymm0));

			}
			for (; i != ret_vec.m_nsize; ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] * c;
			}
#endif
			
			return (ret_vec);
		}

		template<int32_t N> SVec1DYMM8r4<N>
		inline operator/(const SVec1DYMM8r4<N> &x,
				 const SVec1DYMM8r4<N> &y) {
			if (x.m_nsize != y.m_nsize) { return (SVec1DYMM8r4<N>{}); }
			SVec1DYMM8r4<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 8) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.

				const __m256 ymm0 ( _mm256_load_ps(&x.m_Re[i + 0]));
				const __m256 ymm1 ( _mm256_load_ps(&y.m_Im[i + 0]));
				const __m256 ymm2 ( _mm256_load_ps(&x.m_Im[i + 0]));
				const __m256 re_term1 ( _mm256_add_ps(
					_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2, ymm1)));
				const __m256 re_term2 ( _mm256_add_ps(
					_mm256_mul_ps(ymm2, ymm1), _mm256_mul_ps(ymm0, ymm1)));
				const __m256 ymm3 ( _mm256_load_ps(&y.m_Re[i + 0]));
				const __m256 den_term ( _mm256_add_ps(
					_mm256_mul_ps(ymm3, ymm3), _mm256_mul_ps(ymm1, ymm1)));

				_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_div_ps(re_term1, den_term));
				_mm256_stream_ps(&ret_vec.m_Im[i + 0], _mm256_div_ps(re_term2, den_term));
			}
			_mm_sfence();
			for (; i != ret_vec.m_nsize; ++i) {
				const float tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
				const float tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
				const float den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
				ret_vec.m_Re[i] = tre / den;
				ret_vec.m_Im[i] = tim / den;
			}
#else
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 8) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.

				const __m256 ymm0(_mm256_load_ps(&x.m_Re[i + 0]));
				const __m256 ymm1(_mm256_load_ps(&y.m_Im[i + 0]));
				const __m256 ymm2(_mm256_load_ps(&x.m_Im[i + 0]));
				const __m256 re_term1(_mm256_add_ps(
					_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2, ymm1)));
				const __m256 re_term2(_mm256_add_ps(
					_mm256_mul_ps(ymm2, ymm1), _mm256_mul_ps(ymm0, ymm1)));
				const __m256 ymm3(_mm256_load_ps(&y.m_Re[i + 0]));
				const __m256 den_term(_mm256_add_ps(
					_mm256_mul_ps(ymm3, ymm3), _mm256_mul_ps(ymm1, ymm1)));

				_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_div_ps(re_term1, den_term));
				_mm256_store_ps(&ret_vec.m_Im[i + 0], _mm256_div_ps(re_term2, den_term));
			}
			
			for (; i != ret_vec.m_nsize; ++i) {
				const float tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
				const float tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
				const float den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
				ret_vec.m_Re[i] = tre / den;
				ret_vec.m_Im[i] = tim / den;
			}
#endif		
			return (ret_vec);
		 }

		 template<int32_t N> SVec1DYMM8r4<N>
		 inline operator/(const SVec1DYMM8r4<N> &x,
				  const float  Re[N]) {
			using namespace gms::common;
			if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM8r4<N>{}); }
			SVec1DYMM8r4<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 8) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.
				const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i]);
				const __m256 ymm1 = _mm256_load_ps(&Re[i]);
				_mm256_stream_ps(&ret_vec.m_Re[i], _mm256_div_ps(ymm0, ymm1));
				const __m256 ymm2 = _mm256_load_ps(&x.m_Im[i]);
				_mm256_stream_ps(&ret_vec.m_Im[i], _mm256_div_ps(ymm2, ymm1));
			}
			_mm_sfence();
			for (; i != x.size(); ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] / Re[i];
				ret_vec.m_Im[i] = x.m_Im[i] / Re[i];
			}
#else
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 8) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.
				const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i]);
				const __m256 ymm1 = _mm256_load_ps(&Re[i]);
				_mm256_store_ps(&ret_vec.m_Re[i], _mm256_div_ps(ymm0, ymm1));
				const __m256 ymm2 = _mm256_load_ps(&x.m_Im[i]);
				_mm256_store_ps(&ret_vec.m_Im[i], _mm256_div_ps(ymm2, ymm1));
			}
			for (; i != x.size(); ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] / Re[i];
				ret_vec.m_Im[i] = x.m_Im[i] / Re[i];
			}
#endif
			return (ret_vec);
		}

		template<int32_t N> SVec1DYMM8r4<N>
		inline operator==(const SVec1DYMM8r4<N> &x,
				  const SVec1DYMM8r4<N> &y) {
			using namespace gms::common;
			if (x.m_nsize != y.m_nsize) { return (SVec1DYMM8r4<N>{}); }

			__ATTR_ALIGN__(32) struct {
                                __m256 ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
			}ca;

			ca.ymm0 = _mm256_setzero_ps(); 
			ca.ymm1 = _mm256_setzero_ps();
			ca.ymm2 = _mm256_setzero_ps();
			ca.ymm3 = _mm256_setzero_ps();
			ca.ymm4 = _mm256_setzero_ps();
			ca.ymm5 = _mm256_setzero_ps();
			ca.ymm6 = _mm256_setzero_ps();
			ca.ymm7 = _mm256_setzero_ps();
			SVec1DYMM8r4<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
			for(i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {
				 ca.ymm0(_mm256_load_ps(&x.m_Re[i+0]));
				 ca.ymm1(_mm256_load_ps(&y.m_Re[i+0]));
				 _mm256_stream_ps(&ret_vec.m_Re[i+0], _mm256_cmp_ps(ca.ymm0,ca.ymm1,_CMP_EQ_OQ));
				 ca.ymm2(_mm256_load_ps(&x.m_Re[i+8]));
				 ca.ymm3(_mm256_load_ps(&y.m_Re[i+8]));
				 _mm256_stream_ps(&ret_vec.m_Re[i+8], _mm256_cmp_ps(ca.ymm2,ca.ymm3,_CMP_EQ_OQ));
				 ca.ymm4(_mm256_load_ps(&x.m_Im[i+0]));
				 ca.ymm5(_mm256_load_ps(&y.m_Im[i+0]));
				 _mm256_stream_ps(&ret_vec.m_Im[i+0], _mm256_cmp_ps(ca.ymm4,ca.ymm5,_CMP_EQ_OQ));
				 ca.ymm6(_mm256_load_ps(&x.m_Im[i+8]));
				 ca.ymm7(_mm256_load_ps(&y.m_Im[i+8]));
				 _mm256_stream_ps(&ret_vec.m_Im[i+8],_mm256_cmp_ps(ca.ymm6,ca.ymm7,_CMP_EQ_OQ));
			}
			_mm_sfence();
			for(; i != ret_vec.m_nsize; ++i) {
				if(approximately_equalf32(x.m_Re[i],
					        y.m_Re[i],std::numeric_limits<float>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				 }
				 else {
					 ret_vec.m_Re[i] = 0.0;
				 }
				 if (approximately_equalf32(x.m_Im[i],
						y.m_Im[i],std::numeric_limits<float>::epsilon())) {
					 ret_vec.m_Im[i] = 1.0;
				}
				else {
					ret_vec.m_Im[i] = 0.0;
				}
			}
#else
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {
				ca.ymm0(_mm256_load_ps(&x.m_Re[i + 0]));
				ca.ymm1(_mm256_load_ps(&y.m_Re[i + 0]));
				_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_cmp_ps(ca.ymm0, ca.ymm1, _CMP_EQ_OQ));
				ca.ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
				ca.ymm3(_mm256_load_ps(&y.m_Re[i + 8]));
				_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_cmp_ps(ca.ymm2, ca.ymm3, _CMP_EQ_OQ));
				ca.ymm4(_mm256_load_ps(&x.m_Im[i + 0]));
				ca.ymm5(_mm256_load_ps(&y.m_Im[i + 0]));
				_mm256_store_ps(&ret_vec.m_Im[i + 0], _mm256_cmp_ps(ca.ymm4, ca.ymm5, _CMP_EQ_OQ));
				ca.ymm6(_mm256_load_ps(&x.m_Im[i + 8]));
				ca.ymm7(_mm256_load_ps(&y.m_Im[i + 8]));
				_mm256_store_ps(&ret_vec.m_Im[i + 8], _mm256_cmp_ps(ca.ymm6, ca.ymm7, _CMP_EQ_OQ));
			}
			
			for (; i != ret_vec.m_nsize; ++i) {
				if (approximately_equalf32(x.m_Re[i],
					y.m_Re[i], std::numeric_limits<float>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				}
				else {
					ret_vec.m_Re[i] = 0.0;
				}
				if (approximately_equalf32(x.m_Im[i],
					y.m_Im[i], std::numeric_limits<float>::epsilon())) {
					ret_vec.m_Im[i] = 1.0;
				}
				else {
					ret_vec.m_Im[i] = 0.0;
				}
		  }
#endif
				return (ret_vec);
		}

		template<int32_t N> SVec1DYMM8r4<N>
		inline operator!=(const SVec1DYMM8r4<N> &x,
				  const SVec1DYMM8r4<N> &y) {
			using namespace gms::common;
		        if (x.m_nsize != y.m_nsize) { return (SVec1DYMM8r4<N>{}); }

			__ATTR_ALIGN__(32) struct {
                                __m256 ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
			}ca;

			ca.ymm0 = _mm256_setzero_ps();
			ca.ymm1 = _mm256_setzero_ps();
			ca.ymm2 = _mm256_setzero_ps();
			ca.ymm3 = _mm256_setzero_ps();
			ca.ymm4 = _mm256_setzero_ps();
			ca.ymm5 = _mm256_setzero_ps();
			ca.ymm6 = _mm256_setzero_ps();
			ca.ymm7 = _mm256_setzero_ps();
			SVec1DYMM8r4<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM8R4_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {
				ca.ymm0(_mm256_load_ps(&x.m_Re[i + 0]));
				ca.ymm1(_mm256_load_ps(&y.m_Re[i + 0]));
				_mm256_stream_ps(&ret_vec.m_Re[i + 0], _mm256_cmp_ps(ca.ymm0, ca.ymm1, _CMP_NEQ_OQ));
				ca.ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
				ca.ymm3(_mm256_load_ps(&y.m_Re[i + 8]));
				_mm256_stream_ps(&ret_vec.m_Re[i + 8], _mm256_cmp_ps(ca.ymm2, ca.ymm3, _CMP_NEQ_OQ));
				ca.ymm4(_mm256_load_ps(&x.m_Im[i + 0]));
				ca.ymm5(_mm256_load_ps(&y.m_Im[i + 0]));
				_mm256_stream_ps(&ret_vec.m_Im[i + 0], _mm256_cmp_ps(ca.ymm4, ca.ymm5, _CMP_NEQ_OQ));
				ca.ymm6(_mm256_load_ps(&x.m_Im[i + 8]));
				ca.ymm7(_mm256_load_ps(&y.m_Im[i + 8]));
				_mm256_stream_ps(&ret_vec.m_Im[i + 8], _mm256_cmp_ps(ca.ymm6, ca.ymm7, _CMP_NEQ_OQ));
			}
			_mm_sfence();
			for (; i != ret_vec.m_nsize; ++i) {
				if (!approximately_equalf32(x.m_Re[i],
					y.m_Re[i], std::numeric_limits<float>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				}
				else {
					ret_vec.m_Re[i] = 0.0;
				}
				if (!approximately_equalf32(x.m_Im[i],
					y.m_Im[i], std::numeric_limits<float>::epsilon())) {
					ret_vec.m_Im[i] = 1.0;
				}
				else {
					ret_vec.m_Im[i] = 0.0;
				}
			}
#else
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {
				ca.ymm0(_mm256_load_ps(&x.m_Re[i + 0]));
				ca.ymm1(_mm256_load_ps(&y.m_Re[i + 0]));
				_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_cmp_ps(ca.ymm0, ca.ymm1, _CMP_NEQ_OQ));
				ca.ymm2(_mm256_load_ps(&x.m_Re[i + 8]));
				ca.ymm3(_mm256_load_ps(&y.m_Re[i + 8]));
				_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_cmp_ps(ca.ymm2, ca.ymm3, _CMP_NEQ_OQ));
				ca.ymm4(_mm256_load_ps(&x.m_Im[i + 0]));
				ca.ymm5(_mm256_load_ps(&y.m_Im[i + 0]));
				_mm256_store_ps(&ret_vec.m_Im[i + 0], _mm256_cmp_ps(ca.ymm4, ca.ymm5, _CMP_NEQ_OQ));
				ca.ymm6(_mm256_load_ps(&x.m_Im[i + 8]));
				ca.ymm7(_mm256_load_ps(&y.m_Im[i + 8]));
				_mm256_store_ps(&ret_vec.m_Im[i + 8], _mm256_cmp_ps(ca.ymm6, ca.ymm7, _CMP_NEQ_OQ));
			}

			for (; i != ret_vec.m_nsize; ++i) {
				if (!approximately_equalf32(x.m_Re[i],
					y.m_Re[i], std::numeric_limits<float>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				}
				else {
					ret_vec.m_Re[i] = 0.0;
				}
				if (!approximately_equalf32(x.m_Im[i],
					y.m_Im[i], std::numeric_limits<float>::epsilon())) {
					ret_vec.m_Im[i] = 1.0;
				}
				else {
					ret_vec.m_Im[i] = 0.0;
				}
			}
#endif
			return (ret_vec);
		  }


		  template<int32_t N>
		  void cnormalize_product_ymm8r4(SVec1DYMM8r4<N> &out,
					       const SVec1DYMM8r4<N> &v1,
					       const SVec1DYMM8r4<N> &v2) {
			  avx256_cnormalize_prod<SVec1DYMM8r4<N>>(out,v1,v2);
		  }

		  template<int32_t N>
		  void cmean_product_ymm8r4(std::complex<float> &mean,
					  const SVec1DYMM8r4<N> &v1,
					  const SVec1DYMM8r4<N> &v2) {
			  avx256_cmean_prod<SVec1DYMM8r4<N>>(mean,v1,v2);
		  }

		  template<int32_t N>
		  void cmean_quotient_ymm8r4( std::complex<float> &mean,
					    const SVec1DYMM8r4<N> &v1,
					    const SVec1DYMM8r4<N> &v2) {
			  avx256_cmean_quot<SVec1DYMM8r4<N>>(mean,v1,v2);
		  }

		  template<int32_t N>
		  void cconj_product_ymm8r4(SVec1DYMM8r4<N> &out,
					  const SVec1DYMM8r4<N> &v1,
					  const SVec1DYMM8r4<N> &v2,
					  const bool do_nt_store) {
			  avx256_cconj_prod<SVec1DYMM8r4<N>>(out,v1,v2,do_nt_store);
		  }

		  template<int32_t N>
		  void cnorm_conjprod_ymm8r4(SVec1DYMM8r4<N> &out,
					   const SVec1DYMM8r4<N> &v1,
					   const SVec1DYMM8r4<N> &v2,
					   const bool do_nt_store) {
			  avx256_cnorm_conjprod<SVec1DYMM8r4<N>>(out,v1,v2,do_nt_store);
		  }

		  template<int32_t N>
		  void cmean_conjprod_ymm8r4( SVec1DYMM8r4<N> &out,
					    const SVec1DYMM8r4<N> &x,
					    const SVec1DYMM8r4<N> &y) {
			  avx256_cmean_conjprod<SVec1DYMM8r4<N>>(out,v1,v2);
		  }

		  template<int32_t N>
		  void carith_mean_ymm8r4(std::complex<float> &mean,
					const SVec1DYMM8r4<N> &v) {
			  avx256_arith_mean<SVec1DYMM8r4<N>>(mean,v);
		  }

		  template<int32_t N>
		  void cnormalize_ymm8r4(SVec1DYMM8r4<N> &norm,
					const SVec1DYMM8r4<N> &v,
					const SVec1DYMM8r4<N> &cv,
					const bool do_nt_store) {
			  avx256_cnormalize<SVec1DYMM8r4<N>>(norm,v,cv,do_nt_store);
		  }

		  template<int32_t N>
		  void cmagnitude_ymm8r4( SVec1DYMM8r4<N> &vmag,
					 const SVec1DYMM8r4<N> &v,
					 const SVec1DYMM8r4<N> &v2) {
			  avx256_cmagnitude<SVec1DYMM8r4<N>>(vmag,v,v2);
		  }
	}
}


#endif /*__GMS_STATIC_CVEC1D_YMM8R4_H__*/
