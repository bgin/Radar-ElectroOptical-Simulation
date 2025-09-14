
#ifndef __GMS_STATIC_CVEC1D_YMM4R8_H__
#define __GMS_STATIC_CVEC1D_YMM4R8_H__ 061020191135

namespace file_info {

	
	const unsigned int gGMS_STATIC_CVEC1D_YMM4R8_MAJOR = 1U;

	const unsigned int gGMS_STATIC_CVEC1D_YMM4R8_MINOR = 0U;

	const unsigned int gGMS_STATIC_CVEC1D_YMM4R8_MICRO = 0U;

	const unsigned int gGMS_STATIC_CVEC1D_YMM4R8_FULLVER = 
		1000U*gLAM_AVXCOMPLEX_SMALLV_MAJOR+100U*gLAM_AVXCOMPLEX_SMALLV_MINOR+10U*gLAM_AVXCOMPLEX_SMALLV_MICRO;

	const char * const pgGMS_STATIC_CVEC1D_YMM4R8_CREATE_DATE = "06-10-2019 11:35 + 00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_STATIC_CVEC1D_YMM4R8_BUILD_DATE = __DATE__ ":" __TIME__;

	const char * const pgGMS_STATIC_CVEC1D_YMM4R8_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_STATIC_CVEC1D_YMM4R8_SYNOPSIS = "AVX complex vector (1D) stack-allocated storage.";
}

#include <cstdint>
#include <iostream>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_common.h"
#include "GMS_complex_common_ymm4r8.h"

#if !defined (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) // Streaming stores defined per this struct (default set to 0)
#define USE_STATIC_CVEC1D_YMM4R8_NT_STORES 0
#endif

namespace gms {
	namespace math {

#if !defined (STATIC_CVEC1D_YMM4R8_LOAD_YMM)
#define STATIC_CVEC1D_YMM4R8_LOAD_YMM(reg1,reg2,reg3,reg4,v1,v2,idx,off) \
	(reg1) = _mm256_load_pd(&(v1).m_Re[(idx)+(off)]);            \
	(reg2) = _mm256_load_pd(&(v2).m_Re[(idx)+(off)]);            \
	(reg3) = _mm256_load_pd(&(v1).m_Im[(idx)+(off)]);			  \
	(reg4) = _mm256_load_pd(&(v2).m_Im[(idx)+(off)]);
#endif

		
		template<int32_t N>
		struct __ATTR_ALIGN__(32) SVec1DYMM4r8 {

                       __ATTR_ALIGN__(32) double m_Re[(m_nsize == 0) ? 4 : N];
		       __ATTR_ALIGN__(32) double m_Im[(m_nsize == 0) ? 4 : N];
                        int32_t m_nsize = N;
			
			SVec1DYMM4r8() noexcept(true) {
				m_Re[N];
				m_Im[N];
			}

			SVec1DYMM4r8(const double Re[N],
				       const double Im[N]) {
				using namespace gms::common;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				avx256_uncached_memmove(&m_Re[0], &Re[0], N);
				avx256_uncached_memmove(&m_Im[0], &Im[0], N);
#else
				avx256_cached_memmove(&m_Re[0], &Re[0], N);
				avx256_cached_memmove(&m_Im[0], &Im[0], N);
#endif
			}

			SVec1DYMM4r8(const SVec1DYMM4r8 &x) {
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

			~SVec1DYMM4r8() = default;

			SVec1DYMM4r8 &
			operator=(const SVec1DYMM4r8 &x) {
				using namespace gms::common;
				if (this == &x || m_nsize != x.m_nsize) 
				    { return (*this); }
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
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
			   const SVec1DYMM4r8<N> &x) {
			for (int32_t i = 0; i != x.m_nsize; ++i) {
				os << std::fixed << std::showpoint << std::setprecision(15) <<
					std::setw(4)  << "Re: " << "{" << x.m_Re[i] << "}" <<
					std::setw(12) << "Im: " << "{" << x.m_Im[i] << "}" << std::endl;
			}
			 return (os);
		 }


		 template<int32_t N> SVec1DYMM4r8<N>
		 inline operator+(const SVec1DYMM4r8<N> &x,
				  const SVec1DYMM4r8<N> &y) {
			 if (x.m_nsize != y.m_nsize) 
			     { return (SVec1DYMM4r8<N>{}); }
				 SVec1DYMM4r8<N> ret_vec;
				 int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				 for (i = 0; i != ROUND_TO_FOUR(x.m_nsize, 4); i += 8) {
					 // Linearly growing indices, no need for software prefetch.
					 // HW prefetch will kick in after 2 maybe 3 cache misses.
					 const __m256d ymm0(_mm256_load_pd(&x.m_Re[i + 0]));
					 const __m256d ymm1(_mm256_load_pd(&y.m_Re[i + 0]));
					 _mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					 const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
					 const __m256d ymm3(_mm256_load_pd(&y.m_Re[i + 4]));
					 _mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));
					 const __m256d ymm4(_mm256_load_pd(&x.m_Im[i + 0]));
					 const __m256d ymm5(_mm256_load_pd(&y.m_Im[i + 0]));
					 _mm256_stream_pd(&ret_vec.m_Im[i + 0], _mm256_add_pd(ymm4, ymm5));
					 const __m256d ymm6(_mm256_load_pd(&x.m_Im[i + 4]));
					 const __m256d ymm7(_mm256_load_pd(&y.m_Im[i + 4]));
					 _mm256_stream_pd(&ret_vec.m_Im[i + 4], _mm256_add_pd(ymm6, ymm7));
				}	
				 _mm_sfence();
				 for (; i != ret_vec.m_nsize; ++i) {
					 ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
					 ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
				 }
#else					
				 for (i = 0; i != ROUND_TO_FOUR(x.m_nsize, 4); i += 8) {
					 // Linearly growing indices, no need for software prefetch.
					 // HW prefetch will kick in after 2 maybe 3 cache misses.
					 const __m256d ymm0(_mm256_load_pd(&x.m_Re[i + 0]));
					 const __m256d ymm1(_mm256_load_pd(&y.m_Re[i + 0]));
					 _mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					 const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
					 const __m256d ymm3(_mm256_load_pd(&y.m_Re[i + 4]));
					 _mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));
					 const __m256d ymm4(_mm256_load_pd(&x.m_Im[i + 0]));
					 const __m256d ymm5(_mm256_load_pd(&y.m_Im[i + 0]));
					 _mm256_store_pd(&ret_vec.m_Im[i + 0], _mm256_add_pd(ymm4, ymm5));
					 const __m256d ymm6(_mm256_load_pd(&x.m_Im[i + 4]));
					 const __m256d ymm7(_mm256_load_pd(&y.m_Im[i + 4]));
					 _mm256_store_pd(&ret_vec.m_Im[i + 4], _mm256_add_pd(ymm6, ymm7));
				 }
				 for (; i != ret_vec.m_nsize; ++i) {
					 ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
					 ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
				 }
#endif					
				 return (ret_vec)
			}		

				 
				



				
					 

					 

					

				 
				 
				 
			 

			 template<int32_t N> SVec1DYMM4r8<N>	
			 inline operator+(const SVec1DYMM4r8<N> &x, 
					  const double __restrict Re[N]) {   // If Re is not equal to x --> udefined behaviour.
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM4r8<N>{}); }
				SVec1DYMM4r8<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM4r8<N>
		        inline operator+(const double __restrict Re[N],
					 const SVec1DYMM4r8<N> &x) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM4r8<N>{}); }
				SVec1DYMM4r8<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_load_pd(&Re[i + 0]));
					const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
					_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&Re[i + 4]));
					const __m256d ymm3(_mm256_load_pd(&x.m_Re[i + 4]));
					_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = Re[i] + x.m_Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_load_pd(&Re[i + 0]));
					const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
					_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&Re[i + 4]));
					const __m256d ymm3(_mm256_load_pd(&x.m_Re[i + 4]));
					_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = Re[i] + x.m_Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM4r8<N>
			inline operator+(SVec1DYMM4r8<N> &x,
					 const double c) {
				SVec1DYMM4r8<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
					_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
					_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm0, ymm2));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + c;
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
					_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
					_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_add_pd(ymm0, ymm2));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] + c;
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM4r8<N>
			inline operator-( const SVec1DYMM4r8<N> &x,
					  const SVec1DYMM4r8<N> &y) {
				using namespace gms::common;
				if (x.m_nsize != y.m_nsize) { return (SVec1DYMM4r8<N>{}); }
				int32_t i;
				SVec1DYMM4r8<N> ret_vec;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&y.m_Re[i + 0]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&y.m_Re[i + 4]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.m_Im[i + 0]);
					const __m256d ymm5 = _mm256_load_pd(&y.m_Im[i + 0]);
					_mm256_stream_pd(&ret_vec.m_Im[i + 0], _mm256_sub_pd(ymm4, ymm5));
					const __m256d ymm6 = _mm256_load_pd(&x.m_Im[i + 4]);
					const __m256d ymm7 = _mm256_load_pd(&y.m_Im[i + 4]);
					_mm256_stream_pd(&ret_vec.m_Im[i + 4], _mm256_sub_pd(ymm6, ymm7));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&y.m_Re[i + 0]);
					_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&y.m_Re[i + 4]);
					_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.m_Im[i + 0]);
					const __m256d ymm5 = _mm256_load_pd(&y.m_Im[i + 0]);
					_mm256_store_pd(&ret_vec.m_Im[i + 0], _mm256_sub_pd(ymm4, ymm5));
					const __m256d ymm6 = _mm256_load_pd(&x.m_Im[i + 4]);
					const __m256d ymm7 = _mm256_load_pd(&y.m_Im[i + 4]);
					_mm256_store_pd(&ret_vec.m_Im[i + 4], _mm256_sub_pd(ymm6, ymm7));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM4r8<N>
			inline operator-(const SVec1DYMM4r8<N> &x,
					 const double __restrict Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM4r8<N>{}) };
				SVec1DYMM4r8<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM4r8<N>
			inline operator-(const SVec1DYMM4r8<N> &x,
					 const double c) {
				SVec1DYMM4r8<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
					_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(ymm1, ymm0));
					const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
					_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm0));


				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - c;
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
					_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(ymm1, ymm0));
					const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
					_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm0));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] - c;
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> SVec1DYMM4r8<N>
			inline operator*(const SVec1DYMM4r8<N> &x,
					 const SVec1DYMM4r8<N> &y) {
				if (x.m_nsize != y.m_nsize) { return (SVec1DYMM4r8<N>{}); }
				SVec1DYMM4r8<N> ret_vec;

				__ATTR_ALIGN__(32) struct {
                                        __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
				}ca;

				ca.ymm0 = _mm256_setzero_pd();
				ca.ymm1 = _mm256_setzero_pd();
				ca.ymm2 = _mm256_setzero_pd();
				ca.ymm3 = _mm256_setzero_pd();
				ca.ymm4 = _mm256_setzero_pd();
				ca.ymm5 = _mm256_setzero_pd();
				ca.ymm6 = _mm256_setzero_pd();
				ca.ymm7 = _mm256_setzero_pd();
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm0, ca.ymm1, ca.ymm2, ca.ymm3, x, y, i, 0)
						_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm0, ca.ymm1), _mm256_mul_pd(ca.ymm2, ca.ymm3)));
					_mm256_stream_pd(&ret_vec.m_Im[i + 0], _mm256_add_pd(
						_mm256_mul_pd(ca.ymm2, ca.ymm1), _mm256_mul_pd(ca.ymm0, ca.ymm3)));

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm4, ca.ymm5, ca.ymm6, ca.ymm7, x, y, i, 4)
						_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm4, ca.ymm5), _mm256_mul_pd(ca.ymm6, ca.ymm7)));
					_mm256_stream_pd(&ret_vec.m_Im[i + 4], _mm256_add_pd(
						_mm512_mul_pd(ca.ymm6, ca.ymm5), _mm512_mul_pd(ca.ymm4, ca.ymm7)));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
					ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm0, ca.ymm1, ca.ymm2, ca.ymm3, x, y, i, 0)
						_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm0, ca.ymm1), _mm256_mul_pd(ca.ymm2, ca.ymm3)));
					_mm256_store_pd(&ret_vec.m_Im[i + 0], _mm256_add_pd(
						_mm256_mul_pd(ca.ymm2, ca.ymm1), _mm256_mul_pd(ca.ymm0, ca.ymm3)));

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm4, ca.ymm5, ca.ymm6, ca.ymm7, x, y, i, 4)
						_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm4, ca.ymm5), _mm256_mul_pd(ca.ymm6, ca.ymm7)));
					_mm256_store_pd(&ret_vec.m_Im[i + 4], _mm256_add_pd(
						_mm512_mul_pd(ca.ymm6, ca.ymm5), _mm512_mul_pd(ca.ymm4, ca.ymm7)));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
					ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
				}
#endif
				
				return (ret_vec);
		  }

		  template<int32_t N> SVec1DYMM4r8<N>
		  inline operator*(const SVec1DYMM4r8<N> &x,
				   const double  __restrict Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM4r8<N>{}); }
				SVec1DYMM4r8<N> ret_vec;
				int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_mul_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.m_Im[i + 0]);
					_mm256_stream_pd(&ret_vec.m_Im[i + 0], _mm256_mul_pd(ymm4, ymm1));
					const __m256d ymm5 = _mm256_load_pd(&x.m_Im[i + 4]);
					_mm256_stream_pd(&ret_vec.m_Im[i + 4], _mm256_mul_pd(ymm5, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] * Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] * Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_mul_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.m_Im[i + 0]);
					_mm256_store_pd(&ret_vec.m_Im[i + 0], _mm256_mul_pd(ymm4, ymm1));
					const __m256d ymm5 = _mm256_load_pd(&x.m_Im[i + 4]);
					_mm256_store_pd(&ret_vec.m_Im[i + 4], _mm256_mul_pd(ymm5, ymm3));

				}
				for (; i != ret_vec.m_nsize; ++i) {
					ret_vec.m_Re[i] = x.m_Re[i] * Re[i];
					ret_vec.m_Im[i] = x.m_Im[i] * Re[i];
				}
#endif
				
				return (ret_vec);
		  }

		  template<int32_t N> SVec1DYMM4r8<N>
		  inline operator*(const SVec1DYMM4r8<N> &x,
				   const double c) {
			SVec1DYMM4r8<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

				const __m256d ymm0(_mm256_set1_pd(c));
				const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
				_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_mul_pd(ymm1, ymm0));
				const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
				_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm0));

			}
			_mm_sfence();
			for (; i != ret_vec.m_nsize; ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] * c;
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {

				const __m256d ymm0(_mm256_set1_pd(c));
				const __m256d ymm1(_mm256_load_pd(&x.m_Re[i + 0]));
				_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_mul_pd(ymm1, ymm0));
				const __m256d ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
				_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm0));

			}
			for (; i != ret_vec.m_nsize; ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] * c;
			}
#endif
			
			return (ret_vec);
		}

		template<int32_t N> SVec1DYMM4r8<N>
		inline operator/(const SVec1DYMM4r8<N> &x,
				 const SVec1DYMM4r8<N> &y) {
			if (x.m_nsize != y.m_nsize) { return (SVec1DYMM4r8<N>{}); }
			SVec1DYMM4r8<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.

				const __m256d ymm0 ( _mm256_load_pd(&x.m_Re[i + 0]));
				const __m256d ymm1 ( _mm256_load_pd(&y.m_Im[i + 0]));
				const __m256d ymm2 ( _mm256_load_pd(&x.m_Im[i + 0]));
				const __m256d re_term1 ( _mm256_add_pd(
					_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm1)));
				const __m256d re_term2 ( _mm256_add_pd(
					_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0, ymm1)));
				const __m256d ymm3 ( _mm256_load_pd(&y.m_Re[i + 0]));
				const __m256d den_term ( _mm256_add_pd(
					_mm256_mul_pd(ymm3, ymm3), _mm256_mul_pd(ymm1, ymm1)));

				_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_div_pd(re_term1, den_term));
				_mm256_stream_pd(&ret_vec.m_Im[i + 0], _mm256_div_pd(re_term2, den_term));
			}
			_mm_sfence();
			for (; i != ret_vec.m_nsize; ++i) {
				const double tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
				const double tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
				const double den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
				ret_vec.m_Re[i] = tre / den;
				ret_vec.m_Im[i] = tim / den;
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.

				const __m256d ymm0(_mm256_load_pd(&x.m_Re[i + 0]));
				const __m256d ymm1(_mm256_load_pd(&y.m_Im[i + 0]));
				const __m256d ymm2(_mm256_load_pd(&x.m_Im[i + 0]));
				const __m256d re_term1(_mm256_add_pd(
					_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm1)));
				const __m256d re_term2(_mm256_add_pd(
					_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0, ymm1)));
				const __m256d ymm3(_mm256_load_pd(&y.m_Re[i + 0]));
				const __m256d den_term(_mm256_add_pd(
					_mm256_mul_pd(ymm3, ymm3), _mm256_mul_pd(ymm1, ymm1)));

				_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_div_pd(re_term1, den_term));
				_mm256_store_pd(&ret_vec.m_Im[i + 0], _mm256_div_pd(re_term2, den_term));
			}
			
			for (; i != ret_vec.m_nsize; ++i) {
				const double tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
				const double tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
				const double den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
				ret_vec.m_Re[i] = tre / den;
				ret_vec.m_Im[i] = tim / den;
			}
#endif		
			return (ret_vec);
		 }

		 template<int32_t N> SVec1DYMM4r8<N>
		 inline operator/(const SVec1DYMM4r8<N> &x,
				  const double __restrict Re[N]) {
			using namespace gms::common;
			if (!Is_ptr_aligned32(Re)) { return (SVec1DYMM4r8<N>{}); }
			SVec1DYMM4r8<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.
				const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i]);
				const __m256d ymm1 = _mm256_load_pd(&Re[i]);
				_mm256_stream_pd(&ret_vec.m_Re[i], _mm256_div_pd(ymm0, ymm1));
				const __m256d ymm2 = _mm256_load_pd(&x.m_Im[i]);
				_mm256_stream_pd(&ret_vec.m_Im[i], _mm256_div_pd(ymm2, ymm1));
			}
			_mm_sfence();
			for (; i != x.size(); ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] / Re[i];
				ret_vec.m_Im[i] = x.m_Im[i] / Re[i];
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.
				const __m256d ymm0 = _mm256_load_pd(&x.m_Re[i]);
				const __m256d ymm1 = _mm256_load_pd(&Re[i]);
				_mm256_store_pd(&ret_vec.m_Re[i], _mm256_div_pd(ymm0, ymm1));
				const __m256d ymm2 = _mm256_load_pd(&x.m_Im[i]);
				_mm256_store_pd(&ret_vec.m_Im[i], _mm256_div_pd(ymm2, ymm1));
			}
			for (; i != x.size(); ++i) {
				ret_vec.m_Re[i] = x.m_Re[i] / Re[i];
				ret_vec.m_Im[i] = x.m_Im[i] / Re[i];
			}
#endif
			return (ret_vec);
		}

		template<int32_t N> SVec1DYMM4r8<N>
		inline operator==(const SVec1DYMM4r8<N> &x,
				  const SVec1DYMM4r8<N> &y) {
			using namespace gms::common;
			if (x.m_nsize != y.m_nsize) { return (SVec1DYMM4r8<N>{}); }

			__ATTR_ALIGN__(32) struct {
                                __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
			}ca;

			ca.ymm0 = _mm256_setzero_pd(); 
			ca.ymm1 = _mm256_setzero_pd();
			ca.ymm2 = _mm256_setzero_pd();
			ca.ymm3 = _mm256_setzero_pd();
			ca.ymm4 = _mm256_setzero_pd();
			ca.ymm5 = _mm256_setzero_pd();
			ca.ymm6 = _mm256_setzero_pd();
			ca.ymm7 = _mm256_setzero_pd();
			SVec1DYMM4r8<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
			for(i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize,4); i += 8) {
				 ca.ymm0(_mm256_load_pd(&x.m_Re[i+0]));
				 ca.ymm1(_mm256_load_pd(&y.m_Re[i+0]));
				 _mm256_stream_pd(&ret_vec.m_Re[i+0], _mm256_cmp_pd(ca.ymm0,ca.ymm1,_CMP_EQ_OQ));
				 ca.ymm2(_mm256_load_pd(&x.m_Re[i+4]));
				 ca.ymm3(_mm256_load_pd(&y.m_Re[i+4]));
				 _mm256_stream_pd(&ret_vec.m_Re[i+4], _mm256_cmp_pd(ca.ymm2,ca.ymm3,_CMP_EQ_OQ));
				 ca.ymm4(_mm256_load_pd(&x.m_Im[i+0]));
				 ca.ymm5(_mm256_load_pd(&y.m_Im[i+0]));
				 _mm256_stream_pd(&ret_vec.m_Im[i+0], _mm256_cmp_pd(ca.ymm4,ca.ymm5,_CMP_EQ_OQ));
				 ca.ymm6(_mm256_load_pd(&x.m_Im[i+4]));
				 ca.ymm7(_mm256_load_pd(&y.m_Im[i+1]));
				 _mm256_stream_pd(&ret_vec.m_Im[i+4],_mm256_cmp_pd(ca.ymm6,ca.ymm7,_CMP_EQ_OQ));
			}
			_mm_sfence();
			for(; i != ret_vec.m_nsize; ++i) {
				if(approximately_equalf64(x.m_Re[i],
					        y.m_Re[i],std::numeric_limits<double>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				 }
				 else {
					 ret_vec.m_Re[i] = 0.0;
				 }
				 if (approximately_equalf64(x.m_Im[i],
						y.m_Im[i],std::numeric_limits<double>::epsilon())) {
					 ret_vec.m_Im[i] = 1.0;
				}
				else {
					ret_vec.m_Im[i] = 0.0;
				}
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {
				ca.ymm0(_mm256_load_pd(&x.m_Re[i + 0]));
				ca.ymm1(_mm256_load_pd(&y.m_Re[i + 0]));
				_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_cmp_pd(ca.ymm0, ca.ymm1, _CMP_EQ_OQ));
				ca.ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
				ca.ymm3(_mm256_load_pd(&y.m_Re[i + 4]));
				_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_cmp_pd(ca.ymm2, ca.ymm3, _CMP_EQ_OQ));
				ca.ymm4(_mm256_load_pd(&x.m_Im[i + 0]));
				ca.ymm5(_mm256_load_pd(&y.m_Im[i + 0]));
				_mm256_store_pd(&ret_vec.m_Im[i + 0], _mm256_cmp_pd(ca.ymm4, ca.ymm5, _CMP_EQ_OQ));
				ca.ymm6(_mm256_load_pd(&x.m_Im[i + 4]));
				ca.ymm7(_mm256_load_pd(&y.m_Im[i + 4]));
				_mm256_store_pd(&ret_vec.m_Im[i + 4], _mm256_cmp_pd(ca.ymm6, ca.ymm7, _CMP_EQ_OQ));
			}
			
			for (; i != ret_vec.m_nsize; ++i) {
				if (approximately_equalf64(x.m_Re[i],
					y.m_Re[i], std::numeric_limits<double>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				}
				else {
					ret_vec.m_Re[i] = 0.0;
				}
				if (approximately_equalf64(x.m_Im[i],
					y.m_Im[i], std::numeric_limits<double>::epsilon())) {
					ret_vec.m_Im[i] = 1.0;
				}
				else {
					ret_vec.m_Im[i] = 0.0;
				}
		  }
#endif
				return (ret_vec);
		}

		template<int32_t N> SVec1DYMM4r8<N>
		inline operator!=(const SVec1DYMM4r8<N> &x,
				  const SVec1DYMM4r8<N> &y) {
			using namespace gms::common;
		        if (x.m_nsize != y.m_nsize) { return (SVec1DYMM4r8<N>{}); }

			__ATTR_ALIGN__(32) struct {
                                __m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
			}ca;

			ca.ymm0 = _mm256_setzero_pd();
			ca.ymm1 = _mm256_setzero_pd();
			ca.ymm2 = _mm256_setzero_pd();
			ca.ymm3 = _mm256_setzero_pd();
			ca.ymm4 = _mm256_setzero_pd();
			ca.ymm5 = _mm256_setzero_pd();
			ca.ymm6 = _mm256_setzero_pd();
			ca.ymm7 = _mm256_setzero_pd();
			SVec1DYMM4r8<N> ret_vec;
			int32_t i;
#if (USE_STATIC_CVEC1D_YMM4R8_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {
				ca.ymm0(_mm256_load_pd(&x.m_Re[i + 0]));
				ca.ymm1(_mm256_load_pd(&y.m_Re[i + 0]));
				_mm256_stream_pd(&ret_vec.m_Re[i + 0], _mm256_cmp_pd(ca.ymm0, ca.ymm1, _CMP_NEQ_OQ));
				ca.ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
				ca.ymm3(_mm256_load_pd(&y.m_Re[i + 4]));
				_mm256_stream_pd(&ret_vec.m_Re[i + 4], _mm256_cmp_pd(ca.ymm2, ca.ymm3, _CMP_NEQ_OQ));
				ca.ymm4(_mm256_load_pd(&x.m_Im[i + 0]));
				ca.ymm5(_mm256_load_pd(&y.m_Im[i + 0]));
				_mm256_stream_pd(&ret_vec.m_Im[i + 0], _mm256_cmp_pd(ca.ymm4, ca.ymm5, _CMP_NEQ_OQ));
				ca.ymm6(_mm256_load_pd(&x.m_Im[i + 4]));
				ca.ymm7(_mm256_load_pd(&y.m_Im[i + 1]));
				_mm256_stream_pd(&ret_vec.m_Im[i + 4], _mm256_cmp_pd(ca.ymm6, ca.ymm7, _CMP_NEQ_OQ));
			}
			_mm_sfence();
			for (; i != ret_vec.m_nsize; ++i) {
				if (!approximately_equalf64(x.m_Re[i],
					y.m_Re[i], std::numeric_limits<double>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				}
				else {
					ret_vec.m_Re[i] = 0.0;
				}
				if (!approximately_equalf64(x.m_Im[i],
					y.m_Im[i], std::numeric_limits<double>::epsilon())) {
					ret_vec.m_Im[i] = 1.0;
				}
				else {
					ret_vec.m_Im[i] = 0.0;
				}
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.m_nsize, 4); i += 8) {
				ca.ymm0(_mm256_load_pd(&x.m_Re[i + 0]));
				ca.ymm1(_mm256_load_pd(&y.m_Re[i + 0]));
				_mm256_store_pd(&ret_vec.m_Re[i + 0], _mm256_cmp_pd(ca.ymm0, ca.ymm1, _CMP_NEQ_OQ));
				ca.ymm2(_mm256_load_pd(&x.m_Re[i + 4]));
				ca.ymm3(_mm256_load_pd(&y.m_Re[i + 4]));
				_mm256_store_pd(&ret_vec.m_Re[i + 4], _mm256_cmp_pd(ca.ymm2, ca.ymm3, _CMP_NEQ_OQ));
				ca.ymm4(_mm256_load_pd(&x.m_Im[i + 0]));
				ca.ymm5(_mm256_load_pd(&y.m_Im[i + 0]));
				_mm256_store_pd(&ret_vec.m_Im[i + 0], _mm256_cmp_pd(ca.ymm4, ca.ymm5, _CMP_NEQ_OQ));
				ca.ymm6(_mm256_load_pd(&x.m_Im[i + 4]));
				ca.ymm7(_mm256_load_pd(&y.m_Im[i + 4]));
				_mm256_store_pd(&ret_vec.m_Im[i + 4], _mm256_cmp_pd(ca.ymm6, ca.ymm7, _CMP_NEQ_OQ));
			}

			for (; i != ret_vec.m_nsize; ++i) {
				if (!approximately_equalf64(x.m_Re[i],
					y.m_Re[i], std::numeric_limits<double>::epsilon())) {
					ret_vec.m_Re[i] = 1.0;
				}
				else {
					ret_vec.m_Re[i] = 0.0;
				}
				if (!approximately_equalf64(x.m_Im[i],
					y.m_Im[i], std::numeric_limits<double>::epsilon())) {
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
		  void cnormalize_product_ymm4r8(SVec1DYMM4r8<N> &out,
					       const SVec1DYMM4r8<N> &v1,
					       const SVec1DYMM4r8<N> &v2) {
			  avx256_cnormalize_prod<SVec1DYMM4r8<N>>(out,v1,v2);
		  }

		  template<int32_t N>
		  void cmean_product_ymm4r8(std::complex<double> &mean,
					  const SVec1DYMM4r8<N> &v1,
					  const SVec1DYMM4r8<N> &v2) {
			  avx256_cmean_prod<SVec1DYMM4r8<N>>(mean,v1,v2);
		  }

		  template<int32_t N>
		  void cmean_quotient_ymm4r8( std::complex<double> &mean,
					    const SVec1DYMM4r8<N> &v1,
					    const SVec1DYMM4r8<N> &v2) {
			  avx256_cmean_quot<SVec1DYMM4r8<N>>(mean,v1,v2);
		  }

		  template<int32_t N>
		  void cconj_product_ymm4r8(SVec1DYMM4r8<N> &out,
					  const SVec1DYMM4r8<N> &v1,
					  const SVec1DYMM4r8<N> &v2,
					  const bool do_nt_store) {
			  avx256_cconj_prod<SVec1DYMM4r8<N>>(out,v1,v2,do_nt_store);
		  }

		  template<int32_t N>
		  void cnorm_conjprod_ymm4r8(SVec1DYMM4r8<N> &out,
					   const SVec1DYMM4r8<N> &v1,
					   const SVec1DYMM4r8<N> &v2,
					   const bool do_nt_store) {
			  avx256_cnorm_conjprod<SVec1DYMM4r8<N>>(out,v1,v2,do_nt_store);
		  }

		  template<int32_t N>
		  void cmean_conjprod_ymm4r8( SVec1DYMM4r8<N> &out,
					    const SVec1DYMM4r8<N> &x,
					    const SVec1DYMM4r8<N> &y) {
			  avx256_cmean_conjprod<SVec1DYMM4r8<N>>(out,v1,v2);
		  }

		  template<int32_t N>
		  void carith_mean_ymm4r8(std::complex<double> &mean,
					const SVec1DYMM4r8<N> &v) {
			  avx256_arith_mean<SVec1DYMM4r8<N>>(mean,v);
		  }

		  template<int32_t N>
		  void cnormalize_ymm4r8(SVec1DYMM4r8<N> &norm,
					const SVec1DYMM4r8<N> &v,
					const SVec1DYMM4r8<N> &cv,
					const bool do_nt_store) {
			  avx256_cnormalize<SVec1DYMM4r8<N>>(norm,v,cv,do_nt_store);
		  }

		  template<int32_t N>
		  void cmagnitude_ymm4r8( SVec1DYMM4r8<N> &vmag,
					 const SVec1DYMM4r8<N> &v,
					 const SVec1DYMM4r8<N> &v2) {
			  avx256_cmagnitude<SVec1DYMM4r8<N>>(vmag,v,v2);
		  }
	}
}


#endif /*__GMS_STATIC_CVEC1D_YMM4R8_H__*/
