
#ifndef __GMS_AVXCOMPLEX_SHORTVEC_H__
#define __GMS_AVXCOMPLEX_SHORTVEC_H__ 061020191135

namespace file_info {

	
	const unsigned int gGMS_AVXCOMPLEX_SHORTVEC_MAJOR = 1;
	const unsigned int gGMS_AVXCOMPLEX_SHORTVEC_MINOR = 0;
	const unsigned int gGMS_AVXCOMPLEX_SHORTVEC_MICRO = 1;

	const unsigned int gGMS_AVXCOMPLEX_SHORTVEC_FULLVER = 
		1000U*gGMS_AVXCOMPLEX_SHORTVECV_MAJOR+100U*gGMS_AVXCOMPLEX_SHORTVEC_MINOR+10U*gGMS_AVXCOMPLEX_SHORTVEC_MICRO;

	const char * const pgGMS_AVXCOMPLEX_SHORTVEC_CREATE_DATE = "06-10-2019 11:35 + 00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_AVXCOMPLEX_SHORTVEC_BUILD_DATE = __DATE__":"__TIME__;

	const char * const pgGMS_AVXCOMPLEX_SHORTVEC_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVXCOMPLEX_SHORTVEC_SYNOPSIS = "AVX complex vector (1D) stack-allocated storage.";
}

#include <cstdint>
#include <iostream>
#include "GMS_config.h"
#include "GMS_common.h"
#include "GMS_avxcomplex_common.h"
#include "GMS_constants.h"


#if !defined (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) // Streaming stores defined per this struct (default set to 0)
#define USE_AVXCOMPLEX_SHORTVEC_NT_STORES 0
#endif

namespace gms {
	namespace math {

#if !defined (AVXCOMPLEX_SHORTVEC_LOAD_YMM)
#define AVXCOMPLEX_SHORTVEC_LOAD_YMM(reg1,reg2,reg3,reg4,v1,v2,idx,off) \
	(reg1) = _mm256_load_pd(&(v1).data.m_Re[(idx)+(off)]);            \
	(reg2) = _mm256_load_pd(&(v2).data.m_Re[(idx)+(off)]);            \
	(reg3) = _mm256_load_pd(&(v1).data.m_Im[(idx)+(off)]);			  \
	(reg4) = _mm256_load_pd(&(v2).data.m_Im[(idx)+(off)]);
#endif

		template<int32_t N>
		struct AVXSCVData {

			static constexpr int32_t MAX_SIZE = 4096;
			int32_t m_nsize = N;
			PAD_TO(1,4)
		        PAD_TO(2,4)
			static_assert(N <= MAX_SIZE, "Invalid size of AVXSCVData -- has been passed!!");


		        __attribute__((align(64))) double m_Re[(N == 0) ? 4 : N];

		        __attribute__((align(64))) double m_Im[(N == 0) ? 4 : N];

		};

		template<int32_t N>
		struct __attribute__((align(64))) AVXShortCVec1D {


			__attribute__((align(64))) AVXSCVData<N> data;

			AVXShortCVec1D() noexcept(true) {
				data.m_Re[N] = {};
				data.m_Im[N] = {};
			}

			AVXShortCVec1D(const double Re[N],
				       const double Im[N]) {
				using namespace gms::common;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				avx256_uncached_memmove(&data.m_Re[0], &Re[0], N);
				avx256_uncached_memmove(&data.m_Im[0], &Im[0], N);
#else
				avx256_cached_memmove(&data.m_Re[0], &Re[0], N);
				avx256_cached_memmove(&data.m_Im[0], &Im[0], N);
#endif
			}

			AVXShortCVec1D(const AVXSmallCVec1D &x) {
				using namespace gms::common;
				data.m_nsize = x.data.m_nsize;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				avx256_uncached_memmove(&data.m_Re[0], &x.data.m_Re[0], x.data.m_nsize);
				avx256_uncached_memmove(&data.m_Im[0], &x.data.m_Im[0], x.data.m_nsize);
#else
				avx256_cached_memmove(&data.m_Re[0], &x.data.m_Re[0], x.data.m_nsize);
				avx256_cached_memmove(&data.m_Im[0], &x.data.m_Im[0], x.data.m_nsize);
#endif
			}

			~AVXShortCVec1D() = default;

			AVXShortCVec1D &
			operator=(const AVXShortCVec1D &x) {
				using namespace gms::common;
				if (this == &x || data.m_nsize != x.data.m_nsize) 
				    { return (*this); }
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				avx256_uncached_memmove(&data.m_Re[0], &x.data.m_Re[0],x.data.m_nsize);
				avx256_uncached_memmove(&data.m_Im[0], &x.data.m_Im[0],x.data.m_nsize);
#else
				avx256_cached_memmove(&data.m_Re[0], &x.data.m_Re[0], x.data.m_nsize);
				avx256_cached_memmove(&data.m_Im[0], &x.data.m_Im[0], x.data.m_nsize);
#endif				
				return (*this);
				}
		};

		template<int32_t N> std::ostream &
		operator<<(std::ostream &os,
			   const AVXShortCVec1D<N> &x) {
			for (int32_t i = 0; i != x.data.m_nsize; ++i) {
				os << std::fixed << std::showpoint << std::setprecision(15) <<
					std::setw(4)  << "Re: " << "{" << x.data.m_Re[i] << "}" <<
					std::setw(12) << "Im: " << "{" << x.data.m_Im[i] << "}" << std::endl;
			}
			 return (os);
		 }


		 template<int32_t N> AVXShortCVec1D<N>
		 inline operator+(const AVXShortCVec1D<N> &x,
				  const AVXShortCVec1D<N> &y) {
			 if (x.data.m_nsize != y.data.m_nsize) 
			     { return (AVXShortCVec1D<N>{}); }
				 AVXShortCVec1D<N> ret_vec;
				 int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				 for (i = 0; i != ROUND_TO_FOUR(x.data.m_nsize, 4); i += 8) {
					 // Linearly growing indices, no need for software prefetch.
					 // HW prefetch will kick in after 2 maybe 3 cache misses.
					 const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i + 0]));
					 const __m256d ymm1(_mm256_load_pd(&y.data.m_Re[i + 0]));
					 _mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					 const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
					 const __m256d ymm3(_mm256_load_pd(&y.data.m_Re[i + 4]));
					 _mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));
					 const __m256d ymm4(_mm256_load_pd(&x.data.m_Im[i + 0]));
					 const __m256d ymm5(_mm256_load_pd(&y.data.m_Im[i + 0]));
					 _mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_add_pd(ymm4, ymm5));
					 const __m256d ymm6(_mm256_load_pd(&x.data.m_Im[i + 4]));
					 const __m256d ymm7(_mm256_load_pd(&y.data.m_Im[i + 4]));
					 _mm256_stream_pd(&ret_vec.data.m_Im[i + 4], _mm256_add_pd(ymm6, ymm7));
				}	
				 _mm_sfence();
				 for (; i != ret_vec.data.m_nsize; ++i) {
					 ret_vec.data.m_Re[i] = x.data.m_Re[i] + y.data.m_Re[i];
					 ret_vec.data.m_Im[i] = x.data.m_Im[i] + y.data.m_Re[i];
				 }
#else					
				 for (i = 0; i != ROUND_TO_FOUR(x.data.m_nsize, 4); i += 8) {
					 // Linearly growing indices, no need for software prefetch.
					 // HW prefetch will kick in after 2 maybe 3 cache misses.
					 const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i + 0]));
					 const __m256d ymm1(_mm256_load_pd(&y.data.m_Re[i + 0]));
					 _mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					 const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
					 const __m256d ymm3(_mm256_load_pd(&y.data.m_Re[i + 4]));
					 _mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));
					 const __m256d ymm4(_mm256_load_pd(&x.data.m_Im[i + 0]));
					 const __m256d ymm5(_mm256_load_pd(&y.data.m_Im[i + 0]));
					 _mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_add_pd(ymm4, ymm5));
					 const __m256d ymm6(_mm256_load_pd(&x.data.m_Im[i + 4]));
					 const __m256d ymm7(_mm256_load_pd(&y.data.m_Im[i + 4]));
					 _mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_add_pd(ymm6, ymm7));
				 }
				 for (; i != ret_vec.data.m_nsize; ++i) {
					 ret_vec.data.m_Re[i] = x.data.m_Re[i] + y.data.m_Re[i];
					 ret_vec.data.m_Im[i] = x.data.m_Im[i] + y.data.m_Re[i];
				 }
#endif					
				 return (ret_vec)
			}		

				 
				



				
					 

					 

					

				 
				 
				 
			 

			 template<int32_t N> AVXShortCVec1D<N>	
			 inline operator+(const AVXShortCVec1D<N> &x, 
					  const double __restrict Re[N]) {   // If Re is not equal to x --> udefined behaviour.
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (AVXShortCVec1D<N>{}); }
				AVXShortCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] + Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] + Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> AVXShortCVec1D<N>
		        inline operator+(const double __restrict Re[N],
					 const AVXShortCVec1D<N> &x) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (AVXShortCVec1D<N>{}); }
				AVXShortCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_load_pd(&Re[i + 0]));
					const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&Re[i + 4]));
					const __m256d ymm3(_mm256_load_pd(&x.data.m_Re[i + 4]));
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = Re[i] + x.data.m_Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_load_pd(&Re[i + 0]));
					const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
					_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&Re[i + 4]));
					const __m256d ymm3(_mm256_load_pd(&x.data.m_Re[i + 4]));
					_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = Re[i] + x.data.m_Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> AVXShortCVec1D<N>
			inline operator+(AVXShortCVec1D<N> &x,
					 const double c) {
				AVXShortCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm0, ymm2));

				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] + c;
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
					_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
					const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
					_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm0, ymm2));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] + c;
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> AVXShortCVec1D<N>
			inline operator-( const AVXShortCVec1D<N> &x,
					  const AVXSmallCVec1D<N> &y) {
				using namespace gms::common;
				if (x.data.m_nsize != y.data.m_nsize) { return (AVXSmallCVec1D<N>{}); }
				int32_t i;
				AVXShortCVec1D<N> ret_vec;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&y.data.m_Re[i + 0]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&y.data.m_Re[i + 4]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.data.m_Im[i + 0]);
					const __m256d ymm5 = _mm256_load_pd(&y.data.m_Im[i + 0]);
					_mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_sub_pd(ymm4, ymm5));
					const __m256d ymm6 = _mm256_load_pd(&x.data.m_Im[i + 4]);
					const __m256d ymm7 = _mm256_load_pd(&y.data.m_Im[i + 4]);
					_mm256_stream_pd(&ret_vec.data.m_Im[i + 4], _mm256_sub_pd(ymm6, ymm7));

				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] - y.data.m_Re[i];
					ret_vec.data.m_Im[i] = x.data.m_Im[i] - y.data.m_Im[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&y.data.m_Re[i + 0]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&y.data.m_Re[i + 4]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.data.m_Im[i + 0]);
					const __m256d ymm5 = _mm256_load_pd(&y.data.m_Im[i + 0]);
					_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_sub_pd(ymm4, ymm5));
					const __m256d ymm6 = _mm256_load_pd(&x.data.m_Im[i + 4]);
					const __m256d ymm7 = _mm256_load_pd(&y.data.m_Im[i + 4]);
					_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_sub_pd(ymm6, ymm7));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] - y.data.m_Re[i];
					ret_vec.data.m_Im[i] = x.data.m_Im[i] - y.data.m_Im[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> AVXShortCVec1D<N>
			inline operator-(const AVXShortCVec1D<N> &x,
					 const double __restrict Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (AVXShortCVec1D<N>{}) };
				AVXShortCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] - Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] - Re[i];
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> AVXShortCVec1D<N>
			inline operator-(const AVXShortCVec1D<N> &x,
					 const double c) {
				AVXShortCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm1, ymm0));
					const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm0));


				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] - c;
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0(_mm256_set1_pd(c));
					const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
					_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm1, ymm0));
					const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
					_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm0));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] - c;
				}
#endif
				
				return (ret_vec);
			}

			template<int32_t N> AVXShortCVec1D<N>
			inline operator*(const AVXShortCVec1D<N> &x,
					 const AVXShortCVec1D<N> &y) {
				if (x.data.m_nsize != y.data.m_nsize) { return (AVXShortCVec1D<N>{}); }
				AVXShortCVec1D<N> ret_vec;


				__attribute__((align(64))) struct {
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
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm0, ca.ymm1, ca.ymm2, ca.ymm3, x, y, i, 0)
						_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm0, ca.ymm1), _mm256_mul_pd(ca.ymm2, ca.ymm3)));
					_mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_add_pd(
						_mm256_mul_pd(ca.ymm2, ca.ymm1), _mm256_mul_pd(ca.ymm0, ca.ymm3)));

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm4, ca.ymm5, ca.ymm6, ca.ymm7, x, y, i, 4)
						_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm4, ca.ymm5), _mm256_mul_pd(ca.ymm6, ca.ymm7)));
					_mm256_stream_pd(&ret_vec.data.m_Im[i + 4], _mm256_add_pd(
						_mm512_mul_pd(ca.ymm6, ca.ymm5), _mm512_mul_pd(ca.ymm4, ca.ymm7)));

				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = (x.data.m_Re[i] * y.data.m_Re[i]) - (x.data.m_Im[i] * y.data.m_Im[i]);
					ret_vec.data.m_Im[i] = (x.data.m_Im[i] * y.data.m_Re[i]) + (x.data.m_Re[i] * y.data.m_Im[i]);
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm0, ca.ymm1, ca.ymm2, ca.ymm3, x, y, i, 0)
						_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm0, ca.ymm1), _mm256_mul_pd(ca.ymm2, ca.ymm3)));
					_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_add_pd(
						_mm256_mul_pd(ca.ymm2, ca.ymm1), _mm256_mul_pd(ca.ymm0, ca.ymm3)));

					AVXCOMPLEX_SMALLV_LOAD_YMM(ca.ymm4, ca.ymm5, ca.ymm6, ca.ymm7, x, y, i, 4)
						_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(
						_mm256_mul_pd(ca.ymm4, ca.ymm5), _mm256_mul_pd(ca.ymm6, ca.ymm7)));
					_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_add_pd(
						_mm512_mul_pd(ca.ymm6, ca.ymm5), _mm512_mul_pd(ca.ymm4, ca.ymm7)));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = (x.data.m_Re[i] * y.data.m_Re[i]) - (x.data.m_Im[i] * y.data.m_Im[i]);
					ret_vec.data.m_Im[i] = (x.data.m_Im[i] * y.data.m_Re[i]) + (x.data.m_Re[i] * y.data.m_Im[i]);
				}
#endif
				
				return (ret_vec);
		  }

		  template<int32_t N> AVXShortCVec1D<N>
		  inline operator*(const AVXShortCVec1D<N> &x,
				   const double  __restrict Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned32(Re)) { return (AVXShortCVec1D<N>{}); }
				AVXShortCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_mul_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.data.m_Im[i + 0]);
					_mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_mul_pd(ymm4, ymm1));
					const __m256d ymm5 = _mm256_load_pd(&x.data.m_Im[i + 4]);
					_mm256_stream_pd(&ret_vec.data.m_Im[i + 4], _mm256_mul_pd(ymm5, ymm3));

				}
				_mm_sfence();
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] * Re[i];
					ret_vec.data.m_Im[i] = x.data.m_Im[i] * Re[i];
				}
#else
				for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

					const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i + 0]);
					const __m256d ymm1 = _mm256_load_pd(&Re[i + 0]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_mul_pd(ymm0, ymm1));
					const __m256d ymm2 = _mm256_load_pd(&x.data.m_Re[i + 4]);
					const __m256d ymm3 = _mm256_load_pd(&Re[i + 4]);
					_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm3));
					const __m256d ymm4 = _mm256_load_pd(&x.data.m_Im[i + 0]);
					_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_mul_pd(ymm4, ymm1));
					const __m256d ymm5 = _mm256_load_pd(&x.data.m_Im[i + 4]);
					_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_mul_pd(ymm5, ymm3));

				}
				for (; i != ret_vec.data.m_nsize; ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] * Re[i];
					ret_vec.data.m_Im[i] = x.data.m_Im[i] * Re[i];
				}
#endif
				
				return (ret_vec);
		  }

		  template<int32_t N> AVXShortCVec1D<N>
		  inline operator*(const AVXShortCVec1D<N> &x,
				   const double c) {
			AVXSmallCVec1D<N> ret_vec;
			int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

				const __m256d ymm0(_mm256_set1_pd(c));
				const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
				_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_mul_pd(ymm1, ymm0));
				const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
				_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm0));

			}
			_mm_sfence();
			for (; i != ret_vec.data.m_nsize; ++i) {
				ret_vec.data.m_Re[i] = x.data.m_Re[i] * c;
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {

				const __m256d ymm0(_mm256_set1_pd(c));
				const __m256d ymm1(_mm256_load_pd(&x.data.m_Re[i + 0]));
				_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_mul_pd(ymm1, ymm0));
				const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
				_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm0));

			}
			for (; i != ret_vec.data.m_nsize; ++i) {
				ret_vec.data.m_Re[i] = x.data.m_Re[i] * c;
			}
#endif
			
			return (ret_vec);
		}

		template<int32_t N> AVXShortCVec1D<N>
		inline operator/(const AVXShortCVec1D<N> &x,
				 const AVXShortCVec1D<N> &y) {
			if (x.data.m_nsize != y.data.m_nsize) { return (AVXShortCVec1D<N>{}); }
			AVXShortCVec1D<N> ret_vec;
			int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.

				const __m256d ymm0 ( _mm256_load_pd(&x.data.m_Re[i + 0]));
				const __m256d ymm1 ( _mm256_load_pd(&y.data.m_Im[i + 0]));
				const __m256d ymm2 ( _mm256_load_pd(&x.data.m_Im[i + 0]));
				const __m256d re_term1 ( _mm256_add_pd(
					_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm1)));
				const __m256d re_term2 ( _mm256_add_pd(
					_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0, ymm1)));
				const __m256d ymm3 ( _mm256_load_pd(&y.data.m_Re[i + 0]));
				const __m256d den_term ( _mm256_add_pd(
					_mm256_mul_pd(ymm3, ymm3), _mm256_mul_pd(ymm1, ymm1)));

				_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_div_pd(re_term1, den_term));
				_mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_div_pd(re_term2, den_term));
			}
			_mm_sfence();
			for (; i != ret_vec.data.m_nsize; ++i) {
				const double tre = (x.data.m_Re[i] * y.data.m_Im[i]) + (x.data.m_Im[i] * y.data.m_Im[i]);
				const double tim = (x.data.m_Im[i] * y.data.m_Im[i]) - (x.data.m_Re[i] * y.data.m_Im[i]);
				const double den = (y.data.m_Re[i] * y.data.m_Re[i]) + (y.data.m_Im[i] * y.data.m_Im[i]);
				ret_vec.data.m_Re[i] = tre / den;
				ret_vec.data.m_Im[i] = tim / den;
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.

				const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i + 0]));
				const __m256d ymm1(_mm256_load_pd(&y.data.m_Im[i + 0]));
				const __m256d ymm2(_mm256_load_pd(&x.data.m_Im[i + 0]));
				const __m256d re_term1(_mm256_add_pd(
					_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm1)));
				const __m256d re_term2(_mm256_add_pd(
					_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0, ymm1)));
				const __m256d ymm3(_mm256_load_pd(&y.data.m_Re[i + 0]));
				const __m256d den_term(_mm256_add_pd(
					_mm256_mul_pd(ymm3, ymm3), _mm256_mul_pd(ymm1, ymm1)));

				_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_div_pd(re_term1, den_term));
				_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_div_pd(re_term2, den_term));
			}
			
			for (; i != ret_vec.data.m_nsize; ++i) {
				const double tre = (x.data.m_Re[i] * y.data.m_Im[i]) + (x.data.m_Im[i] * y.data.m_Im[i]);
				const double tim = (x.data.m_Im[i] * y.data.m_Im[i]) - (x.data.m_Re[i] * y.data.m_Im[i]);
				const double den = (y.data.m_Re[i] * y.data.m_Re[i]) + (y.data.m_Im[i] * y.data.m_Im[i]);
				ret_vec.data.m_Re[i] = tre / den;
				ret_vec.data.m_Im[i] = tim / den;
			}
#endif		
			return (ret_vec);
		 }

		 template<int32_t N> AVXShortCVec1D<N>
		 inline operator/(const AVXShortCVec1D<N> &x,
				  const double __restrict Re[N]) {
			using namespace gms::common;
			if (!Is_ptr_aligned32(Re)) { return (AVXShortCVec1D<N>{}); }
			AVXShortCVec1D<N> ret_vec;
			int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.
				const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i]);
				const __m256d ymm1 = _mm256_load_pd(&Re[i]);
				_mm256_stream_pd(&ret_vec.data.m_Re[i], _mm256_div_pd(ymm0, ymm1));
				const __m256d ymm2 = _mm256_load_pd(&x.data.m_Im[i]);
				_mm256_stream_pd(&ret_vec.data.m_Im[i], _mm256_div_pd(ymm2, ymm1));
			}
			_mm_sfence();
			for (; i != x.size(); ++i) {
				ret_vec.data.m_Re[i] = x.data.m_Re[i] / Re[i];
				ret_vec.data.m_Im[i] = x.data.m_Im[i] / Re[i];
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 4) {
				// Will unrolling 2x not saturate divider unit.
				// We have two parallel division so at least second
				// operation will be pipelined at divider level.
				const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i]);
				const __m256d ymm1 = _mm256_load_pd(&Re[i]);
				_mm256_store_pd(&ret_vec.data.m_Re[i], _mm256_div_pd(ymm0, ymm1));
				const __m256d ymm2 = _mm256_load_pd(&x.data.m_Im[i]);
				_mm256_store_pd(&ret_vec.data.m_Im[i], _mm256_div_pd(ymm2, ymm1));
			}
			for (; i != x.size(); ++i) {
				ret_vec.data.m_Re[i] = x.data.m_Re[i] / Re[i];
				ret_vec.data.m_Im[i] = x.data.m_Im[i] / Re[i];
			}
#endif
			return (ret_vec);
		}

		template<int32_t N> AVXShortCVec1D<N>
		inline operator==(const AVXShortCVec1D<N> &x,
				  const AVXShortCVec1D<N> &y) {
			using namespace gms::common;
			using namespace gms::math::constants;
			if (x.data.m_nsize != y.data.m_nsize) { return (AVXShortCVec1D<N>{}); }

			__attribute__((align(64))) struct {
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
			AVXShortCVec1D<N> ret_vec;
			int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
			for(i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize,4); i += 8) {
				 ca.ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
				 ca.ymm1(_mm256_load_pd(&y.data.m_Re[i+0]));
				 _mm256_stream_pd(&ret_vec.data.m_Re[i+0], _mm256_cmp_pd(ca.ymm0,ca.ymm1,_CMP_EQ_OQ));
				 ca.ymm2(_mm256_load_pd(&x.data.m_Re[i+4]));
				 ca.ymm3(_mm256_load_pd(&y.data.m_Re[i+4]));
				 _mm256_stream_pd(&ret_vec.data.m_Re[i+4], _mm256_cmp_pd(ca.ymm2,ca.ymm3,_CMP_EQ_OQ));
				 ca.ymm4(_mm256_load_pd(&x.data.m_Im[i+0]));
				 ca.ymm5(_mm256_load_pd(&y.data.m_Im[i+0]));
				 _mm256_stream_pd(&ret_vec.data.m_Im[i+0], _mm256_cmp_pd(ca.ymm4,ca.ymm5,_CMP_EQ_OQ));
				 ca.ymm6(_mm256_load_pd(&x.data.m_Im[i+4]));
				 ca.ymm7(_mm256_load_pd(&y.data.m_Im[i+1]));
				 _mm256_stream_pd(&ret_vec.data.m_Im[i+4],_mm256_cmp_pd(ca.ymm6,ca.ymm7,_CMP_EQ_OQ));
			}
			_mm_sfence();
			for(; i != ret_vec.data.m_nsize; ++i) {
				if(approximately_equalf64(x.data.m_Re[i],
					               y.data.m_Re[i],DEPS)) {
					ret_vec.data.m_Re[i] = 1.0;
				 }
				 else {
					 ret_vec.data.m_Re[i] = 0.0;
				 }
				 if (approximately_equalf64(x.data.m_Im[i],
								   y.data.m_Im[i],DEPS)) {
					 ret_vec.data.m_Im[i] = 1.0;
				}
				else {
					ret_vec.data.m_Im[i] = 0.0;
				}
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {
				ca.ymm0(_mm256_load_pd(&x.data.m_Re[i + 0]));
				ca.ymm1(_mm256_load_pd(&y.data.m_Re[i + 0]));
				_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_cmp_pd(ca.ymm0, ca.ymm1, _CMP_EQ_OQ));
				ca.ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
				ca.ymm3(_mm256_load_pd(&y.data.m_Re[i + 4]));
				_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_cmp_pd(ca.ymm2, ca.ymm3, _CMP_EQ_OQ));
				ca.ymm4(_mm256_load_pd(&x.data.m_Im[i + 0]));
				ca.ymm5(_mm256_load_pd(&y.data.m_Im[i + 0]));
				_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_cmp_pd(ca.ymm4, ca.ymm5, _CMP_EQ_OQ));
				ca.ymm6(_mm256_load_pd(&x.data.m_Im[i + 4]));
				ca.ymm7(_mm256_load_pd(&y.data.m_Im[i + 4]));
				_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_cmp_pd(ca.ymm6, ca.ymm7, _CMP_EQ_OQ));
			}
			
			for (; i != ret_vec.data.m_nsize; ++i) {
				if (approximately_equalf64(x.data.m_Re[i],
					y.data.m_Re[i], DEPS)) {
					ret_vec.data.m_Re[i] = 1.0;
				}
				else {
					ret_vec.data.m_Re[i] = 0.0;
				}
				if (approximately_equalf64(x.data.m_Im[i],
					y.data.m_Im[i], DEPS)) {
					ret_vec.data.m_Im[i] = 1.0;
				}
				else {
					ret_vec.data.m_Im[i] = 0.0;
				}
		  }
#endif
				return (ret_vec);
		}

		template<int32_t N> AVXShortCVec1D<N>
		inline operator!=(const AVXShortCVec1D<N> &x,
				  const AVXShortCVec1D<N> &y) {
			using namespace gms::common;
			using namespace gms::math::constants;
			if (x.data.m_nsize != y.data.m_nsize) { return (AVXShortCVec1D<N>{}); }

			__attribute__((align(64))) struct {
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
			AVXShortCVec1D<N> ret_vec;
			int32_t i;
#if (USE_AVXCOMPLEX_SHORTVEC_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {
				ca.ymm0(_mm256_load_pd(&x.data.m_Re[i + 0]));
				ca.ymm1(_mm256_load_pd(&y.data.m_Re[i + 0]));
				_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_cmp_pd(ca.ymm0, ca.ymm1, _CMP_NEQ_OQ));
				ca.ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
				ca.ymm3(_mm256_load_pd(&y.data.m_Re[i + 4]));
				_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_cmp_pd(ca.ymm2, ca.ymm3, _CMP_NEQ_OQ));
				ca.ymm4(_mm256_load_pd(&x.data.m_Im[i + 0]));
				ca.ymm5(_mm256_load_pd(&y.data.m_Im[i + 0]));
				_mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_cmp_pd(ca.ymm4, ca.ymm5, _CMP_NEQ_OQ));
				ca.ymm6(_mm256_load_pd(&x.data.m_Im[i + 4]));
				ca.ymm7(_mm256_load_pd(&y.data.m_Im[i + 1]));
				_mm256_stream_pd(&ret_vec.data.m_Im[i + 4], _mm256_cmp_pd(ca.ymm6, ca.ymm7, _CMP_NEQ_OQ));
			}
			_mm_sfence();
			for (; i != ret_vec.data.m_nsize; ++i) {
				if (!approximately_equalf64(x.data.m_Re[i],
					y.data.m_Re[i], DEPS)) {
					ret_vec.data.m_Re[i] = 1.0;
				}
				else {
					ret_vec.data.m_Re[i] = 0.0;
				}
				if (!approximately_equalf64(x.data.m_Im[i],
					y.data.m_Im[i], DEPS)) {
					ret_vec.data.m_Im[i] = 1.0;
				}
				else {
					ret_vec.data.m_Im[i] = 0.0;
				}
			}
#else
			for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {
				ca.ymm0(_mm256_load_pd(&x.data.m_Re[i + 0]));
				ca.ymm1(_mm256_load_pd(&y.data.m_Re[i + 0]));
				_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_cmp_pd(ca.ymm0, ca.ymm1, _CMP_NEQ_OQ));
				ca.ymm2(_mm256_load_pd(&x.data.m_Re[i + 4]));
				ca.ymm3(_mm256_load_pd(&y.data.m_Re[i + 4]));
				_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_cmp_pd(ca.ymm2, ca.ymm3, _CMP_NEQ_OQ));
				ca.ymm4(_mm256_load_pd(&x.data.m_Im[i + 0]));
				ca.ymm5(_mm256_load_pd(&y.data.m_Im[i + 0]));
				_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_cmp_pd(ca.ymm4, ca.ymm5, _CMP_NEQ_OQ));
				ca.ymm6(_mm256_load_pd(&x.data.m_Im[i + 4]));
				ca.ymm7(_mm256_load_pd(&y.data.m_Im[i + 4]));
				_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_cmp_pd(ca.ymm6, ca.ymm7, _CMP_NEQ_OQ));
			}

			for (; i != ret_vec.data.m_nsize; ++i) {
				if (!approximately_equalf64(x.data.m_Re[i],
					y.data.m_Re[i], DEPS)) {
					ret_vec.data.m_Re[i] = 1.0;
				}
				else {
					ret_vec.data.m_Re[i] = 0.0;
				}
				if (!approximately_equalf64(x.data.m_Im[i],
					y.data.m_Im[i], DEPS)) {
					ret_vec.data.m_Im[i] = 1.0;
				}
				else {
					ret_vec.data.m_Im[i] = 0.0;
				}
			}
#endif
			return (ret_vec);
		  }


		  template<int32_t N>
		  void v256scnormalize_product(AVXShortCVec1D<N> &out,
					       const AVXShortCVec1D<N> &v1,
					       const AVXShortCVec1D<N> &v2) {
			  avx256_cnormalize_prod<AVXShortCVec1D<N>>(out,v1,v2);
		  }

		  template<int32_t N>
		  void v256scmean_product(std::complex<double> &mean,
					  const AVXShortCVec1D<N> &v1,
					  const AVXShortCVec1D<N> &v2) {
			  avx256_cmean_prod<AVXShortCVec1D<N>>(mean,v1,v2);
		  }

		  template<int32_t N>
		  void v256scmean_quotient( std::complex<double> &mean,
					    const AVXShortCVec1D<N> &v1,
					    const AVXShortCVec1D<N> &v2) {
			  avx256_cmean_quot<AVXShortCVec1D<N>>(mean,v1,v2);
		  }

		  template<int32_t N>
		  void v256scconj_product(AVXShortCVec1D<N> &out,
					  const AVXShortCVec1D<N> &v1,
					  const AVXShortCVec1D<N> &v2,
					  const bool do_nt_store) {
			  avx256_cconj_prod<AVXShortCVec1D<N>>(out,v1,v2,do_nt_store);
		  }

		  template<int32_t N>
		  void v256scnorm_conjprod(AVXShortCVec1D<N> &out,
					   const AVXShortCVec1D<N> &v1,
					   const AVXShortCVec1D<N> &v2,
					   const bool do_nt_store) {
			  avx256_cnorm_conjprod<AVXShortCVec1D<N>>(out,v1,v2,do_nt_store);
		  }

		  template<int32_t N>
		  void v256scmean_conjprod( AVXShortCVec1D<N> &out,
					    const AVXShortCVec1D<N> &x,
					    const AVXShortCVec1D<N> &y) {
			  avx256_cmean_conjprod<AVXShortCVec1D<N>>(out,v1,v2);
		  }

		  template<int32_t N>
		  void v256sc_arithmean(std::complex<double> &mean,
					const AVXShortCVec1D<N> &v) {
			  avx256_arith_mean<AVXShortCVec1D<N>>(mean,v);
		  }

		  template<int32_t N>
		  void v256sc_normalize(AVXShortCVec1D<N> &norm,
					const AVXShortCVec1D<N> &v,
					const AVXShortCVec1D<N> &cv,
					const bool do_nt_store) {
			  avx256_cnormalize<AVXShortCVec1D<N>>(norm,v,cv,do_nt_store);
		  }

		  template<int32_t N>
		  void v256sc_magnitude( AVXShortCVec1D<N> &vmag,
					 const AVXShortCVec1D<N> &v,
					 const AVXShortCVec1D<N> &v2) {
			  avx256_cmagnitude<AVXShortCVec1D<N>>(vmag,v,v2);
		  }
	}
}


#endif /*__GMS_AVXCOMPLEX_SHORTVEC_H__*/
