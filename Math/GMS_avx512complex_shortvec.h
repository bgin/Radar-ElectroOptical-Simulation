
#ifndef __GMS_AVX512COMPLEX_SHORTVEC_H__
#define __GMS_AVX512COMPLEX_SHORTVEC_H__

namespace file_info {
#if defined _WIN64 
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

	const unsigned int gGMS_AVX512COMPLEX_SHORTVEC_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_AVX512COMPLEX_SHORTVEC_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_AVX512COMPLEX_SHORTVEC_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_AVX512COMPLEX_SHORTVEC_FULLVER = 
		1000U*gGMS_AVX512COMPLEX_SHORTVEC_MAJOR+100U*gGMS_AVX512COMPLEX_SHORTVEC_MINOR+10U*gGMS_AVX512COMPLEX_SHORTVEC_MICRO;

	const char * const pgGMS_AVX512COMPLEX_SHORTVEC_CREATE_DATE = "12-10-2019 13:12 +00200 (SAT 12 OCT 2019 GMT+2)";

	const char * const pgGMS_AVX512COMPLEX_SHORTVEC_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_AVX512COMPLEX_SHORTVEC_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVX512COMPLEX_SHORTVEC_SYNOPSIS = "AVX512 complex vector (1D) stack-allocated storage.";
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
#if defined _WIN64
    #include "../GMS_config.h"
    #include "../GMS_common.h"
#elif defined __linux
    #include "GMS_config.h"
    #include "GMS_common.h"
#endif
#include "LAM_avx512complex_common.h"
#if defined _WIN64
    #include "../Math/GMS_constants.h"
#elif defined __linux
    #include "GMS_constants.h"
#endif

#if !defined (USE_AVX512COMPLEX_SHORTV_NT_STORES) // Streaming stores defined per this struct (default set to 0)
#define USE_AVX512COMPLEX_SHORTV_NT_STORES 0
#endif

namespace gms {
	namespace math {

#if !defined (AVX512COMPLEX_SHORTV_LOAD_ZMM)
#define AVX3COMPLEX_SHORTV_LOAD_ZMM(reg1,reg2,reg3,reg4,v1,v2,idx,off) \
	(reg1) = _mm512_load_pd(&(v1).data.m_Re[(idx)+(off)]);			   \
	(reg2) = _mm512_load_pd(&(v2).data.m_Re[(idx)+(off)]);			   \
	(reg3) = _mm512_load_pd(&(v1).data.m_Im[(idx)+(off)]);			   \
	(reg4) = _mm512_load_pd(&(v2).data.m_Im[(idx)+(off)]);
#endif
#if defined _WIN64
		template<int32_t N>
		  __declspec(align(64))  struct AVX512SCVData {
#elif defined
		template<int32_t N>
		__attribute__((align(64))) struct AVX512SCVData {
#endif
			static  constexpr int32_t MAX_SIZE = 4096;
			int32_t m_nsize = N;
			PAD_TO(1,4)
			PAD_TO(2,4)
			
			static_assert(N <= MAX_SIZE, "Invalid size of AVX512SCVData -- passed!!");
#if defined _WIN64
			__declspec(align(64))	double m_Re[(N == 0) ? 8 : N]; // Smallest size 8 elements

			__declspec(align(64))   double m_Im[(N == 0) ? 8 : N];
#elif defined __linux
		        __attribute__((align(64))) double m_Re[(N == 0) ? 8 : N];
			__attribute__((align(64))) double m_Im[(N == 0) ? 8 : N];
#endif
		};

#if defined _WIN64
		template<int32_t N> 
		struct  __declspec(align(64)) AVX512SmallCVec1D {
#elif defined __linux
		template<int32_t N>
		__attribute__((align(64)))  struct AVX512SmallCVec1D {
#endif
		        AVX512SCVData<N> data;
			
			AVX512SmallCVec1D() noexcept(true) {

				data.m_Re[N] = {};
				data.m_Im[N] = {};
				data.m_ntstore = false;
			}

			AVX512SmallCVec1D(const double (&Re)[N],
				          const double (&Im)[N]
							  ) {
				using namespace gms::common;
				
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1	
				avx512_uncached_memmove(&data.m_Re[0],Re,N);	
				avx512_uncached_memmove(&data.m_Im[0],Im,N);
#else
				avx512_cached_memmove(&data.m_Re[0],Re,N);
				avx512_cached_memmove(&data.m_Im[0],Im,N);
#endif					
			}

			AVX512SmallCVec1D(const AVX512SmallCVec1D &x)  {
				using namespace gms::common;
				data.m_nsize = x.data.m_nsize;
				
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
				avx512_uncached_memmove(&data.m_Re[0], &x.data.m_Re[0], N);
				avx512_uncached_memmove(&data.m_Im[0], &x.data.m_Im[0], N);
#else
				avx512_cached_memmove(&data.m_Re[0], &x.data.m_Re[0], N);
				avx512_cached_memmove(&data.m_Im[0], &x.data.m_Im[0], N);
#endif
			}

		
		    // operator= not needed here, although implemented
			AVX512SmallCVec1D & 
			operator=(const AVX512SmallCVec1D &x){
			    using namespace GMS::common;
				if (this == &x || data.m_nsize != x.data.m_nsize){
				    return (*this);
				}
				// Destructive copy
				
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
				avx512_uncached_memmove(&data.m_Re[0],&x.data.m_Re[0], N);
				avx512_cached_memmove(&data.m_Im[0],&x.data.m_Im[0], N);
#else
				avx512_cached_memmove(&data.m_Re[0], &x.data.m_Re[0], N);
				avx512_cached_memmove(&data.m_Im[0], &x.data.m_Im[0], N);
#endif
				return (*this);
			}

			inline int32_t  size() { return (data.m_nsize); };

		};

		// Global operators

		template<int32_t N> std::ostream & 
		operator<<(std::ostream &os,
			   const AVX512SmallCVec1D<N> &x) {
			for (int32_t i = 0; i != x.data.m_nsize; ++i) {
				os << std::fixed << std::showpoint << std::setprecision(15) <<
					std::setw(4)  <<  "Re: " << "{"  << x.data.m_Re[i] << "}" <<
					std::setw(12) << "Im: "  << "{"  << x.data.m_Im[i] << "}" << std::endl;
			}
			return (os);
		}

		template<int32_t N> AVX512SmallCVec1D<N>
		inline operator+(const AVX512SmallCVec1D<N> &x,
				 const AVX512SmallCVec1D<N> &y) {
			if (x.size() != y.size()) { return (AVX512SmallCVec1D<N>{}); }
			AVX512SmallCVec1D<N> ret_vec;
			int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.size(), 8); i += 16) {
				// Linearly growing indices, no need for software prefetch.
				// HW prefetch will kick in after 2 maybe 3 cache misses.
				const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
				const __m512d zmm1 = _mm512_load_pd(&y.data.m_Re[i + 0]);
				_mm512_stream_pd(&ret_vec.data.m_Re[i], _mm512_add_pd(zmm0, zmm1));
				const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
				const __m512d zmm3 = _mm512_load_pd(&y.data.m_Re[i + 8]);
				_mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
				const __m512d zmm4 = _mm512_load_pd(&x.data.m_Im[i + 0]);
				const __m512d zmm5 = _mm512_load_pd(&y.data.m_Im[i + 0]);
				_mm512_stream_pd(&ret_vec.data.m_Im[i + 0], _mm512_add_pd(zmm4, zmm5));
				const __m512d zmm6 = _mm512_load_pd(&x.data.m_Im[i + 8]);
				const __m512d zmm7 = _mm512_load_pd(&y.data.m_Im[i + 8]);
				_mm512_stream_pd(&ret_vec.data.m_Im[i + 8], _mm512_add_pd(zmm6, zmm7));
			}
			_mm_sfence();
			for (; i != ret_vec.data.m_nsize; ++i) {
				ret_vec.data.m_Re[i] = x.data.m_Re[i] + y.data.m_Re[i];
				ret_vec.data.m_Im[i] = x.data.m_Im[i] + y.data.m_Re[i];
			}
#else			
			for (i = 0; i != ROUND_TO_EIGHT(ret_vec.size(), 8); i += 16) {
				// Linearly growing indices, no need for software prefetch.
				// HW prefetch will kick in after 2 maybe 3 cache misses.
				const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
				const __m512d zmm1 = _mm512_load_pd(&y.data.m_Re[i + 0]);
				_mm512_store_pd(&ret_vec.data.m_Re[i], _mm512_add_pd(zmm0, zmm1));
				const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
				const __m512d zmm3 = _mm512_load_pd(&y.data.m_Re[i + 8]);
				_mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
				const __m512d zmm4 = _mm512_load_pd(&x.data.m_Im[i + 0]);
				const __m512d zmm5 = _mm512_load_pd(&y.data.m_Im[i + 0]);
				_mm512_store_pd(&ret_vec.data.m_Im[i + 0], _mm512_add_pd(zmm4, zmm5));
				const __m512d zmm6 = _mm512_load_pd(&x.data.m_Im[i + 8]);
				const __m512d zmm7 = _mm512_load_pd(&y.data.m_Im[i + 8]);
				_mm512_store_pd(&ret_vec.data.m_Im[i + 8], _mm512_add_pd(zmm6, zmm7));
			}
			for (; i != ret_vec.data.m_nsize; ++i) {
				ret_vec.data.m_Re[i] = x.data.m_Re[i] + y.data.m_Re[i];
				ret_vec.data.m_Im[i] = x.data.m_Im[i] + y.data.m_Re[i];
			}
#endif				
			return (ret_vec);
				
		}
				


			




			

				


				

				


				


			
			
			
		

		template<int32_t N> AVX512SmallCVec1D<N>
		inline   operator+(const AVX512SmallCVec1D<N> &x,
				  const double __restrict Re[N]) {         // If Re is not equal to x --> udefined behaviour.
				  using namespace gms::common;
				  if (!Is_ptr_aligned64(Re)) { return AVX512SmallCVec1D<N>{}; }
				  AVX512SmallCVec1D<N> ret_vec;
				  int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
					  // Linearly growing indices, no need for software prefetch.
					  // HW prefetch will kick in after 2 maybe 3 cache misses.
					  const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
					  const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					  _mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
					  const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					  _mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
					} 
				  _mm_sfence();
				  for (; i != ret_vec.size(); ++i) {
					  ret_vec.data.m_Re[i] = x.data.m_Re[i] + Re[i];
				  }
#else
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
					  // Linearly growing indices, no need for software prefetch.
					  // HW prefetch will kick in after 2 maybe 3 cache misses.

					  const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
					  const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					  _mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
					  const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					  _mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));

				  }
				  for (; i != ret_vec.size(); ++i) {
					  ret_vec.data.m_Re[i] = x.data.m_Re[i] + Re[i];
				  }
#endif
				  return (ret_vec);
		 }
				 
				 
		

		template<int32_t N> AVX512SmallCVec1D<N>
		inline operator+(const double __restrict Re[N],
				 const AVX512SmallCVec1D<N> &x) {
				  using namespace gms::common;
				  if (!Is_ptr_aligned64(Re)) { return (AVX512SmallCVec1D<N>{}); }
				  AVX512SmallCVec1D<N> ret_vec;
				  int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					  const __m512d zmm0(_mm512_load_pd(&Re[i + 0]));
					  const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
					  _mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2(_mm512_load_pd(&Re[i + 8]));
					  const __m512d zmm3(_mm512_load_pd(&x.data.m_Re[i + 8]));
					  _mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
					}
				  _mm_sfence();
				  for (; i != x.size(); ++i) {
					  ret_vec.data.m_Re[i] = Re[i] + x.data.m_Re[i];
				  }
#else
				  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
					  const __m512d zmm0(_mm512_load_pd(&Re[i + 0]));
					  const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
					  _mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
					  const __m512d zmm2(_mm512_load_pd(&Re[i + 8]));
					  const __m512d zmm3(_mm512_load_pd(&x.data.m_Re[i + 8]));
					  _mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
					} 
				  for (; i != x.size(); ++i) {
					  ret_vec.data.m_Re[i] = Re[i] + x.data.m_Re[i];
				  }
#endif				  
				  return (ret_vec);
	       }
				 
				 
				
		  

		  template<int32_t> AVX512SmallCVec1D<N>
		  inline operator+(const AVX512SmallCVec1D<N> &x,
				   const double c) {
					AVX512SmallCVec1D<N> ret_vec;
					int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0(_mm512_set1_pd(c));
						const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
						_mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
						const __m512d zmm2(_mm512_load_pd(&x.data.m_Re[i+8]));
						_mm512_stream_pd(&ret_vec.data.m_Re[i+8],_mm512_add_pd(zmm0,zmm2));
						}
					_mm_sfence();
					for (; i != x.size(); ++i) {
						ret_vec.data.m_Re[i] = x.data.m_Re[i] + c;
					}
#else
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0(_mm512_set1_pd(c));
						const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
						_mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
						const __m512d zmm2(_mm512_load_pd(&x.data.m_Re[i + 8]));
						_mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_add_pd(zmm0, zmm2));
					}

					for (; i != x.size(); ++i) {
						ret_vec.data.m_Re[i] = x.data.m_Re[i] + c;
					}
#endif
					return (ret_vec);
			}		
		   

		   template<int32_t N> AVX512SmallCVec1D<N>
		   inline operator-(const AVX512SmallCVec1D<N> &x,
				    const AVX512SmallCVec1D<N> &y) {
					 using namespace gms::common;
					 if (x.size() != y.size()) { return (AVX512SmallCVec1D<N>{}); }
					 AVX512SmallCVec1D<N> ret_vec;
					 int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
					 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						 const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
						 const __m512d zmm1 = _mm512_load_pd(&y.data.m_Re[i + 0]);
						 _mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
						 const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
						 const __m512d zmm3 = _mm512_load_pd(&y.data.m_Re[i + 8]);
						 _mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
						 const __m512d zmm4 = _mm512_load_pd(&x.data.m_Im[i + 0]);
						 const __m512d zmm5 = _mm512_load_pd(&y.data.m_Im[i + 0]);
						 _mm512_stream_pd(&ret_vec.data.m_Im[i + 0], _mm512_sub_pd(zmm4, zmm5));
						 const __m512d zmm6 = _mm512_load_pd(&x.data.m_Im[i + 8]);
						 const __m512d zmm7 = _mm512_load_pd(&y.data.m_Im[i + 8]);
						 _mm512_stream_pd(&ret_vec.data.m_Im[i + 8], _mm512_sub_pd(zmm6, zmm7));
					 }
					 _mm_sfence();
					 for (; i != x.size(); ++i) {
						 ret_vec.data.m_Re[i] = x.data.m_Re[i] - y.data.m_Re[i];
						 ret_vec.data.m_Im[i] = x.data.m_Im[i] - y.data.m_Im[i];
					 }
#else
						 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

							 const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
							 const __m512d zmm1 = _mm512_load_pd(&y.data.m_Re[i + 0]);
							 _mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
							 const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
							 const __m512d zmm3 = _mm512_load_pd(&y.data.m_Re[i + 8]);
							 _mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
							 const __m512d zmm4 = _mm512_load_pd(&x.data.m_Im[i + 0]);
							 const __m512d zmm5 = _mm512_load_pd(&y.data.m_Im[i + 0]);
							 _mm512_store_pd(&ret_vec.data.m_Im[i + 0], _mm512_sub_pd(zmm4, zmm5));
							 const __m512d zmm6 = _mm512_load_pd(&x.data.m_Im[i + 8]);
							 const __m512d zmm7 = _mm512_load_pd(&y.data.m_Im[i + 8]);
							 _mm512_store_pd(&ret_vec.data.m_Im[i + 8], _mm512_sub_pd(zmm6, zmm7));
						}
						 for (; i != x.size(); ++i) {
							 ret_vec.data.m_Re[i] = x.data.m_Re[i] - y.data.m_Re[i];
							 ret_vec.data.m_Im[i] = x.data.m_Im[i] - y.data.m_Im[i];
						 }
#endif
					 return (ret_vec);
			}		 
					
		  

		  template<int32_t N> AVX512SmallCVec1D<N>
		  inline operator-(AVX512SmallCVec1D<N> &x,
				   const double __restrict Re[N]) {
					using namespace gms::common;
					if (!Is_ptr_aligned64(Re)) { return (AVX512SmallCVec1D<N>{}); }
					AVX512SmallCVec1D<N> ret_vec;
					int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
						const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
						_mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
						const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
						const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
						_mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));

					}
					_mm_sfence();
					for (; i != x.size(); ++i) {
						ret_vec.data.m_Re[i] = x.data.m_Re[i] - Re[i];
					}
#else
					for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
						const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
						_mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
						const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
						const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
						_mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
					}

					for (; i != x.size(); ++i) {
						ret_vec.data.m_Re[i] = x.data.m_Re[i] - Re[i];
					}
#endif
					
					return (ret_vec);
			}

			template<int32_t N> AVX512SmallCVec1D<N>
			inline operator-(  const AVX512SmallCVec1D<N> &x,
					   const double c) {
					  AVX512SmallCVec1D<N> ret_vec;
					  int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
					  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						  const __m512d zmm0(_mm512_set1_pd(c));
						  const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
						  _mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(zmm1, zmm0));
						  const __m512d zmm2(_mm512_load_pd(&x.data.m_Re[i+8]));
						  _mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(zmm2,zmm0));
						  }
					  _mm_sfence();
					  for (; i != x.size(); ++i) {
						  ret_vec.data.m_Re[i] = x.data.m_Re[i] - c;
					  }
#else
					  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

						  const __m512d zmm0(_mm512_set1_pd(c));
						  const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
						  _mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(zmm1, zmm0));
						  const __m512d zmm2(_mm512_load_pd(&x.data.m_Re[i + 8]));
						  _mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm0));
					  }

					  for (; i != x.size(); ++i) {
						  ret_vec.data.m_Re[i] = x.data.m_Re[i] - c;
					  }
#endif
					 
					  return (ret_vec);
		   }

		   template<int32_t N> AVX512SmallCVec1D<N>
		   inline operator*(const AVX512SmallCVec1D<N> &x,
				    const AVX512SmallCVec1D<N> &y) {
			   if (x.size() != y.size()) { return (AVX512SmallCVec1D<N>{}); }
			    __m512d zmm0,zmm1,zmm2,zmm3,zmm4,zmm5,zmm6,zmm7;
			   AVX512SmallCVec1D<N> ret_vec;
			   int32_t i; 
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
			   for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm0, zmm1, zmm2, zmm3, x, y, i, 0)
					   _mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(
					   _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm3)));
				   _mm512_stream_pd(&ret_vec.data.m_Im[i + 0], _mm512_add_pd(
					   _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));

				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm4, zmm5, zmm6, zmm7, x, y, i, 8)
					   _mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(
					   _mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6, zmm7)));
				   _mm512_stream_pd(&ret_vec.data.m_Im[i + 8], _mm512_add_pd(
					   _mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));
			   }    
			   _mm_sfence();
			   for (; i != x.size(); ++i) {
				   ret_vec.data.m_Re[i] = (x.data.m_Re[i] * y.data.m_Re[i]) - (x.data.m_Im[i] * y.data.m_Im[i]);
				   ret_vec.data.m_Im[i] = (x.data.m_Im[i] * y.data.m_Re[i]) + (x.data.m_Re[i] * y.data.m_Im[i]);
			   }
#else
			   for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {
				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm0, zmm1, zmm2, zmm3, x, y, i, 0)
					   _mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_sub_pd(
					   _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm3)));
				   _mm512_store_pd(&ret_vec.data.m_Im[i + 0], _mm512_add_pd(
					   _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));

				   AVX3COMPLEX_SMALLV_LOAD_ZMM(zmm4, zmm5, zmm6, zmm7, x, y, i, 8)
					   _mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_sub_pd(
					   _mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6, zmm7)));
				   _mm512_store_pd(&ret_vec.data.m_Im[i + 8], _mm512_add_pd(
					   _mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));
				}
			   for (; i != x.size(); ++i) {
				   ret_vec.data.m_Re[i] = (x.data.m_Re[i] * y.data.m_Re[i]) - (x.data.m_Im[i] * y.data.m_Im[i]);
				   ret_vec.data.m_Im[i] = (x.data.m_Im[i] * y.data.m_Re[i]) + (x.data.m_Re[i] * y.data.m_Im[i]);
			   }
#endif				 
			   return (ret_vec);
		}	  

			 
			  
		 


			  
			 
			 

			 template<int32_t N> AVX512SmallCVec1D<N>
			 inline operator*(const AVX512SmallCVec1D<N> &x,
					  const double __restrict Re[N]) {
				 using namespace gms::common;
				 if (!Is_ptr_aligned64(Re)) { return (AVX512SmallCVec1D<N>{}); }
				 AVX512SmallCVec1D<N> ret_vec;
				 int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
				 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					 const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
					 const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					 _mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
					 const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
					 const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					 _mm512_stream_pd(&ret_vec.data.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
					 const __m512d zmm4 = _mm512_load_pd(&x.data.m_Im[i + 0]);
					 _mm512_stream_pd(&ret_vec.data.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
					 const __m512d zmm5 = _mm512_load_pd(&x.data.m_Im[i + 8]);
					 _mm512_stream_pd(&ret_vec.data.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
					 }
				 _mm_sfence();
				 for (; i != ret_vec.data.m_nsize; ++i) {
					 ret_vec.data.m_Re[i] = x.data.m_Re[i] * Re[i];
					 ret_vec.data.m_Im[i] = x.data.m_Im[i] * Re[i];
				 }
#else
				 for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					 const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
					 const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
					 _mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
					 const __m512d zmm2 = _mm512_load_pd(&x.data.m_Re[i + 8]);
					 const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
					 _mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
					 const __m512d zmm4 = _mm512_load_pd(&x.data.m_Im[i + 0]);
					 _mm512_store_pd(&ret_vec.data.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
					 const __m512d zmm5 = _mm512_load_pd(&x.data.m_Im[i + 8]);
					 _mm512_store_pd(&ret_vec.data.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
				 }
				 for (; i != ret_vec.data.m_nsize; ++i) {
					 ret_vec.data.m_Re[i] = x.data.m_Re[i] * Re[i];
					 ret_vec.data.m_Im[i] = x.data.m_Im[i] * Re[i];
				 }
#endif
				 return (ret_vec);
		 }
				
				 
			

			template<int32_t N> AVX512SmallCVec1D<N>
			inline operator*(const AVX512SmallCVec1D<N> &x,
					 const double c) {
				AVX512SmallCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVX512COMPLEX_SMALLV_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					const __m512d zmm0(_mm512_set1_pd(c));
					const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
					_mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_mul_pd(zmm1, zmm0));
					const __m512d zmm2(_mm512_load_pd(&x.data.m_Re[i+8]));
					_mm512_stream_pd(&ret_vec.data.m_Re[i+8], _mm512_mul_pd(zmm2,zmm0));

				}
				_mm_sfence();
				for (; i != x.size(); ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] * c;
				}
#else
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 16) {

					const __m512d zmm0(_mm512_set1_pd(c));
					const __m512d zmm1(_mm512_load_pd(&x.data.m_Re[i + 0]));
					_mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_mul_pd(zmm1, zmm0));
					const __m512d zmm2(_mm512_load_pd(&x.data.m_Re[i + 8]));
					_mm512_store_pd(&ret_vec.data.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm0));
				}

				for (; i != x.size(); ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] * c;
				}
#endif
				return (ret_vec);
			}	
		  

		  template<int32_t N> AVX512SmallCVec1D<N>
		  inline operator/(const AVX512SmallCVec1D<N> &x,
				   const AVX512SmallCVec1D<N> &y) {
			  if (x.size() != y.size()) { return (AVX512SmallCVec1D<N>{}); }
			  AVX512SmallCVec1D<N> ret_vec;
			  int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
			  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
				  // Will unrolling 2x not saturate divider unit.
				  // We have two parallel division so at least second
				  // operation will be pipelined at divider level.

				  const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
				  const __m512d zmm1 = _mm512_load_pd(&y.data.m_Im[i + 0]);
				  const __m512d zmm2 = _mm512_load_pd(&x.data.m_Im[i + 0]);
				  const __m512d re_term1 = _mm512_add_pd(
					  _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm1));
				  const __m512d re_term2 = _mm512_add_pd(
					  _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm1));
				  const __m512d zmm3 = _mm512_load_pd(&y.data.m_Re[i + 0]);
				  const __m512d den_term = _mm512_add_pd(
					  _mm512_mul_pd(zmm3, zmm3), _mm512_mul_pd(zmm1, zmm1));

				  _mm512_stream_pd(&ret_vec.data.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
				  _mm512_stream_pd(&ret_vec.data.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));
				  }
			  _mm_sfence();
			  for (; i != x.size(); ++i) {
				  const double tre = (x.data.m_Re[i] * y.data.m_Im[i]) + (x.data.m_Im[i] * y.data.m_Im[i]);
				  const double tim = (x.data.m_Im[i] * y.data.m_Im[i]) - (x.data.m_Re[i] * y.data.m_Im[i]);
				  const double den = (y.data.m_Re[i] * y.data.m_Re[i]) + (y.data.m_Im[i] * y.data.m_Im[i]);
				  ret_vec.data.m_Re[i] = tre / den;
				  ret_vec.data.m_Im[i] = tim / den;
			  }
#else
			  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
				  // Will unrolling 2x not saturate divider unit.
				  // We have two parallel division so at least second
				  // operation will be pipelined at divider level.

				  const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i + 0]);
				  const __m512d zmm1 = _mm512_load_pd(&y.data.m_Im[i + 0]);
				  const __m512d zmm2 = _mm512_load_pd(&x.data.m_Im[i + 0]);
				  const __m512d re_term1 = _mm512_add_pd(
					  _mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm1));
				  const __m512d re_term2 = _mm512_add_pd(
					  _mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm1));
				  const __m512d zmm3 = _mm512_load_pd(&y.data.m_Re[i + 0]);
				  const __m512d den_term = _mm512_add_pd(
					  _mm512_mul_pd(zmm3, zmm3), _mm512_mul_pd(zmm1, zmm1));

				  _mm512_store_pd(&ret_vec.data.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
				  _mm512_store_pd(&ret_vec.data.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));
			  }
			  for (; i != x.size(); ++i) {
				  const double tre = (x.data.m_Re[i] * y.data.m_Im[i]) + (x.data.m_Im[i] * y.data.m_Im[i]);
				  const double tim = (x.data.m_Im[i] * y.data.m_Im[i]) - (x.data.m_Re[i] * y.data.m_Im[i]);
				  const double den = (y.data.m_Re[i] * y.data.m_Re[i]) + (y.data.m_Im[i] * y.data.m_Im[i]);
				  ret_vec.data.m_Re[i] = tre / den;
				  ret_vec.data.m_Im[i] = tim / den;
			  }
#endif
			  return (ret_vec);
		}
			 
			 
		 	 
			 
		 template<int32_t N> AVX512SmallCVec1D<N>
		 inline operator/(const AVX512SmallCVec1D<N> &x,
				  const double  __restrict Re[N]) {
				using namespace gms::common;
				if (!Is_ptr_aligned64(Re)) { return (AVX512SmallCVec1D<N>{}); }
				AVX512SmallCVec1D<N> ret_vec;
				int32_t i;
#if (USE_AVX512COMPLEX_SHORTV_NT_STORES) == 1
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
					// Will unrolling 2x not saturate divider unit.
					// We have two parallel division so at least second
					// operation will be pipelined at divider level.
					const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i]);
					const __m512d zmm1 = _mm512_load_pd(&Re[i]);
					_mm512_stream_pd(&ret_vec.data.m_Re[i], _mm512_div_pd(zmm0, zmm1));
					const __m512d zmm2 = _mm512_load_pd(&x.data.m_Im[i]);
					_mm512_stream_pd(&ret_vec.data.m_Im[i], _mm512_div_pd(zmm2, zmm1));
				}	
				_mm_sfence();
				for (; i != x.size(); ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] / Re[i];
					ret_vec.data.m_Im[i] = x.data.m_Im[i] / Re[i];
				}	
#else
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
					// Will unrolling 2x not saturate divider unit.
					// We have two parallel division so at least second
					// operation will be pipelined at divider level.
					const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i]);
					const __m512d zmm1 = _mm512_load_pd(&Re[i]);
					_mm512_store_pd(&ret_vec.data.m_Re[i], _mm512_div_pd(zmm0, zmm1));
					const __m512d zmm2 = _mm512_load_pd(&x.data.m_Im[i]);
					_mm512_store_pd(&ret_vec.data.m_Im[i], _mm512_div_pd(zmm2, zmm1));
				}
				for (; i != x.size(); ++i) {
					ret_vec.data.m_Re[i] = x.data.m_Re[i] / Re[i];
					ret_vec.data.m_Im[i] = x.data.m_Im[i] / Re[i];
				}
#endif
				return (ret_vec);
			}
				
					
				

				
				
				
		  

		 
		  template<int32_t N>
		  std::pair<std::array<__mmask8,Size>,
		  std::array<__mmask8,Size>>
		  inline   operator==(const AVX512SmallCVec1D<N> &x,
				      const AVX512SmallCVec1D<N> &y) {
			  using namespace gms::common;
			  using namespace gms::math::constants;
			 
			  if (x.size() != y.size()) {
				  return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
			  }
			  int32_t i;
			  int32_t k1 = 0, k2 = 0;
			  std::pair<std::array<__mmask8,N>,
						std::array<__mmask8,N>> ret_val;
			  for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
				  const __m512d zmm0   = _mm512_load_pd(&x.data.m_Re[i]);
				  const __m512d zmm1   = _mm512_load_pd(&y.data.m_Re[i]);
				  ret_val.first[++k1]  = _mm512_cmp_pd_mask(zmm0, zmm1, _CMP_EQ_OQ);
				  const __m512d zmm2   = _mm512_load_pd(&x.data.m_Im[i]);
				  const __m512d zmm3   = _mm512_load_pd(&y.data.m_Im[i]);
				  ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2, zmm3, _CMP_EQ_OQ);
			  }
			  unsigned char t1 = 0x00, t2 = 0x00;
			  k1 += 1;
			  k2 += 1;
			  for (; i != x.size(); ++i) {
				  if (approximately_equalf64(x.data.m_Re[i],
					  y.data.m_Re[i], DEPS)) {
					  t1 |= 1 << i;
				  }
				  if (approximately_equalf64(x.data.m_Im[i],
					  y.data.m_Im[i], DEPS)) {
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
		    inline operator!=(const AVX512SmallCVec1D<N> &x,
				      const AVX512SmallCVec1D<N> &y) {
				using namespace gms::common;
				using namespace gms::math::constants;
				std::pair<std::array<__mmask8,N>,
				std::array<__mmask8,N>> ret_val;
				if (x.size() != y.size()) {
					return (std::make_pair(std::array<__mmask8, N>{}, std::array<__mmask8, N>{}));
				}
				int32_t i;
				int32_t k1 = 0, k2 = 0;
				for (i = 0; i != ROUND_TO_EIGHT(x.size(), 8); i += 8) {
					const __m512d zmm0 = _mm512_load_pd(&x.data.m_Re[i]);
					const __m512d zmm1 = _mm512_load_pd(&y.data.m_Re[i]);
					ret_val.first[++k1] = _mm512_cmp_pd_mask(zmm0, zmm1, _CMP_NEQ_OQ);
					const __m512 zmm2 = _mm512_load_pd(&x.data.m_Im[i]);
					const __m512 zmm3 = _mm512_load_pd(&y.data.m_Im[i]);
					ret_val.second[++k2] = _mm512_cmp_pd_mask(zmm2, zmm3, _CMP_NEQ_OQ);
				}
				unsigned char t1 = 0x00, t2 = 0x00;
				k1 += 1;
				k2 += 1;
				for (; i != x.size(); ++i) {
					if (!approximately_equalf64(x.data.m_Re[i],
						y.data.m_Re[i], DEPS)) {
						t1 |= 1 << i;
					}
					if (!approximately_equalf64(x.data.m_Im[i],
						y.data.m_Im[i], DEPS)) {
						t2 |= 1 << i;
					}
				}
				ret_val.first[k1] = t1;
				ret_val.second[k2] = t2;
				return (ret_val);
			}

			template<int32_t N>
			void v512scnormalize_product(AVX512SmallCVec1D<N> &vout,
						     const AVX512SmallCVec1D<N> &v1,
						     const AVX512SmallCVec1D<N> &v2,
						     const bool do_nt_stream) {

				avx512_cnormalize_prod<AVX512SmallCVec1D<N>>(vout,v1,v2,do_nt_stream);
			}

			template<int32_t N>
			void v512scmean_product(std::complex<double> &mean,
						const AVX512SmallCVec1D<N> &v1,
						const AVX512SmallCVec1D<N> &v2) {
				
				avx512_cmean_prod<AVX512SmallCVec1D<N>>(mean,v1,v2);
			}

			template<int32_t N>
			void v512scmean_quotient(std::complex<double> &mean,
						 const AVX512SmallCVec1D<N> &v1,
						 const AVX512SmallCVec1D<N> &v2) {

				avx512_cmean_quot<AVX512SmallCVec1D<N>>(mean,v1,v2);
			}

			template<int32_t N>
			void v512scconj_product(AVX512SmallCVec1D<N> &vout,
						const AVX512SmallCVec1D<N> &v1,
						const AVX512SmallCVec1D<N> &v2,
						const bool do_nt_store) {
				
				avx512_cconj_prod<AVX512SmallCVec1D<N>>(vout,v1,v2,do_nt_store);
			}

			template<int32_t N>
			void v512scnorm_conjprod(AVX512SmallCVec1D<N> &vout,
						 const AVX512SmallCVec1D<N> &v1,
						 const AVX512SmallCVec1D<N> &v2,
						 const bool do_nt_store) {
				
				avx512_cnorm_conjprod<AVX512SmallCVec1D<N>>(vout,v1,v2,do_nt_store);
			}

			template<int32_t N>
			void v512scmean_conjprod(std::complex<double> &mean,
						 const AVX512SmallCVec1D<N> &v1,
						 const AVX512SmallCVec1D<N> &v2) {
				
				avx512_cmean_conjprod<AVX512SmallCVec1D<N>>(mean,v1,v2);
			}

			template<int32_t N>
			void v512sc_arithmean(std::complex<double> &mean,
					      const AVX512SmallCVec1D<N> &v) {
				
				avx512_arith_mean<AVX512SmallCVec1D<N>>(mean,v);
			}

			template<int32_t N>
			void v512sc_normalize(AVX512SmallCVec1D<N> &norm,
					      const AVX512SmallCVec1D<N> &v,
					      const AVX512SmallCVec1D<N> &cv,
					      const bool do_nt_store) {
					
				avx512_cnormalize<AVX512SmallCVec1D<N>>(norm,v,cv,do_nt_store);
			}

			template<int32_t N>
			void v512sc_magnitude(const double * __restrict vmag,
					      const AVX512SmallCVec1D<N> &v,
					      const AVX512SmallCVec1D<N> &cv) {
				
				avx512_cmagnitude<AVX512SmallCVec1D<N>>(vmag,v,cv);
			}

	}
}




#endif /*__GMS_AVX512COMPLEX_SHORTVEC_H__*/
