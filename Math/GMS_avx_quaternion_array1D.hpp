

#ifndef __GMS_AVX_QUATERNION_ARRAY1D_HPP__
#define __GMS_AVX_QUATERNION_ARRAY1D_HPP__


namespace file_info {

const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_MAJOR = 1U;
const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_MINOR = 0U;
const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_MICRO = 0U;
const unsigned int gGMS_AVX_QUATERNION_ARRAY1D_FULLVER =
   1000U*gGMS_AVX_QUATERNION_ARRAY1D_MAJOR+
   100U*gGMS_AVX_QUATERNION_ARRAY1D_MINOR+
   10U*gGMS_AVX_QUATERNION_ARRAY1D_MICRO;

const char * const pgGMS_AVX_QUATERNION_ARRAY1D_CREATION_DATE = "01-08-2021 09:58 +00200 (SUN 01 AUG 2021 GMT+2)";
const char * const pgGMS_AVX_QUATERNION_ARRAY1D_BUILD_DATE = __DATE__ ":" __TIME__;
const char * const pgGMS_AVX_QUATERNION_ARRAY1D_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_AVX_QUATERNION_ARRAY1D_DESCRIPTION = "Single-dimensional array of quaternionic vectors."
}


#include <cstdint>
#include <omp.h>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_malloc.h"
#include "GMS_common.h"
#include "GMS_constants.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined(USE_AVXQUATERNION_ARRAY1D_NT_STORES)
    #define USE_AVXQUATERNION_ARRAY1D_NT_STORES 0
#endif

// Enable/disable software prefetching.
#if !defined(AVXQUATERNION_ARRAY1D_PREFETCH_FROM_OBJ)
#define AVXQUATERNION_ARRAY1D_PREFETCH_FROM_OBJ(obj,idx,off,hint)  \
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_x[idx+off],hint);
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_y[idx+off],hint);
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_z[idx+off],hint);
   _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_w[idx+off],hint);
#endif

#if !defined(AVXQUATERNION_ARRAY1D_PREFETCH_FROM_PTR)
#define AVXQUATERNION_ARRAY1D_PREFETCH_FROM_PTR(ptr,idx,off,hint)  \
   _mm_prefetch(reinterpret_cast<const char*>(&(ptr)[idx+off],hint);
#endif


namespace gms {

          namespace math {


              struct AVXQuatArray1D __ATTR_ALIGN__(32) {

                     __ATTR_ALIGN__(8) double * __restrict m_x; // x - part
		     __ATTR_ALIGN__(8) double * __restrict m_y; // y - part
		     __ATTR_ALIGN__(8) double * __restrict m_z; // z - part
		     __ATTR_ALIGN__(8) double * __restrict m_w; // w - part
		     int32_t m_size; // size of quaternion array 1D.
#if(USE_STRUCT_PADDING) == 1
                     PAD_TO_ALIGNED(8,0,28)
#endif
                     /*
                          Default Ctor
                      */
		     __ATTR_ALWAYS_INLINE__
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
                     AVXQuatArray1D() {

                         m_x = nullptr;
			 m_y = nullptr;
			 m_z = nullptr;
			 m_w = nullptr;
			 m_size = 0;
		     }

		     __ATTR_ALWAYS_INLINE__
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     AVXQuatArray1D(const int32_t size,
		                    const bool lockm) {
                        using namespace gms::common;
			m_size = size;
			const std::size_t ullsize = static_cast<m_size>;
#if (USE_MMAP_2MiB) == 1
#pragma omp parallel sections
{
                #pragma omp section
		{
                        m_x = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
		}
		#pragma omp section
		{
			m_y = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
		}
		#pragma omp section
		{
			m_z = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
		}
		#pragma omp section
		{
			m_w = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
		}
}
#else
#pragma omp parallel sections
{
                #pragma omp section
		{
                        m_x = gms_dmalloca(ullsize,align32B);
		}
		#pragma omp section
		{
			m_y = gms_dmalloca(ullsize,align32B);
		}
		#pragma omp section
			m_z = gms_dmalloca(ullsize,align32B);
		}
		#pragma omp section
		{
			m_w = gms_dmalloca(ullsize,align32B);
		}
}
#endif
                        // Zero initialization.
#if (ZERO_INIT_ARRAYS) == 1

             	      avx256_init_unroll8x_pd(&m_x[0],m_size,0.0);
	      	      avx256_init_unroll8x_pd(&m_y[0],m_size,0.0);
	      	      avx256_init_unroll8x_pd(&m_z[0],m_size,0.0);
	      	      avx256_init_unroll8x_pd(&m_w[0],m_size,0.0);

#endif
		 }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     AVXQuatArray1D(const double * __restrict __ATTR_ALIGN__(32) x,
		                    const double * __restrict __ATTR_ALIGN__(32) y,
				    const double * __restrict __ATTR_ALIGN__(32) z,
				    const double * __restrict __ATTR_ALIGN__(32) w,
				    const int32_t size,
				    const bool lockm) {
                          using namespace gms::common;
			  m_size = size;
                          const std::size_t ullsize = static_cast<m_size>;
#if (USE_MMAP_2MiB) == 1
                        m_x = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
			m_y = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
			m_z = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
			m_w = gms_edmmap_2MiB(ullsize,
			                      PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
					      -1,0,lockm);
#else
                        m_x = gms_dmalloca(ullsize,align32B);
			m_y = gms_dmalloca(ullsize,align32B);
			m_z = gms_dmalloca(ullsize,align32B);
			m_w = gms_dmalloca(ullsize,align32B);
#endif
#if (USE_AVXQUATERNION_ARRAY1D_NT_STORES) == 1
                        avx256_uncached_memmove(&m_x[0],&x[0],m_size);
			avx256_uncached_memmove(&m_y[0],&y[0],m_size);
			avx256_uncached_memmove(&m_z[0],&z[0],m_size);
			avx256_uncached_memmove(&m_w[0],&w[0],m_size);
#else
                        avx256_cached_memmove(&m_x[0],&x[0],m_size);
			avx256_cached_memmove(&m_y[0],&y[0],m_size);
			avx256_cached_memmove(&m_z[0],&z[0],m_size);
			avx256_cached_memmove(&m_w[0],&w[0],m_size);
#endif
		    }


		    AVXQuatArray1D(const AVXQuatArray1D &x) = delete;

		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    AVXQuatArray1D(AVXQuatArray1D &&rhs) {

                        m_x = &rhs.m_x[0];
			m_y = &rhs.m_y[0];
			m_z = &rhs.m_z[0];
			m_w = &rhs.m_w[0];
			m_size = rhs.m_size;
			rhs.m_x = NULL;
			rhs.m_y = NULL;
			rhs.m_z = NULL;
			rhs.m_w = NULL;
			rhs.m_size = 0;
		    }


		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    ~AVXQuatArray1D() {
                       if(NULL != m_x) _mm_free(m_x); m_x = NULL;
		       if(NULL != m_y) _mm_free(m_y); m_y = NULL;
		       if(NULL != m_z) _mm_free(m_z); m_z = NULL;
		       if(NULL != m_w) _mm_free(m_w); m_w = NULL;
                       
		    }

		    AVXQuatArray1D & operator=(const AVXQuatArray1D &rhs) = delete;

		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    AVXQuatArray1D &
		    operator=(AVXQuatArray1D &&rhs) {
                        if(this == &rhs) return (*this);
			_mm_free(m_x);
			_mm_free(m_y);
			_mm_free(m_z);
			_mm_free(m_w);
			m_x = &rhs.m_x[0];
			rhs.m_x = NULL;
			m_y = &rhs.m_y[0];
			rhs.m_y = NULL;
			m_z = &rhs.m_z[0];
			rhs.m_z = NULL;
			m_w = &rhs.m_w[0];
			rhs.m_w = NULL;
			m_size = rhs.m_size;
			rhs.m_size = 0;
			return (*this);
		    }

		     
		      
	      }; // AVXQuatArray1D


	       __ATTR_ALWAYS_INLINE__
	       __ATTR_HOT__
	       __ATTR_ALIGN__(32)
	       static inline
	       void AVXQuat_Add_AVXQuat(AVXQuatArray1D &c,
	                                const AVXQuatArray1D &b,
					const AVXQuatArray1D &a) {
		      if(b.m_size != a.m_size) { return;}
                      int32_t i;
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma code_align(32)
#endif		      
                      for(i = 0; i != ROUND_TO_FOUR(c.m_size,4); i += 16) {
		      
                          const __m256d ymm0(_mm256_load_pd(&b.m_x[i+0]));
			  const __m256d ymm1(_mm256_load_pd(&a.m_x[i+0]));
			  _mm256_store_pd(&c.m_x[i+0],
			                  _mm256_add_pd(ymm0,ymm1));
			  const __m256d ymm2(_mm256_load_pd(&b.m_y[i+0]));
			  const __m256d ymm3(_mm256_load_pd(&a.m_y[i+0]));
			  _mm256_store_pd(&c.m_y[i+0],
			                  _mm256_add_pd(ymm2,ymm3));
			  const __m256d ymm4(_mm256_load_pd(&b.m_z[i+0]));
			  const __m256d ymm5(_mm256_load_pd(&a.m_z[i+0]));
			  _mm256_store_pd(&c.m_z[i+0],
			                  _mm256_add_pd(ymm4,ymm5));
			  const __m256d ymm6(_mm256_load_pd(&b.m_w[i+0]));
			  const __m256d ymm7(_mm256_load_pd(&a.m_w[i+0]));
			  _mm256_store_pd(&c.m_w[i+0],
			                  _mm256_add_pd(ymm6,ymm7));
			  const __m256d ymm8(_mm256_load_pd(&b.m_x[i+4]));
			  const __m256d ymm9(_mm256_load_pd(&a.m_x[i+4]));
			  _mm256_store_pd(&c.m_x[i+4],
			                  _mm256_add_pd(ymm8,ymm9));
			  const __m256d ymm10(_mm256_load_pd(&b.m_y[i+4]));
			  const __m256d ymm11(_mm256_load_pd(&a.m_y[i+4]));
			  _mm256_store_pd(&c.m_y[i+4],
			                  _mm256_add_pd(ymm10,ymm11));
			  const __m256d ymm12(_mm256_load_pd(&b.m_z[i+4]));
			  const __m256d ymm13(_mm256_load_pd(&a.m_z[i+4]));
			  _mm256_store_pd(&c.m_z[i+4],
			                  _mm256_add_pd(ymm12,ymm13));
			  const __m256d ymm14(_mm256_load_pd(&b.m_w[i+4]));
			  const __m256d ymm15(_mm256_load_pd(&a.m_w[i+4]));
			  _mm256_store_pd(&c.m_w[i+4],
			                  _mm256_add_pd(ymm14,ymm15));
			  const __m256d ymm16(_mm256_load_pd(&b.m_x[i+8]));
			  const __m256d ymm17(_mm256_load_pd(&a.m_x[i+8]));
			  _mm256_store_pd(&c.m_x[i+8],
			                  _mm256_add_pd(ymm16,ymm17));
			  const __m256d ymm18(_mm256_load_pd(&b.m_y[i+8]));
			  const __m256d ymm19(_mm256_load_pd(&a.m_y[i+8]));
			  _mm256_store_pd(&c.m_y[i+8],
			                  _mm256_add_pd(ymm18,ymm19));
			  const __m256d ymm20(_mm256_load_pd(&b.m_z[i+8]));
			  const __m256d ymm21(_mm256_load_pd(&a.m_z[i+8]));
			  _mm256_store_pd(&c.m_z[i+8],
			                  _mm256_add_pd(ymm20,ymm21));
			  const __m256d ymm22(_mm256_load_pd(&b.m_w[i+8]));
			  const __m256d ymm23(_mm256_load_pd(&a.m_w[i+8]));
			  _mm256_store_pd(&c.m_w[i+8],
			                  _mm256_add_pd(ymm22,ymm23));
			  const __m256d ymm24(_mm256_load_pd(&b.m_x[i+12]));
			  const __m256d ymm25(_mm256_load_pd(&a.m_x[i+12]));
			  _mm256_store_pd(&c.m_x[i+12],
			                  _mm256_add_pd(ymm24,ymm25));
			  const __m256d ymm26(_mm256_load_pd(&b.m_y[i+12]));
			  const __m256d ymm27(_mm256_load_pd(&a.m_y[i+12]));
			  _mm256_store_pd(&c.m_y[i+12],
			                  _mm256_add_pd(ymm26,ymm27));
			  const __m256d ymm28(_mm256_load_pd(&b.m_z[i+12]));
			  const __m256d ymm29(_mm256_load_pd(&a.m_z[i+12]));
			  _mm256_store_pd(&c.m_z[i+12],
			                  _mm256_add_pd(ymm28,ymm29));
			  const __m256d ymm30(_mm256_load_pd(&b.m_w[i+12]));
			  const __m256d ymm31(_mm256_load_pd(&a.m_w[i+12]));
			  _mm256_store_pd(&c.m_w[i+12],
			                  _mm256_add_pd(ymm30,ymm31));
			  
		      }
#if defined (__ICC) || defined (__INTEL_COMPILER)
#pragma loop_count min(1),avg(2),max(3)
#endif
                      for(; i != c.m_size; ++i) {
                           c.m_x[i] = b.m_x[i]+a.m_x[i];
			   c.m_y[i] = b.m_y[i]+a.m_y[i];
			   c.m_z[i] = b.m_z[i]+a.m_z[i];
			   c.m_w[i] = b.m_w[i]+a.m_w[i];
		      }
	       }
	      
	  
     } // math


}// gms






#endif /*__GMS_AVX_QUATERNION_ARRAY1D__HPP__*/
