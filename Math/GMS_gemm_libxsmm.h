

#ifndef __GMS_GEMM_LIBXSMM_H__
#define __GMS_GEMM_LIBXSMM_H__


#include <libxsmm.h>
#include <cstdint>
#include "GSM_config.h"


namespace math {

       namespace math {

                __ATTR_ALWAYS_INLINE__
		inline
		static
		void sgemm_libxsmm_NN(const float * __restrict __ATTR_ALIGN__(64) a,
		                      const float * __restrict __ATTR_ALIGN__(64) b,
				      float * __restrict __ATTR_ALIGN__(64) c,
				      const int32_t M,
				      const int32_t K,
				      const int32_t N) {
 
                  const float alpha = 1.0f;
		  const float beta  = 0.0f;
		  const char transa = 'N';
		  const char transb = 'N';
		  libxsmm_sgemm(&transa,&transb,
		                N,M,K,&alpha,
				b,&N,a,&K,&beta,
				c,&N);

	     }


	        __ATTR_ALWAYS_INLINE__
		inline
		static
		void dgemm_libxsmm_NN(const double * __restrict __ATTR_ALIGN__(64) a,
		                      const double * __restrict __ATTR_ALIGN__(64) b,
				      double * __restrict __ATTR_ALIGN__(64) c,
				      const int32_t M,
				      const int32_t K,
				      const int32_t N) {
 
                  const double alpha = 1.0;
		  const double beta  = 0.0;
		  const char transa = 'N';
		  const char transb = 'N';
		  libxsmm_dgemm(&transa,&transb,
		                N,M,K,&alpha,
				b,&N,a,&K,&beta,
				c,&N);

	     }


	       __ATTR_ALWAYS_INLINE__
		inline
		static
		void dgemm_libxsmm_NY(const double * __restrict __ATTR_ALIGN__(64) a,
		                      const double * __restrict __ATTR_ALIGN__(64) b,
				      double * __restrict __ATTR_ALIGN__(64) c,
				      const int32_t M,
				      const int32_t K,
				      const int32_t N) {
 
                  const double alpha = 1.0;
		  const double beta  = 0.0;
		  const char transa = 'N';
		  const char transb = 'Y';
		  libxsmm_dgemm(&transa,&transb,
		                N,M,K,&alpha,
				b,&N,a,&M,&beta,
				c,&N);

	     }


	        __ATTR_ALWAYS_INLINE__
		inline
		static
		void dgemm_libxsmm_YN(const double * __restrict __ATTR_ALIGN__(64) a,
		                      const double * __restrict __ATTR_ALIGN__(64) b,
				      double * __restrict __ATTR_ALIGN__(64) c,
				      const int32_t M,
				      const int32_t K,
				      const int32_t N) {
 
                  const double alpha = 1.0;
		  const double beta  = 0.0;
		  const char transa = 'Y';
		  const char transb = 'N';
		  libxsmm_dgemm(&transa,&transb,
		                N,M,K,&alpha,
				b,&K,a,&K,&beta,
				c,&N);

	     }


	        __ATTR_ALWAYS_INLINE__
		inline
		static
		void dgemm_libxsmm_YY(const double * __restrict __ATTR_ALIGN__(64) a,
		                      const double * __restrict __ATTR_ALIGN__(64) b,
				      double * __restrict __ATTR_ALIGN__(64) c,
				      const int32_t M,
				      const int32_t K,
				      const int32_t N) {
 
                  const double alpha = 1.0;
		  const double beta  = 0.0;
		  const char transa = 'Y';
		  const char transb = 'Y';
		  libxsmm_dgemm(&transa,&transb,
		                N,M,K,&alpha,
				b,&K,a,&M,&beta,
				c,&N);

	     }


	     




	     
   }


}















#endif /*__GMS_GEMM_LIBXSMM_H__*/
