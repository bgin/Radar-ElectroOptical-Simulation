
#ifndef __GMS_SMALL_CVECTORS_HPP__
#define __GMS_SMALL_CVECTORS_HPP__

#if !defined (LAM_SMALL_CVECTORS_MAJOR)
#define LAM_SMALL_CVECTORS_MAJOR 1
#endif

#if !defined (LAM_SMALL_CVECTORS_MINOR)
#define LAM_SMALL_CVECTORS_MINOR 0
#endif

#if !defined (LAM_SMALL_CVECTORS_MICRO)
#define LAM_SMALL_CVECTORS_MICRO 0
#endif

#if !defined (LAM_SMALL_CVECTORS_FULLVER)
#define LAM_SMALL_VECTORS_FULLVER 1000
#endif

#if !defined (LAM_SMALL_CVECTORS_CREATE_DATE)
#define LAM_SMALL_VECTORS_CREATE_DATE "23-08-2017 15:28 +00200 (WED 23 AUG 2017 GMT+2)"
#endif

/*
	Set this value to latest build date/time
*/
#if !defined (LAM_SMALL_CVECTORS_BUILD_DATE)
#define LAM_SMALL_CVECTORS_BUILD_DATE " "
#endif

#if !defined (LAM_SMALL_CVECTORS_AUTHOR)
#define LAM_SMALL_CVECTORS_AUTHOR "Programmer: Bernard Gingold e-mail: beniekg@gmail.com"
#endif

#if !defined (LAM_SMALL_CVECTORS_DESCRIPT)
#define LAM_SMALL_CVECTORS_DESCRIPT "Compile-time complex short vectors and their basic arithmetic operators."
#endif

namespace file_info {

#if defined _WIN64  
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

   const unsigned int gGMS_SMALL_CVECTORS_MAJOR =  gms::common::gVersionInfo.m_VersionMajor;
   const unsigned int gGMS_SMALL_CVECTORS_MINOR =  gms::common::gVersionInfo.m_VersionMinor;
   const unsigned int gGMS_SMALL_CVECTORS_MICRO =  gms::common::gVersionInfo.m_VersionMicro;
   const unsigned int gGMS_SMALL_CVECTORS_FULLVER =
     1000U*gGMS_SMALL_CVECTORS_MAJOR+100U*gGMS_SMALL_CVECTORS_MINOR+10U*gGMS_SMALL_CVECTORS_MICRO;
  const char * const pgGMS_SMALL_CVECTORS_CREATE_DATE = "23-08-2017 15:28 +00200 (WED 23 AUG 2017 GMT+2)";
  const char * const pgGMS_SMALL_CVECTORS_BUILD_DATE  = "00-00-0000 00:00";
  const char * const pgGMS_SMALL_CVECTORS_AUTHOR      = "Programmer: Bernard Gingold contact: beniekg@gmail.com";
  const char * const pgGMS_SMALL_CVECTORS_SYNOPSIS    = "Compile-time complex short vectors and their basic arithmetic operators.";
}


#include <cstdint>
#include <type_traits>
#if defined _WIN64
    #include "../GMS_config.h"
    #include "../GMS_simd_defs.h"
#elif defined __linux
    #include "GMS_config.h"
    #include "GMS_simd_defs.h"
#endif

namespace gms  {
	namespace math {


		/*

		     Meta Vectors for small-sized complex-valued vectors.

		*/

		/*
		    @brief   Complex addition of two vectors
					 oRe = iRe1+iRe2
					 oIm = iIm1+iIm2
			@Result:
					 _Inout_ vectors: oRe,oIm
		*/
		template<typename Real_t, uint64_t N> struct CVecAdd1 {

			 
			void  operator()(Real_t* __restrict oRe,
					 Real_t* __restrict oIm,
				         const Real_t* __restrict iRe1,
				         const Real_t* __restrict iIm1,
				         const Real_t* __restrict iRe2,
				         const Real_t* __restrict iIm2) {

						 *oRe = *iRe1 + *iRe2;
						 *oIm = *iIm1 + *iIm2;
						 CVecAdd1<Real_t,N-1ULL> vcadd;
						 vcadd.operator()(oRe+1,oIm+1,iRe1+1,iIm1+1,iRe2+1,iIm2+1);
						 
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecAdd1<Real_t, 0ULL> {

			
			void operator()(Real_t* __restrict oRe,
					Real_t* __restrict oIm,
			                const Real_t* __restrict iRe1,
			                const Real_t* __restrict iIm1,
			                const Real_t* __restrict iRe2,
			                const Real_t* __restrict iIm2){

			}
		};

		/*
			@brief    Complex addition
					  oRe = iRe1+iRe2
		    @Result:  
					  _Inout_ vector: oRe

		*/
		template<typename Real_t, uint64_t N> struct CVecAdd2 {

			
			void operator()(Real_t* __restrict oRe,
				         const Real_t* __restrict iRe1,
				         const Real_t* __restrict iRe2) {
				
				*oRe = *iRe1 + *iRe2;
				CVecAdd2<Real_t, N-1ULL> vcadd2;
				vcadd2.operator()(oRe+1,iRe1+1,iRe2+1);
              
			}
		};

		/*
		@	brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecAdd2<Real_t, 0ULL> {
				
			
			void operator()(Real_t* __restrict oRe,
					  const Real_t* __restrict iRe1,
				        const Real_t* __restrict iRe2) {

			 }
		};

		/*
			@brief   Complex subtraction of two vectors
				     oRe = iRe1-iRe2
		             oIm = iIm1-iIm2
		    @Result:
		            _Inout_ vectors: oRe,oIm
		*/
		template<typename Real_t, uint64_t N> struct CVecSub1 {

			
			void operator()(Real_t* __restrict oRe,
				       Real_t* __restrict oIm,
			          const Real_t* __restrict iRe1,
			           const Real_t* __restrict iRe2,
			           const Real_t* __restrict iIm1,
			           const Real_t* __restrict iIm2) {
				
				
				*oRe = *iRe1 - *iRe2;
				*oIm = *iIm1 - *iIm2;
				CVecSub1<Real_t,N-1ULL> cvsub;
				cvsub.operator()(oRe+1,oIm+1,iRe1+1,iRe2+1,iIm1+1,iIm2+1);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecSub1<Real_t, 0ULL> {

			
			void operator()(Real_t* __restrict oRe,
					   Real_t* __restrict oIm,
			           const Real_t* __restrict iRe1,
			           const Real_t* __restrict iRe2,
			           const Real_t* __restrict iIm1,
			           const Real_t* __restrict iIm2) {

			}
		};

		/*
			@brief    Complex subtraction
		              oRe = iRe1-iRe2
		    @Result:
		             _Inout_ vector: oRe

		*/
		template<typename Real_t, uint64_t N> struct CVecSub2{

			
			void operator()(Real_t* __restrict oRe,
				         const Real_t* __restrict iRe1,
					  const Real_t* __restrict iRe2) {
					
					*oRe = *iRe1-*iRe2;
					CVecSub2<Real_t,N-1ULL> cvsub2;
					cvsub2.operator()(oRe+1,iRe1+1,iRe2+1);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecSub2<Real_t, 0ULL>{

		
			void operator()(Real_t* __restrict oRe,
					  const Real_t* __restrict oRe1,
			            const Real_t* __restrict oRe2) {
 
			}
		};

		/*
			@brief   Complex multiplication of two vectors
					 oRe = iRe1*iRe2
		             oIm = iIm1*iIm2
		    @Result:
		            _Inout_ vectors: oRe,oIm
		*/
		template<typename Real_t, uint64_t N> struct CVecMul1{

			
		 void	operator()(Real_t* __restrict oRe,
				     Real_t* __restrict oIm,
			           const Real_t* __restrict iRe1,
			           const Real_t* __restrict iRe2,
			          const Real_t* __restrict iIm1,
			          const Real_t* __restrict iIm2) {

				*oRe = (*iRe1 * *iRe2) - (*iIm1 * *iIm2);
				*oIm = (*iIm1 * *iRe2) + (*iRe1 * *iIm2);
				CVecMul1<Real_t,N-1ULL> cvmul;
				cvmul.operator()(oRe+1,oIm+1,iRe1+1,iRe2+1,iIm1+1,iIm2+1);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecMul1<Real_t, 0ULL> {

			
			void operator()(Real_t* __restrict oRe,
				         Real_t* __restrict oIm,
					  const Real_t* __restrict iRe1,
			          const Real_t* __restrict iRe2,
			           const Real_t* __restrict iIm1,
			          const Real_t* __restrict iIm2) {

			}
		};

		/*
		    @brief    Complex multiplication
		              oRe = iRe1*iRe
					  oIm = iIm1*iRe
		    @Result:
				     _Inout_ vector: oRe,oIm

		*/
		template<typename Real_t, uint64_t N> struct CVecMul2{

			
		   void	operator()( Real_t* __restrict oRe,
			         Real_t* __restrict oIm,
			         const Real_t* __restrict iRe,
			         const Real_t* __restrict iIm,
			       	   const Real_t* __restrict v) {
				
				*oRe = *iRe * *v;
				*oIm = *iIm * *v;
				CVecMul2<Real_t,N-1ULL> cvecmul;
				cvecmul.operator()(oRe+1,oIm+1,iRe+1,iIm+1,v+1);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecMul2<Real_t, 0ULL> {

			
			void operator()( Real_t* __restrict oRe,
					  Real_t* __restrict oIm,
					  const Real_t* __restrict iRe,
					  const Real_t* __restrict iIm,
					  const Real_t* __restrict v) {

			}
		};

		/*
			@brief   Complex division of two vectors
				     oRe = iRe1/iRe2
				     oIm = iIm1/iIm2
		    @Result:
		            _Inout_ vectors: oRe,oIm
		*/
		template<typename Real_t, uint64_t N> struct CVecDiv1{

			
		  void	operator()(Real_t* __restrict oRe,
				   Real_t* __restrict oIm,
			            const Real_t* __restrict iRe1,
			            const Real_t* __restrict iRe2,
			         const Real_t* __restrict iIm1,
			         const Real_t* __restrict iIm2) {

				 Real_t tr,ti,tmp;
				 tr = (*iRe1 * *iRe2) + (*iIm1 * *iIm2);
				 ti = (*iIm1 * *iRe2) - (*iRe1 * *iRe2);
				 tmp = (*iRe2 * *iRe2) + (*iIm2 * *iIm2);
				 *oRe = tr/tmp;
				 *oIm = tr/tmp;
				 CVecDiv1<Real_t,N-1ULL> cvdiv;
				 cvdiv.operator()(oRe+1,oIm+1,iRe1+1,iRe2+1,iIm1+1,iIm2+1);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecDiv1<Real_t, 0ULL> {

			
		  void	operator()(Real_t* __restrict oRe,
			            Real_t* __restrict oIm,
			            const Real_t* __restrict iRe1,
			            const Real_t* __restrict iRe2,
			           const Real_t* __restrict iIm1,
			            const Real_t* __restrict iIm2) {

			}
		};

		/*
			@brief    Complex multiplication
					  oRe = iRe1*iRe
				      oIm = iIm1*iRe
		    @Result:
		             _Inout_ vector: oRe,oIm

		*/
		template<typename Real_t, uint64_t N> struct CVecDiv2{

			
		 void	operator()(Real_t* __restrict oRe,
			           Real_t* __restrict oIm,
			           const Real_t* __restrict iRe,
			            const Real_t* __restrict iIm,
			           const Real_t* __restrict v) {
				
				*oRe = *iRe / *v;
				*oIm = *iIm / *v;
				CVecDiv2<Real_t,N-1ULL> cvecdiv;
				cvecdiv.operator()(oRe+1,oIm+1,iRe+1,iIm+1,v+1);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<typename Real_t> struct CVecDiv2<Real_t, 0ULL>{

			
			void operator()(Real_t* __restrict oRe,
				         Real_t* __restrict oIm,
			          const Real_t* __restrict iRe,
			           const Real_t* __restrict iIm,
			          const Real_t* __restrict v){

			}
		};


		//                        //
		//  SIMD Vectorization.   //
		//  AVX based versions.  //
		//                       //

		
		template<typename uint64_t N> struct AVXCVecAdd{
			
			static_assert((N % 4ULL) == 0ULL, "N % 4 != 0ULL");
			void operator()(double* __restrict oRe,
				         double* __restrict oIm,
				         const double* __restrict iRe1,
				          const double* __restrict iRe2,
				          const double* __restrict iIm1,
				           const double* __restrict iIm2) {

#if GMS_CACHE_MEM_STORES == 1

				_mm256_storeu_pd(oRe, _mm256_add_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iRe2)));
				_mm256_storeu_pd(oIm, _mm256_add_pd(_mm256_loadu_pd(iIm1), _mm256_loadu_pd(iIm2)));
#else
				_mm256_stream_pd(oRe, _mm256_add_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iRe2)));
				_mm256_stream_pd(oIm, _mm256_add_pd(_mm256_loadu_pd(iIm1), _mm256_loadu_pd(iIm2)));
#endif

				AVXCVecAdd<N-4ULL> avec;
				avec.operator()(oRe+4,oIm+4,iRe1+4,iRe2+4,iIm1+4,iIm2+4);
			}

		};
			
		/*
			@brief     specialization for the terminating condition.
		*/
		template<> struct AVXCVecAdd<0ULL>{

			void operator()(double* __restrict oRe,
				         double* __restrict oIm,
				         const double* __restrict iRe1,
				         const double* __restrict iRe2,
				         const double* __restrict iIm1,
				         const double* __restrict iIm2){

			}
		};

		template<uint64_t N> struct AVXCVecAdd2{

			static_assert((N % 4ULL) == 0, "N % 4 != 0");
			void operator()( double* __restrict oRe,
				         const double* __restrict oRe1,
				         const double* __restrict oRe2) {

#if GMS_CACHE_MEM_STORES == 1
				_mm256_storeu_pd(oRe, _mm256_add_pd(_mm256_loadu_pd(oRe1), _mm256_loadu_pd(oRe2)));
#else
				_mm256_storeu_pd(oRe, _mm256_add_pd(_mm256_loadu_pd(oRe1), _mm256_loadu_pd(oRe2)));
#endif
				AVXCVecAdd2<N-4ULL> avec;
				avec.operator()(oRe+4,oRe1+4,oRe2+4);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/

		template<> struct AVXCVecAdd2<0ULL>{

			void operator()(double* __restrict oRe,
				         const double* __restrict oRe1,
				         const double* __restrict oRe2) {

			}
		};

		template<uint64_t N> struct AVXCVecSub{
			
			static_assert((N % 4ULL) == 0, "N % 4 != 0");
			void operator()(double* __restrict oRe,
				         double* __restrict oIm,
				         const double* __restrict iRe1,
				           const double* __restrict iRe2,
				           const double* __restrict iIm1,
				          const double* __restrict iIm2) {

#if GMS_CACHE_MEM_STORES == 1
				_mm256_storeu_pd(oRe, _mm256_sub_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iRe2)));
				_mm256_storeu_pd(oIm, _mm256_sub_pd(_mm256_loadu_pd(iIm1), _mm256_loadu_pd(iIm2)));
#else
				_mm256_stream_pd(oRe, _mm256_sub_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iRe2)));
				_mm256_stream_pd(oIm, _mm256_sub_pd(_mm256_loadu_pd(iIm1), _mm256_loadu_pd(iIm2)));
#endif
				AVXCVecSub<N-4ULL> avsub;
				avsub.operator()(oRe+4,oIm+4,iRe1+4,iRe2+4,iIm1+4,iIm2+4);

			}

		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<> struct AVXCVecSub<0ULL> {

			void operator()(double* __restrict oRe,
				         double* __restrict oIm,
				         const double* __restrict iRe1,
				         const double* __restrict iRe2,
				         const double* __restrict iIm1,
				         const double* __restrict iIm2){


			}
		};

		template<uint64_t N> struct AVXCVecSub2{

			static_assert((N % 4ULL) == 0, "N % 4 != 0");
			void operator()(double* __restrict oRe,
				         double* __restrict iRe1,
				         double* __restrict iRe2) {

#if GMS_CACHE_MEM_STORES == 1
				_mm256_storeu_pd(oRe, _mm256_sub_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iRe2)));
#else
				_mm256_stream_pd(oRe, _mm256_sub_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iRe2)));
#endif
				AVXCVecSub2<N-4ULL> avsub;
				avsub.operator()(oRe+4,iRe1+4,iRe2+4);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<> struct AVXCVecSub2<0ULL> {

			void operator()(double* __restrict oRe,
				         double* __restrict iRe1,
				         double* __restrict iRe2) {

			}
		};

		template<uint64_t N> struct AVXCVecMul{

			static_assert((N % 4ULL) == 0, "N % 4 != 0");
			void operator()(double* __restrict oRe,
				         double* __restrict oIm,
				         const double* __restrict iRe1,
				         const double* __restrict iRe2,
				         const double* __restrict iIm1,
				         const double* __restrict iIm2) {

#if LAM_CACHE_MEM_STORES == 1
				_mm256_storeu_pd(oRe,_mm256_sub_pd(_mm256_mul_pd(_mm256_loadu_pd(iRe1)),
					_mm256_loadu_pd(iRe2),_mm256_mul_pd(_mm256_loadu_pd(iIm1),_mm256_loadu_pd(iIm2))));
				_mm256_storeu_pd(oIm,_mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(iIm1),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iIm2))));
#else
				_mm256_stream_pd(oRe, _mm256_sub_pd(_mm256_mul_pd(_mm256_loadu_pd(iRe1)),
					_mm256_loadu_pd(iRe2), _mm256_mul_pd(_mm256_loadu_pd(iIm1), _mm256_loadu_pd(iIm2))));
				_mm256_stream_pd(oIm, _mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(iIm1),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iIm2))));
#endif
				AVXCVecMul<N-4ULL> avmul;
				avmul.operator()(oRe+4,oIm+4,iRe1+4,iRe2+4,iIm1+4,iIm2+4);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<> struct AVXCVecMul<0ULL>{

			void operator()(double* __restrict oRe,
					         double* __restrict oIm,
					     const double* __restrict iRe1,
				         const double* __restrict iRe2,
				         const double* __restrict iIm1,
				         const double* __restrict iIm2) {

			}
		};

		template<uint64_t N> struct AVXCVecMul2{

			static_assert((N % 4ULL) == 0, "N % 4 != 0");
			void operator()(double* __restrict oRe,
					         double* __restrict oIm,
					         const double* __restrict iRe,
					         const double* __restrict iIm,
					         const double* __restrict v) {

#if GMS_CACHE_MEM_STORES == 1
				_mm256_storeu_pd(oRe, _mm256_mul_pd(_mm256_loadu_pd(iRe), _mm256_loadu_pd(v)));
				_mm256_storeu_pd(oIm, _mm256_mul_pd(_mm256_loadu_pd(iIm), _mm256_loadu_pd(v)));
#else
				_mm256_stream_pd(oRe, _mm256_mul_pd(_mm256_loadu_pd(iRe), _mm256_loadu_pd(v)));
				_mm256_stream_pd(oIm, _mm256_mul_pd(_mm256_loadu_pd(iIm), _mm256_loadu_pd(v)));
#endif
				AVXCVecMul2<N-4ULL> avmul;
				avmul.operator()(oRe+4,oIm+4,iRe+4,iIm+4,v+4);
			}
		};


		/*
			@brief     specialization for the terminating condition.
		*/
		template<> struct AVXCVecMul2<0ULL> {

			void operator()(double* __restrict oRe,
				         double* __restrict oIm,
				        const double* __restrict iRe,
				         const double* __restrict iIm,
				         const double* __restrict v) {

			}
		};

		template<uint64_t N> struct AVXCVecDiv{

			static_assert((N % 4ULL) == 0, "N % 4 != 0");
			void operator()(double* __restrict oRe,
					         double* __restrict oIm,
					         const double* __restrict iRe1,
					         const double* __restrict iRe2,
					         const double* __restrict iIm1,
					         const double* __restrict iIm2) {

#if GMS_CACHE_MEM_STORES == 1
				__m256d tr(_mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(iRe1),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iIm1),_mm256_loadu_pd(iIm2))));
				__m256d ti(_mm256_sub_pd(_mm256_mul_pd(_mm256_loadu_pd(iIm1),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iRe1),_mm256_loadu_pd(iIm2))));
				__m256d tmp(_mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(iRe2),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iIm2), _mm256_loadu_pd(iIm2))));
				_mm256_storeu_pd(oRe, _mm256_div_pd(tr,tmp));
				_mm256_storeu_pd(oIm, _mm256_div_pd(ti,tmp));

#else
				__m256d tr(_mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(iRe1),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iIm1), _mm256_loadu_pd(iIm2))));
				__m256d ti(_mm256_sub_pd(_mm256_mul_pd(_mm256_loadu_pd(iIm1),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iRe1), _mm256_loadu_pd(iIm2))));
				__m256d tmp(_mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(iRe2),
					_mm256_loadu_pd(iRe2)), _mm256_mul_pd(_mm256_loadu_pd(iIm2), _mm256_loadu_pd(iIm2))));
				_mm256_stream_pd(oRe, _mm256_div_pd(tr, tmp));
				_mm256_stream_pd(oIm, _mm256_div_pd(ti, tmp));
#endif
				AVXCVecDiv<N-4ULL> axvec;
				axvec.operator()(oRe+4,oIm+4,iRe1+4,iRe2+4,iIm1+4,iIm2+4);
			}
		};

		/*
			@brief     specialization for the terminating condition.
		*/
		template<> struct AVXCVecDiv<0ULL>{

			void operator()( double* __restrict oRe,
					         double* __restrict oIm,
						         const double* __restrict iRe1,
							const double* __restrict iRe2,
							 const double* __restrict iIm1,
							 const double* __restrict iIm2) {

			}
		};
	}
}


#endif /*__GMS_SMALL_CVECTORS_HPP__*/
