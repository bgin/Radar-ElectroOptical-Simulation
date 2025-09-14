
#ifndef __GMS_AVXVCPROCEDURES_H__
#define __GMS_AVXVCPROCEDURES_H__ 061020191545





namespace file_info {
  

	
	const unsigned int gGMS_AVXCPROCEDURES_MAJOR = 1;

	const unsigned int gGMS_AVXCPROCEDURES_MINOR = 0;

	const unsigned int gGMS_AVXCPROCEDURES_MICRO = 1;

	const unsigned int gGMS_AVXCPROCEDURES_FULLVER = 
		1000U*gLAM_AVXCOMPLEX_SMALLV_MAJOR+100U*gLAM_AVXCOMPLEX_SMALLV_MINOR+10U*gLAM_AVXCOMPLEX_SMALLV_MICRO;

	const char * const pgGMS_AVXCPROCEDURES_CREATE_DATE = "06-10-2019 15:45 + 00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_AVXCPROCEDURES_BUILD_DATE = __DATE__ ":"__TIME__;

	const char * const pgGMS_AVXCPROCEDURES_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVXCPROCEDURES_SYNOPSIS = "AVX complex vector (1D) stack-allocated storage.";
}

#include "GMS_avxcomplex.h"

namespace gms {
	namespace math {

		//
		//	Collection of optimized void procedures operating on
		//  objects of type AVXVComplex1D.
		//

		/*
		@Purpose:
					Mimicks this operation -- C = A+B (vectors size obviously must be equal to each other).
					where A is of type: AVXVComplex1D
					and   B is of type: AVXVComplex1D
		@Warning:	
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented.
		*/
		void
		avxvcomplex_add(AVXVComplex1D &,const AVXVComplex1D & , 
		                const AVXVComplex1D & )      noexcept(true);

		/*
		@Purpose:
					Mimicks this operation -- C = A+B (vectors size obviously must be equal to each other).
					where A is of type: AVXVComplex1D
					and B is of type: double * __restrict
		@Remark:
					Size of double * argument must be
					equal to size of AVXVComplex1D.m_Re
		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void
		avxvcomplex_add(AVXVComplex1D &,const AVXVComplex1D &, 
			            const double * __restrict)    noexcept(true);

		/*
		@Purpose:
					Mimicks this operation -- C = A-B (vectors size obviously must be equal to each other).
					where A is of type: AVXVComplex1D
					and B is of type:   AVXVComplex1D

		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void
		avxvcomplex_sub(AVXVComplex1D &, const AVXVComplex1D &, 
						const AVXVComplex1D &)  noexcept(true);

		/*
		@Purpose:
					Mimicks this operation -- C = A-B (vectors size obviously must be equal to each other).
					where A is of type: AVXVComplex1D
					and B is of type: double * __restrict
		@Remark:
					Size of double * argument must be
					equal to size of AVXVComplex1D.m_Re

		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void
		avxvcomplex_sub(AVXVComplex1D &,const AVXVComplex1D &, 
					   const double * __restrict)   noexcept(true);

		/*
		@Purpose:
					Mimicks this operation -- C = A*B (vectors size obviously must be equal to each other).
					where A is of type: AVXVComplex1D
					and   B is of type: AVXVComplex1D

		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void avxvcomplex_mul(AVXVComplex1D &, const AVXVComplex1D &,
						      const AVXVComplex1D &)  noexcept(true);

		/*
		@Purpose:
					Mimicks this operation -- C = A*B (vectors size obviously must be equal to each other).
					where A is of type: AVXVComplex1D
					and B is of type: double* __restrict
		@Remark:
					

		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void avxvcomplex_mul(AVXVComplex1D &, const AVXVComplex1D &,
						    const double * __restrict)   noexcept(true);

		/*
		@Purpose:
					Mimicks this operation -- C = A/B (vectors size obviously must be equal to each other).
					where A is of type: AVXVComplex1D
					and   B is of type: AVXVComplex1D
		@Remark:
					In order to enable HW prefetching additional
					argument of type AVXVComplex1D (copy of third argument)
					is passed.

		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void avxcomplex_div(AVXVComplex1D &, const AVXVComplex1D &, 
						    const AVXVComplex1D &, _In_ const AVXVComplex1D &) noexcept(true);

		/*
		@Purpose:
					Mimicks this operation -- C = A/B (vectors size obviously must be equal to each other).
					where A is of type: double * __restrict
					and B is of type:   AVXVComplex1D
		@Remark:
					Size of double * argument must be
					equal to size of AVXVComplex1D.m_Re

		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void avxvcomplex_div(AVXVComplex1D &, const AVXVComplex1D &,
						     const double * __restrict) noexcept(true);

		/*
		@Purpose:
					C = A==B (vectors size obviously must be equal to each other).
					where C is represented as two arrays Re[n] == Re[n] , Im[n] == Im[n]
					A is of type: AVXVComplex1D
					B is of type: AVXVComplex1D
		
		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void avxvcomplex_eq(double * __restrict, double * __restrict,
							const AVXVComplex1D &, _In_ const AVXVComplex1D &) noexcept(true);

		/*
		@Purpose:
					C = A==B (vectors size obviously must be equal to each other).
					where C is represented as two arrays Re[n] != Re[n] , Im[n] != Im[n]
					A is of type: AVXVComplex1D
					B is of type: AVXVComplex1D

		@Warning:
					In order to speed up computation and
					diminish branch predicition misses
					no error checking is implemented (usage of conditional if statements)
		*/
		void avxvcomplex_neq(double * __restrict, double * __restrict,
							 const AVXVComplex1D &, _In_ const AVXVComplex1D &) noexcept(true);
	}
}

#endif /*__GMS_AVXVCPROCEDURES_H__*/
