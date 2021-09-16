
#ifndef __GMS_MKL_EXPONENTIALRNG_H__
#define __GMS_MKL_EXPONENTIALRNG_H__ 240420181525

namespace file_info {



const unsigned int gGMS_MKL_EXPONENTIALRNG_MAJOR = 1;

const unsigned int gGMS_MKL_EXPONENTIALRNG_MINOR = 1;

const unsigned int gGMS_MKL_EXPONENTIALRNG_MICRO = 0;

const unsigned int gGMS_MKL_EXPONENTIALRNG_FULLVER = 
	1000U*gGMS_MKL_EXPONENTIALRNG_MAJOR + 100U*gGMS_MKL_EXPONENTIALRNG_MINOR + 10U*gGMS_MKL_EXPONENTIALRNG_MICRO;

const char * const pgGMS_MKL_EXPONENTIALRNG_CREATE_DATE = "24-04-2018 15:25 +00200 (TUE 24 APR 2018 GMT+2)";

const char * const pgGMS_MKL_EXPONENTIALRNG_BUILD_DATE = __DATE__ ":" __TIME__;

const char * const pgGMS_MKL_EXPONENTIALRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_EXPONENTIALRNG_DESCRIPT = "C++ wrapper around Intel MKL vdRngExponential procedure.";

}



#include <iostream>


#include "GMS_config.h"

#if (USE_MKL) == 1
#include <mkl_vsl.h>
#endif

namespace gms {
	namespace math {
		namespace stat {

			// C++ wrapper around Intel MKL vdRngExponential procedure.

			__attribute__((align(64))) struct MKLExponentialRNG {


			                double * __restrict m_rvec;

					double   m_a;

					double   m_beta;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(1,8)
#endif
					
					MKL_INT m_nvalues;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(2,4)
#endif
					MKL_UINT m_brng;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(3,4)
#endif
					MKL_INT  m_seed;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(4,4)
#endif
					MKL_INT  m_error;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(5,4)
#endif
					

				

					//
					//	Constructor and destructor
					//
					MKLExponentialRNG();

					MKLExponentialRNG(const MKL_INT,
							  const MKL_UINT,
							  const MKL_INT,
							  const double,
							  const double);

					MKLExponentialRNG(const MKLExponentialRNG &);

					MKLExponentialRNG(MKLExponentialRNG &&)noexcept(true);

					~MKLExponentialRNG()noexcept(true);

					MKLExponentialRNG & operator=(const MKLExponentialRNG &);

					MKLExponentialRNG & operator=(MKLExponentialRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);

			};

			std::ostream &
			operator<<(std::ostream &,
				   const MKLExponentialRNG &);
		}
	}
}


#endif /*__GMS_MKL_EXPONENTIALRNG_H__*/
