
#ifndef __GMS_MKL_CAUCHYRNG_H__
#define __GMS_MKL_CAYCHYRNG_H__ 230420180916

namespace file_info {


  const unsigned int gGMS_MKL_CAUCHYRNG_MAJOR = 1;

  const unsigned int gGMS_MKL_CAUCHYRNG_MINOR = 1;

  const unsigned int gGMS_MKL_CAUCHYRNG_MICRO =  0;

const unsigned int gGMS_MKL_CAUCHYRNG_FULLVER = 
	1000U * gGMS_MKL_CAUCHYRNG_MAJOR + 100U*gGMS_MKL_CAUCHYRNG_MINOR + 10U*gGMS_MKL_CAUCHYRNG_MICRO;

const char * const pgGMS_MKL_CAUCHYRNG_CREATE_DATE = "23-04-2018 09:16 +00200 (MON 23 APR 2018 GMT+2)";

const char * const pgGMS_MKL_CAUCHYRNG_BUILD_DATE = __DATE__ ":" __TIME__;

const char * const pgGMS_MKL_CAUCHYRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_CAUCHYRNG_DESCRIPT = "C++ wrapper over MKL vdRngCauchy procedure.";

}





#include <iostream>

#include "GMS_config.h"

#if (USE_MKL) == 1
#include <mkl_vsl.h>
#endif

namespace gms {
	namespace math {
		namespace stat {

			// Wrapper around Mkl_Api(int,vdRngCauchy,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  ))

			     __attribute__((align(64))) MKLCauchyRNG {


				double * __restrict m_rvec;

				double m_a;       // param 'a'

				double m_beta;   // param 'beta'
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,8)
#endif

				MKL_INT m_nvalues;  // number of random values in vector m_rvec
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2,4)
#endif
				MKL_UINT m_brng;    // Basic random number generator
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3,4)
#endif
				MKL_INT m_seed;     // initial seed
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(4,4)
#endif
				MKL_INT m_error;   // error returned by the MKL
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(5,4)
#endif
				

				

				//
				//	Construction and destruction
				//
				MKLCauchyRNG();

				MKLCauchyRNG(const MKL_INT,
					     const MKL_UINT,
					     const MKL_INT,
					     const double,
					     const double);

				MKLCauchyRNG(const MKLCauchyRNG &);

				MKLCauchyRNG(MKLCauchyRNG &&)noexcept(true);

				~MKLCauchyRNG()noexcept(true);

				MKLCauchyRNG & operator=(const MKLCauchyRNG &);

				MKLCauchyRNG & operator=(MKLCauchyRNG &&)noexcept(true);

				void compute_rand_distribution(const MKL_INT);

				void compute_rand_distribution(VSLStreamStatePtr,
							       const MKL_INT);


			};

			std::ostream &
			operator<<(std::ostream &,const MKLCauchyRNG &);
		}
	}
}



#endif /*__GMS_MKL_CAUCHYRNG_H__*/
