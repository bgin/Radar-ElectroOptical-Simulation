
#ifndef __GMS_MKL_LAPLACERNG_H__
#define __GMS_MKL_LAPLACERNG_H__

namespace file_info {



const unsigned int gGMS_MKL_LAPLACERNG_MAJOR =  1;

const unsigned int gGMS_MKL_LAPLACERNG_MINOR = 1;

const unsigned int gGMS_MKL_LAPLACERNG_MICRO =  0;

const unsigned int gGMS_MKL_LAPLACERNG_FULLVER = 
	1000U*gGMS_MKL_LAPLACERNG_MAJOR + 100U*gGMS_MKL_LAPLACERNG_MINOR + 10U*gGMS_MKL_LAPLACERNG_MICRO;

const char * const pgGMS_MKL_LAPLACERNG_CREATE_DATE = "25-04-2018 09:25 +00200 (WED 25 APR 2018 GMT+2)";

const char * const pgGMS_MKL_LAPLACERNG_BUILD_DATE = __DATE__ ":" __TIME__;

const char * const pgGMS_MKL_LAPLACERNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_LAPLACERNG_DESCRIPT = "C++ wrapper for Intel MKL vdRngLaplace procedure.";

}




#include <iostream>

#include "GMS_config.h"

#endif
#if (USE_MKL) == 1
#include <mkl_vsl.h>
#else
#error Intel MKL library required for this compilation
#endif

namespace gms {
	namespace math {
		namespace stat {
				
				//
				//	C++ wrapper around Intel MKL vdRngLaplace procedure.
				//

			  __attribute__((align(64))) struct MKLLaplaceRNG {


			                double * __restrict m_rvec;

					double  m_a;

					double  m_beta;
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
					MKL_INT m_seed;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(4,4)
#endif
					MKL_INT m_error;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(5,4)
#endif

					

					

					//
					//	Construction and destruction
					//

					MKLLaplaceRNG();

					MKLLaplaceRNG(const MKL_INT,
						      const MKL_UINT,
						      const MKL_INT,
						      const double,
						      const double);

					MKLLaplaceRNG(const MKLLaplaceRNG &);

					MKLLaplaceRNG(MKLLaplaceRNG &&)noexcept(true);

					~MKLLaplaceRNG()noexcept(true);

					MKLLaplaceRNG & operator=(const MKLLaplaceRNG &);

					MKLLaplaceRNG & operator=(MKLLaplaceRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,const MKLLaplaceRNG &);
		}
	}
}




#endif /*__GMS_MKL_LAPLACERNG_H__*/
