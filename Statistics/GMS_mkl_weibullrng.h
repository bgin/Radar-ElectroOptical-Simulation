
#ifndef __GMS_MKL_WIEBULLRNG_H__
#define __GMS_MKL_WEIBULLRNG_H__ 260420181417

namespace file_info {

  
const unsigned int gGMS_MKL_WEIBULLRNG_MAJOR = 1;

const unsigned int gGMS_MKL_WEIBULLRNG_MINOR = 1;

const unsigned int gGMS_MKL_WEIBULLRNG_MICRO = 0;

const unsigned int gGMS_MKL_WEIBULLRNG_FULLER = 
	1000U*gGMS_MKL_WEIBULLRNG_MAJOR + 100U*gGMS_MKL_WEIBULLRNG_MINOR + 10U*gGMS_MKL_WEIBULLRNG_MICRO;

const char * const pgGMS_MKL_WEIBULLRNG_CREATE_DATE = "26-04-2018 14:17 +00200 (THR 26 APR 2018 GMT+2)";

const char * const pgGMS_MKL_WEIBULLRNG_BUILD_DATE = __DATE__ ":" __TIME__;

const char * const pgGMS_MKL_WEIBULLRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_WEIBULLRNG_DESCRIPT = "C++ wrapper for Intel MKL vdRngWeibull procedure.";

}

#include <iostream>

 #include "GMS_config.h"

#if (USE_MKL) == 1
    #include <mkl_vsl.h>
#else
    #error Intel MKL library is required for this translation unit compilation
#endif

namespace gms {
	namespace math {
		namespace stat {

			//
			//	C++ wrapper for Intel MKL vdRngWeibull procedure
			//

                        __attribute__(align((64))) struct MKLWRNGData {


                                   double * __restrict m_rvec;

				   double m_alpha;
				   double m_a;
				   double m_beta;
				   MKL_INT m_nvalues;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(1, 4)
#endif
				   MKL_INT m_brng;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(2, 4)
#endif
			       MKL_INT m_seed;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(3, 4)
#endif
				   MKL_INT m_error;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(4, 4)
#endif
			};

			

				   

				  
					

                        __attribute__((align(64))) struct MKLWeibullRNG {

				MKLWRNGData data;
					//
					//	Constuction and destruction
					//

					MKLWeibullRNG();

					MKLWeibullRNG(const MKL_INT,
					              const MKL_UINT,
						      const MKL_INT,
						      const double,
						      const double,
					              const double);

					MKLWeibullRNG(const MKLWeibullRNG &);

					MKLWeibullRNG(MKLWeibullRNG &&)noexcept(true);

					~MKLWeibullRNG()noexcept(true);

					MKLWeibullRNG & operator=(const MKLWeibullRNG &);

					MKLWeibullRNG & operator=(MKLWeibullRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,
			           const MKLWeibullRNG &);


		}
	}
}



#endif /*__GMS_WEIBULLRNG_H__*/
