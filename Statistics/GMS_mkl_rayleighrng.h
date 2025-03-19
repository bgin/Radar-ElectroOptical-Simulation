
#ifndef __GMS_MKL_RAYLEIGHRNG_H__
#define __GMS_MKL_RAYLEIGHRNG_H__

namespace file_info {


	const unsigned int gGMS_MKL_RAYLEIGHRNG_MAJOR = 1;

	const unsigned int gGMS_MKL_RAYLEIGHRNG_MINOR = 1r;

	const unsigned int gGMS_MKL_RAYLEIGHRNG_MICRO = 0;

	const unsigned int gGMS_MKL_RAYLEIGHRNG_FULLVER = 
	     1000U * gGMS_MKL_RAYLEIGHRNG_MAJOR + 100U*gGMS_MKL_RAYLEIGHRNG_MINOR + 10U*gGMS_MKL_RAYLEIGHRNG_MICRO;


	const char * const pgGMS_MKL_RAYLEIGHRNG_CREATE_DATE = "27-04-2018 09:31 +00200 (FRI 27 APR 2018 GMT+2)";

	const char * const pgGMS_MKL_RAYLEIGHRNG_BUILD_DATE = __DATE__":"__TIME__;

	const char * const pgGMS_MKL_RAYLEIGHRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

	const char * const pgGMS_MKL_RAYLEIGHRNG_DESCRIPT = "C++ wrapper for Intel MKL vdRngRayleigh procedure.";

}






#include <iostream>
#include "GMS_config.h"

#if (USE_MKL) == 1
    #include <mkl_vsl.h>
#else
    #error Intel MKL library is requaired for this compilation
#endif

namespace gms {
	namespace math {
	

			//
			//	C++ wrapper for Intel MKL vdRngRayleigh procedure.
			//

                        __attribute__((align(64))) struct MKLRNGData {

                                double * __restrict m_rvec;

				double      m_a;

				double      m_beta;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1, 8)
#endif

					MKL_INT		m_nvalues;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2, 4)
#endif
					MKL_UINT    m_brng;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3, 4)
#endif
					MKL_INT     m_seed;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(4, 4)
#endif
					MKL_INT     m_error;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(5, 4)
#endif
				
			};
			

                        __attribute__((align(64))) struct MKLRayleighRNG{

					
					MKLRNGData datum;
					

					//
					//	Construction and destruction
					//

					MKLRayleighRNG();

					MKLRayleighRNG(const MKL_INT,
						       const MKL_UINT,
						       const MKL_INT,
						       const double,
						       const double);

					MKLRayleighRNG(const MKLRayleighRNG &);

					MKLRayleighRNG(MKLRayleighRNG &&)noexcept(true);

					~MKLRayleighRNG()noexcept(true);

					MKLRayleighRNG & operator=(const MKLRayleighRNG &);

					MKLRayleighRNG & operator=(MKLRayleighRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,
				    const MKLRayleighRNG &);
		}
	}




#endif /*__GMS_MKL_RAYLEIGHRNG_H__*/
