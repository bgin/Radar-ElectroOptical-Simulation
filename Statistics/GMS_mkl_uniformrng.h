
#ifndef __GMS_MKL_UNIFORMRNG_H__
#define __GMS_MKL_UNIFORMRNG_H__

namespace file_info {


  
const unsigned int gGMS_MKL_UNIFORMRNG_MAJOR = 1;

const unsigned int gGMS_MKL_UNIFORMRNG_MINOR = 1;

const unsigned int gGMS_MKL_UNIFORMRNG_MICRO = 0;

const unsigned int gGMS_MKL_UNIFORMRNG_FULLVER = 
	1000U*gGMS_MKL_UNIFORMRNG_MAJOR + 100U*gGMS_MKL_UNIFORMRNG_MINOR + 10U*gGMS_MKL_UNIFORMRNG_MICRO;

const char * const pgGMS_MKL_UNIFORMRNG_CREATE_DATE = "24-04-2018 08:42 +00200 (TUE 24 APR 2018 GMT+2)";

const char * const pgGMS_MKL_UNIFORMRNG_BUILD_DATE = __DATE__ ":"__TIME__;

const char * const pgGMS_MKL_UNIFORMRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_UNIFORMRNG_DESCRIPT = "C++ wrapper over MKL vdRngUniform procedure.";

}

#include <iostream>
#include "GMS_config.h"
#if (USE_MKL) == 1
    #include <mkl_vsl.h>
#else
    #error Intel MKL library is required for compilation of this translation unit
#endif

namespace gms {
	namespace math {
	

			//
			//	C++ wrapper around Intel MKL vdRngUniform procedure.
			//

                        __attribute__((align(64))) struct MKLUniformRNG {


                                    double * __restrict m_rvec;

				    double   m_a; // param 'a'

				    double   m_b;   // param 'b'

#if (USE_STRUCT_PADDING) == 1
					PAD_TO(1,8)
#endif
					
					MKL_INT m_nvalues; // number of random values in vector m_rvec
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(2,4)
#endif
					MKL_UINT m_brng;  // Basic random number generator
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(3,4)
#endif
					MKL_INT  m_seed;  // initial seed 
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(4,4)
#endif
					MKL_INT  m_error; // error returned by the MKL
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(5,4)
#endif
					

					//
					//	Constructor and destructor
					//

					MKLUniformRNG();

					MKLUniformRNG(const MKL_INT,
						      const MKL_UINT,
						      const MKL_INT,
						      const double,
					              const double);

					MKLUniformRNG(const MKLUniformRNG &);

					MKLUniformRNG(MKLUniformRNG &&)noexcept(true);

					~MKLUniformRNG()noexcept(true);

					MKLUniformRNG & operator=(const MKLUniformRNG &);

					MKLUniformRNG & operator=(MKLUniformRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,const MKLUniformRNG &);
		}
	}


#endif /*__GMS_MKL_UNIFORMRNG_H__*/
