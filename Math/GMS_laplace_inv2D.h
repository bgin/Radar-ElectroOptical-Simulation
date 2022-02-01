
#ifndef __GMS_LAPLACE_INV2D_H__
#define __GMS_LAPLACE_INV2D_H__



namespace file_info{

const unsigned int GMS_LAPLACE_INV2D_MAJOR = 1;

const unsigned int GMS_LAPLACE_INV2D_MINOR = 0;

const unsigned int GMS_LAPLACE_INV2D_MICRO = 1;

const unsigned int GMS_LAPLACE_INV2D_FULLVER = 
	1000U*GMS_LAPLACE_INV2D_MAJOR + 100U*GMS_LAPLACE_INV2D_MINOR + 10U*GMS_LAPLACE_INV2D_MICRO;

const char * const GMS_LAPLACE_INV2D_CREATE_DATE = "08-10-2019 19:05 +00200 (TUE 08 OCT 2019 GMT+2)";


const char * const GMS_LAPLACE_INV2D_BUILD_DATE =  __DATE__":"__TIME__;

const char * const GMS_LAPLACE_INV2D_AUTHOR  = "Programmer: Bernard Gingold contact: beniekg@gmail.com";

const char * const GMS_LAPLACE_INV2D_DESCRIPT = "Various implementation of Numerical Inverse Laplace Transform 2D.";

}

#include <cstdint>
#include "GMS_config.h"


namespace gms {

        namespace  math {


	           // Helpers, signum function
		   template<typename T>
		             int32_t sgn(T val) {
                        return (T(0) < val) - (val < T(0));
                   }

	           /*
                           Reference:

                           H. Stehfest, Algorithm 368, Numerical Inversion of Laplace 
                           Transforms, CACM Vol. 13, No.1, 47-49 (1970)
	                   J. Abate, W. Whitt, A Unified Approach for Numerically Inverting 
	                   Laplace Transforms, INFORNS Journal of Computing Vol 18, No. 4,
	                   408-421 (2006) 
                      */

	           bool gaver_stehfest_r8_1(double(*)(double,double),
		                            double * __restrict, // output, matrix of calculated value
				            double * __restrict, // input,  vector Gaver-Stehfest coefficients (abscissas)
				            double * __restrict, // input,  vector Gaver Standing Coefficients (Multipliers)
				            double * __restrict, // input,  vector Vector of t1,t2 values ​​for the f calculated
				            const double,              // input,  coeff c1
				            const double,              // input,  coeff c2
				            const int32_t,             // input,  number of approximations
				            const int32_t,             // input,  half number of coeffs
				            int32_t &)                  // output, number of algorithm iterations
				                  __ATTR_HOT__
						  __ATTR_ALIGN__(32)
						  __ATTR_ALL_TARGETS_CLONES__;

                   bool gaver_stehfest_r4_1(float(*)(float,float),
		                            float * __restrict, // output, matrix of calculated value
				            float * __restrict, // input,  vector Gaver-Stehfest coefficients (abscissas)
				            float * __restrict, // input,  vector Gaver Standing Coefficients (Multipliers)
				            float * __restrict, // input,  vector Vector of t1,t2 values ​​for the f calculated
				            const float,               // input,  coeff c1
				            const float,               // input,  coeff c2
				            const int32_t,             // input,  number of approximations
				            const int32_t,             // input,  half number of coeffs
				            int32_t &)                  // output, number of algorithm iterations
				                  __ATTR_HOT__
						  __ATTR_ALIGN__(32)
						  __ATTR_ALL_TARGETS_CLONES__;


		   bool gaver_stehfest_coff_r8_1(int32_t, // input,  half number of coeffs
						 double * __restrict, // output, vector Gaver-Stehfest coefficients (abscissas)
						 double * __restrict) // output, vector Gaver Standing Coefficients (Multipliers)
                                                           __ATTR_HOT__
						           __ATTR_ALIGN__(32)
						           __ATTR_ALL_TARGETS_CLONES__;
		   
                   bool gaver_stehfest_coff_r4_1(int32_t, // input,  half number of coeffs
						 float * __restrict, // output, vector Gaver-Stehfest coefficients (abscissas)
						 float * __restrict) // output, vector Gaver Standing Coefficients (Multipliers)
                                                           __ATTR_HOT__
						           __ATTR_ALIGN__(32)
						           __ATTR_ALL_TARGETS_CLONES__;
  }// math

}// gms























#endif /*__GMS_LAPLACE_INV2D*/
