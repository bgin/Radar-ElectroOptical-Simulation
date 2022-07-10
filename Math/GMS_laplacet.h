
#ifndef __GMS_LAPLACET_H__
#define __GMS_LAPLACET_H__



namespace file_info{

const unsigned int GMS_LAPLACET_MAJOR = 1;

const unsigned int GMS_LAPLACET_MINOR = 0;

const unsigned int GMS_LAPLACET_MICRO = 1;

const unsigned int GMS_LAPLACET_FULLVER = 
	1000U*GMS_LAPLACET_MAJOR + 100U*GMS_LAPLACET_MINOR + 10U*GMS_LAPLACET_MICRO;

const char * const GMS_LAPLACET_CREATE_DATE = "10-62-2022 08:45 +00200 (SUN 10 JUL 2022 GMT+2)";


const char * const GMS_LAPLACET_BUILD_DATE =  __DATE__":"__TIME__;

const char * const GMS_LAPLACET_AUTHOR  = "Programmer: Bernard Gingold contact: beniekg@gmail.com";

const char * const GMS_LAPLACET_DESCRIPT = "Laplace Integral Transform computed by the Quadpack Integrator.";

}

#include <complex>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

      namespace math {


          typedef struct __ATTR_ALIGN__(64) QuadErrorParams {
                  // Quadpack Integrators return info.
                  double  * __restrict re_abserr;
		  double  * __restrict im_abserr;
		  int32_t * __restrict re_neval;
		  int32_t * __restrict im_neval;
		  int32_t * __restrict re_ier;
		  int32_t * __restrict im_ier;
		  int32_t * __restrict re_last;
		  int32_t * __restrict im_last;
		  bool                 is_allocated;
#if (USE_STRUCT_PADDING) == 1
                  PAD_TO(0,64)
#endif       
	  } QuadErrorParams;

         using c8 = std::complex<double>;
          
	 bool
	 laplace_params_alloc(QuadErrorParams &,
	                     const int32_t) // Laplace number of points
			      __ATTR_HOT__
			      __ATTR_ALIGN__(32);

	 void
	 laplace_params_free(QuadErrorParams &)
	                      __ATTR_HOT__
			      __ATTR_ALIGN__(32);
	  /*
                The integrands [real and imaginary] parts shall be implemented the following way
                void * user_data points to "s-1" complex power
                double re_func(double t, void * user_data) {
                   complex<double> s = *((complex<double>*)user_data); // s-complex power
                   const double x; // = .... Result of integrand evaluation of argument 't', e.g 't=x^2'
                   double re = s.real();
                   double im = s.imag();
                   double ere = std::exp(-re*t);
                   double cim = ere*std::cos(im*t);
                   return (cim*x);
                }

                 double im_func(double t, void * user_data) {
                   complex<double> s = *((complex<double>*)user_data); // s-complex power
                   const double x; // = .... Result of integrand evaluation of argument 't', e.g 't=x^2'
                   double re = s.real();
                   double im = s.imag();
                   double ere = std::exp(-re*t);
                   double sim = -ere*std::sin(im*t);
                   return (sim*x);
                }
            */
/*
   This parameter was removed: 
              //void * __restrict,                   // DQAGI user_data, i.e. complex vector of 's' values
*/

	  // OpenMP versions

	  bool
	  laplacet_dqagi_omp(double(*)(double,void * __restrict), // Laplace Transform integrand real-part
	                     double(*)(double,void * __restrict), // Laplace Transform integrand imaginary-part
	                     const double,                        // DQAGI bound argument
                             const int32_t,                       // DQAGI inf   argument
			     const double,                        // DQAGI epsabs argument
			     const double,                        // DQAGI epsrel argument
  		             const QuadErrorParams &,              // DQAGI aggregated per real and imaginary integrator error results
			     const int,                           // Laplace number of points
			     c8 *   __restrict)                   // Laplace output complex vector transform data
                                                  __ATTR_HOT__
						  __ATTR_ALIGN__(32);         
			    


	  bool
	  laplacet_dqage_omp(double(*)(double, void * __restrict), // Laplace Transform integrand
	                    double(*)(double, void * __restrict),  // Laplace Transform integrand imaginary-part
	                    const double,                         // lower limit of integration.
			    const double,                         // upper limit of integration.
			    const double,                         // absolute accuracy requested.
			    const double,                         // relative accuracy requested.
			    const int32_t,                        // integration rule to be used
			    const QuadErrorParams &,              // DQAGE aggregated per real and imaginary integrator error results
			     const int,                            // Laplace number of points
			    c8 *   __restrict)                    // Laplace output complex vector transform data
                                                  __ATTR_HOT__
						  __ATTR_ALIGN__(32);                                      
			   

	  bool
	  laplacet_dqagp_omp(double(*)(double, void * __restrict), // Laplace Transform integrand
	                    double(*)(double, void * __restrict),  // Laplace Transform integrand imaginary-part
	                    const double,                         // lower limit of integration.
			    const double,                         // upper limit of integration.
			    const int32_t,                        // number equal to 2 more than the number of sinularities.
			    const double * __restrict,            //  vector of dimension npts2, the first (npts2-2) elements  of which are the user provided interior break points.
                            const double,                         // absolute accuracy requested.
			    const double,                         // relative accuracy requested.
			    const QuadErrorParams &,               // DQAGE aggregated per real and imaginary integrator error results
	                     const int,                            // Laplace number of points
			    c8 *   __restrict)
			                           __ATTR_HOT__
						  __ATTR_ALIGN__(32);                     
			   

          // Single-threaded versions
	  bool
	  laplacet_dqagi(double(*)(double,void * __restrict), // Laplace Transform integrand real-part
	                    double(*)(double,void * __restrict), // Laplace Transform integrand imaginary-part
	                    const double,                        // DQAGI bound argument
                            const int32_t,                       // DQAGI inf   argument
			    const double,                        // DQAGI epsabs argument
			    const double,                        // DQAGI epsrel argument
			    const QuadErrorParams &,              // DQAGI aggregated per real and imaginary integrator error results
			    const int,                           // Laplace number of points
			    c8 *   __restrict)                   // Laplace output complex vector transform data
                                                  __ATTR_HOT__
						  __ATTR_ALIGN__(32);                                    
			   


	  bool
	  laplacet_dqage(double(*)(double, void * __restrict), // Laplace Transform integrand
	                    double(*)(double, void * __restrict),  // Laplace Transform integrand imaginary-part
	                    const double,                         // lower limit of integration.
			    const double,                         // upper limit of integration.
			    const double,                         // absolute accuracy requested.
			    const double,                         // relative accuracy requested.
			    const int32_t,                        // integration rule to be used
			    const QuadErrorParams &,              // DQAGE aggregated per real and imaginary integrator error results
			     const int,                            // Laplace number of points
			    c8 *   __restrict)                    // Laplace output complex vector transform data
                                                  __ATTR_HOT__
						  __ATTR_ALIGN__(32);                                       
			   

	  bool
	  laplacet_dqagp(double(*)(double, void * __restrict), // Laplace Transform integrand
	                    double(*)(double, void * __restrict),  // Laplace Transform integrand imaginary-part
	                    const double,                         // lower limit of integration.
			    const double,                         // upper limit of integration.
			    const int32_t,                        // number equal to 2 more than the number of sinularities.
			    const double * __restrict,            //  vector of dimension npts2, the first (npts2-2) elements  of which are the user provided interior break points.
                            const double,                         // absolute accuracy requested.
			    const double,                         // relative accuracy requested.
			    const QuadErrorParams &,              // DQAGE aggregated per real and imaginary integrator error results
			    const int,                            // Laplace number of points
			    c8 *   __restrict)
			                           __ATTR_HOT__
						  __ATTR_ALIGN__(32);                                      


						  
   } // math


} // gms


#endif /*__GMS_LAPLACET_H__*/








