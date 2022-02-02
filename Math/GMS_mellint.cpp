

#include <omp.h>
#include <cmath>
#include "GMS_mellint.h"
#include "GMS_malloc.h"
#include "GMS_cquadpack.h"




bool
gms::math::
mellint_dqagi_omp(double(*re_f)(double t, void * __restrict user_data),
                  double(*im_f)(double t, void * __restrict user_data),
		  const double bound,
		  const int32_t inf,
		  const double epsabs,
		  const double epsrel,
		  const QuadErrorParams &params,
		  void * __restrict input, // values of complex power argument ,i.e. s
		  const int32_t npts,
		  c8 * __restrict output) {
		  
      if(__builtin_expect(!mp.is_allocated,0)) { return (false);}
      const c8 I(0.0,1.0);
      double * __restrict reabser = params.re_abserr;
      double * __restrict imabser = params.im_abserr;
      int32_t* __restrict reneval = params.re_neval;
      int32_t* __restrict imneval = params.im_neval;
      int32_t* __restrict reier   = params.re_ier;
      int32_t* __restrict imier   = params.im_ier;
      double re = 0.0;
      double im = 0.0;

#pragma omp parallel for schedule(static) default(none)         \
        shared(re_f,im_f,reabserr,imabser,reneval,imneval,reier,imier, \
	input,output,npts,bound,inf,epsabs,epsrel,I)           \
	firstprivate(re,im) private(i)
      for(int32_t i = 0; i < npts; ++i) {
          re = dqagi(re_f,bound,inf,epsabs,epsrel,&reabser[i],
	             &reneval[i],&reier[i],&input[i]);
	  im = dqagi(im_f,bound,inf,epsabs,epsrel,&imabser[i],
	             &imneval[i],&imier[i],&input[i]);
	  output[i] = re+I*im;
      }
      return (true);
}


bool
gms::math::
mellint_dqage_omp(double(*re_f)(double t, void * __restrict user_data), // Mellin Transform integrand real-part
                  double(*im_f)(double t, void * __restrict user_data),
	          const double a,                         // lower limit of integration.
		  const double b,                         // upper limit of integration.
		  const double epsabs,                    // absolute accuracy requested.
		  const double epsrel,                    // relative accuracy requested.
		  const int32_t irule,                    // integration rule to be used
		  const QuadErrorParams & params          // DQAGE aggregated per real and imaginary integrator error results
		  void * __restrict input,                // DQAGE user_data, i.e. complex vector of 's' values
		  const int npts,                         // Mellin number of points
		  c8 *   __restrict output) {             // Result of Mellin Transform
      if(__builtin_expect(!mp.is_allocated,0)) { return (false);}
      const c8 I(0.0,1.0);
      double * __restrict reabser = params.re_abserr;
      double * __restrict imabser = params.im_abserr;
      int32_t* __restrict reneval = params.re_neval;
      int32_t* __restrict imneval = params.im_neval;
      int32_t* __restrict reier   = params.re_ier;
      int32_t* __restrict imier   = params.im_ier;
      int32_t* __restrict relast  = params.re_last;
      int32_t* __restrict imlast  = params.im_last;
      double re = 0.0;
      double im = 0.0;
#pragma omp parallel for schedule(static) default(none)            \
        shared(re_f,im_f,reabserr,imabser,reneval,imneval,reier,imier,    \
	       a,b,relast,imlast,input,output,npts,irule,epsabs,epsrel,I) \
	firstprivate(re,im) private(i)
	for(int32_t i = 0; i < npts; ++i) {
            re = dqage(re_f,a,b,epsabs,epsrel,irule,&reabser[i],
	               &reneval[i],&reier[i],&relast[i],&input[i]);
	    im = dqage(im_f,a,b,epsabs,epsrel,irule,&imabser[i],
	               &imneval[i],&imier[i],imlast[i],&input[i]);
	    output[i] = re+I*im;
	}
	return (true);
}


bool
gms::math::
mellint_dqagp_omp(double(*re_f)(double t, void * __restrict user_data), // Mellin Transform integrand
	          double(*im_f)(double t, void * __restrict user_data),  // Mellin Transform integrand imaginary-part
	          const double a,                         // lower limit of integration.
		  const double b,                         // upper limit of integration.
		  const int32_t nsng,                    // number equal to 2 more than the number of sinularities.
		  const double * __restrict sng,         //  vector of dimension nsng, the first (nsng2-2) elements  of which are the user provided interior break points.
                  const double epsabs,                         // absolute accuracy requested.
		  const double epsrel,                         // relative accuracy requested.
		  const QuadErrorParams &params               // DQAGE aggregated per real and imaginary integrator error results
		  void * __restrict input,                    // DQAGE user_data, i.e. complex vector of 's' values
		  const int npts,                            // Mellin number of points
		  c8 *   __restrict output) {

      if(__builtin_expect(!mp.is_allocated,0)) { return (false);}
      const c8 I(0.0,1.0);
      double * __restrict reabser = params.re_abserr;
      double * __restrict imabser = params.im_abserr;
      int32_t* __restrict reneval = params.re_neval;
      int32_t* __restrict imneval = params.im_neval;
      int32_t* __restrict reier   = params.re_ier;
      int32_t* __restrict imier   = params.im_ier;
      double re = 0.0;
      double im = 0.0;
#pragma omp parallel for schedule(static) default(none)            \
        shared(re_f,im_f,reabserr,imabser,reneval,imneval,reier,imier,    \
	       a,b,nsng,sng,input,output,npts,epsabs,epsrel,I)		\
	firstprivate(re,im) private(i)
      	for(int32_t i = 0; i < npts; ++i) {
            re = dqagp(re_f,a,b,nsng,sng,epsabs,epsrel,&reabser[i],
	               &reneval[i],&reier[i],&input[i]);
	    im = dqagp(re_f,im_f,a,b,nsng,sng,epsabs,epsrel,&imabser[i],
	               &imneval[i],&imier[i],&input[i]);
	    output[i] = re+I*im;
	}
	return (true);
}

//    Single-threaded versions

bool
gms::math::
mellint_dqagi(double(*re_f)(double t, void * __restrict user_data),
                  double(*im_f)(double t, void * __restrict user_data),
		  const double bound,
		  const int32_t inf,
		  const double epsabs,
		  const double epsrel,
		  const QuadErrorParams &params,
		  void * __restrict input, // values of complex power argument ,i.e. s
		  const int32_t npts,
		  c8 * __restrict output) {
		  
      if(__builtin_expect(!mp.is_allocated,0)) { return (false);}
      const c8 I(0.0,1.0);
      double * __restrict reabser = params.re_abserr;
      double * __restrict imabser = params.im_abserr;
      int32_t* __restrict reneval = params.re_neval;
      int32_t* __restrict imneval = params.im_neval;
      int32_t* __restrict reier   = params.re_ier;
      int32_t* __restrict imier   = params.im_ier;
      double re = 0.0;
      double im = 0.0;

      for(int32_t i = 0; i < npts; ++i) {
          re = dqagi(re_f,bound,inf,epsabs,epsrel,&reabser[i],
	             &reneval[i],&reier[i],&input[i]);
	  im = dqagi(im_f,bound,inf,epsabs,epsrel,&imabser[i],
	             &imneval[i],&imier[i],&input[i]);
	  output[i] = re+I*im;
      }
      return (true);
}


bool
gms::math::
mellint_dqage(double(*re_f)(double t, void * __restrict user_data), // Mellin Transform integrand real-part
                  double(*im_f)(double t, void * __restrict user_data),
	          const double a,                         // lower limit of integration.
		  const double b,                         // upper limit of integration.
		  const double epsabs,                    // absolute accuracy requested.
		  const double epsrel,                    // relative accuracy requested.
		  const int32_t irule,                    // integration rule to be used
		  const QuadErrorParams & params          // DQAGE aggregated per real and imaginary integrator error results
		  void * __restrict input,                // DQAGE user_data, i.e. complex vector of 's' values
		  const int npts,                         // Mellin number of points
		  c8 *   __restrict output) {             // Result of Mellin Transform
      if(__builtin_expect(!mp.is_allocated,0)) { return (false);}
      const c8 I(0.0,1.0);
      double * __restrict reabser = params.re_abserr;
      double * __restrict imabser = params.im_abserr;
      int32_t* __restrict reneval = params.re_neval;
      int32_t* __restrict imneval = params.im_neval;
      int32_t* __restrict reier   = params.re_ier;
      int32_t* __restrict imier   = params.im_ier;
      int32_t* __restrict relast  = params.re_last;
      int32_t* __restrict imlast  = params.im_last;
      double re = 0.0;
      double im = 0.0;
      for(int32_t i = 0; i < npts; ++i) {
            re = dqage(re_f,a,b,epsabs,epsrel,irule,&reabser[i],
	               &reneval[i],&reier[i],&relast[i],&input[i]);
	    im = dqage(im_f,a,b,epsabs,epsrel,irule,&imabser[i],
	               &imneval[i],&imier[i],imlast[i],&input[i]);
	    output[i] = re+I*im;
	}
	return (true);
}


bool
gms::math::
mellint_dqagp(double(*re_f)(double t, void * __restrict user_data), // Mellin Transform integrand
	          double(*im_f)(double t, void * __restrict user_data),  // Mellin Transform integrand imaginary-part
	          const double a,                         // lower limit of integration.
		  const double b,                         // upper limit of integration.
		  const int32_t nsng,                    // number equal to 2 more than the number of sinularities.
		  const double * __restrict sng,         //  vector of dimension nsng, the first (nsng2-2) elements  of which are the user provided interior break points.
                  const double epsabs,                         // absolute accuracy requested.
		  const double epsrel,                         // relative accuracy requested.
		  const QuadErrorParams &params               // DQAGE aggregated per real and imaginary integrator error results
		  void * __restrict input,                    // DQAGE user_data, i.e. complex vector of 's' values
		  const int npts,                            // Mellin number of points
		  c8 *   __restrict output) {

      if(__builtin_expect(!mp.is_allocated,0)) { return (false);}
      const c8 I(0.0,1.0);
      double * __restrict reabser = params.re_abserr;
      double * __restrict imabser = params.im_abserr;
      int32_t* __restrict reneval = params.re_neval;
      int32_t* __restrict imneval = params.im_neval;
      int32_t* __restrict reier   = params.re_ier;
      int32_t* __restrict imier   = params.im_ier;
      double re = 0.0;
      double im = 0.0;
      for(int32_t i = 0; i < npts; ++i) {
            re = dqagp(re_f,a,b,nsng,sng,epsabs,epsrel,&reabser[i],
	               &reneval[i],&reier[i],&input[i]);
	    im = dqagp(re_f,im_f,a,b,nsng,sng,epsabs,epsrel,&imabser[i],
	               &imneval[i],&imier[i],&input[i]);
	    output[i] = re+I*im;
	}
	return (true);
}





bool
gms::math::
mellin_params_alloc(const QuadErrorParams &mp,
		    const int32_t npts) {
      using namespace gms::common;		    
      if(__builtin_expect(npts<1,0)) {return (false);}
      const std::size_t len = (std::size_t)npts;
      constexpr std::size_t align = (std::size_t)64;
      mp.re_abserr = (double*) gms_mm_malloc(len,align);
      mp.im_abserr = (double*) gms_mm_malloc(len,align);
      mp.re_neval  = (int32_t*)gms_mm_malloc(len,align);
      mp.im_neval  = (int32_t*)gms_mm_malloc(len,align);
      mp.re_ier    = (int32_t*)gms_mm_malloc(len,align);
      mp.im_ier    = (int32_t*)gms_mm_malloc(len,align);
      mp.re_last   = (int32_t*)gms_mm_malloc(len,align);
      mp.im_last   = (int32_t*)gms_mm_malloc(len,align);
      mp.is_allocated = true;
      return (true);
}

void
gms::math::
mellin_params_free(const QuadErrorParams  &mp) {
      using namespace gms::common;
      if(__builtin_expect(!mp.is_allocated,0)) {return;}
      gms_mm_free(mp.im_last);
      gms_mm_free(mp.re_last);
      gms_mm_free(mp.im_ier);
      gms_mm_free(mp.re_ier);
      gms_mm_free(mp.im_neval);
      gms_mm_free(mp.re_neval);
      gms_mm_free(mp.im_abserr);
      gms_mm_free(mp.re_abserr);
      mp.is_allocated = false;
}
