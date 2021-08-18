
#ifndef __GMS_PMC_EVENTS_TIME_SERIES_ANALYSIS_HPP__
#define __GMS_PMC_EVENTS_TIME_SERIES_ANALYSIS_HPP__

#include <cstdint>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "GMS_config.h"
#include "Timsac_iface.h"






__attribute__((aligned(64))) typedef struct {

      double * __restrict pmc0;
      double * __restrict pmc1;
      double * __restrict pmc2;
      double * __restrict pmc3;
      double * __restrict pmc4;
      double * __restrict pmc5;
      double * __restrict pmc6;
      
} PMC06_samples;



/*
    Apply Time-Series analysis (Timsac) subroutine "CANARM".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "CANARM".
    The number of analyzed events is fixed i.e. 7 events.
    Three events are architectural and four events are
    programmable.
    The passed arrays must have the same size!!
   
*/
__attribute__((hot))
__attribute__((aligned(32)))
template<int32_t len, int32_t lagh>
void
pmc06_events_canarm_analysis_omp(const PMC06_samples * __restrict samples,
                               const char * __restrict fname,
			       const std::string * __restrict event_names) {
    static_assert(len <= 100000, "Input data length can not exceed -- **100000** elements!!");
    const std::size_t  lag2len = (std::size_t)(lagh*lagh);
    const std::size_t  lag3len = (std::size_t)lag2len*len;
    const PMC06_samples * __restrict in_data = samples;
    FILE * fp = NULL;
    if(__builtin_except(fopen(&fp,fname,"a+"),0) != 0) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
      }
    __attribute__((aligned(64))) double PMC0_acor[lagh];
    __attribute__((aligned(64))) double PMC0_acov[lagh];
    __attribute__((aligned(64))) double PMC0_xarcoef[lagh];
    __attribute__((aligned(64))) double PMC0_xv[lagh];
    __attribute__((aligned(64))) double PMC0_xaic[lagh];
    __attribute__((aligned(64))) double PMC0_xparcor[lagh];
    __attribute__((aligned(64))) double PMC0_xdicm[lagh];
    __attribute__((aligned(64))) double PMC0_xb[lagh];
    __attribute__((aligned(64))) double PMC0_xa[lagh];
    __attribute__((aligned(64))) int32_t PMC0_xm1[lagh];
    __attribute__((aligned(64))) int32_t PMC0_xm2[lagh];
    __attribute__((aligned(64))) int32_t PMC0_xpo[lagh];
    __attribute__((aligned(64))) double PMC1_acor[lagh];
    __attribute__((aligned(64))) double PMC1_acov[lagh];
    __attribute__((aligned(64))) double PMC1_xarcoef[lagh];
    __attribute__((aligned(64))) double PMC1_xv[lagh];
    __attribute__((aligned(64))) double PMC1_xaic[lagh];
    __attribute__((aligned(64))) double PMC1_xparcor[lagh];
    __attribute__((aligned(64))) double PMC1_xdicm[lagh];
    __attribute__((aligned(64))) double PMC1_xb[lagh];
    __attribute__((aligned(64))) double PMC1_xa[lagh];
    __attribute__((aligned(64))) int32_t PMC1_xm1[lagh];
    __attribute__((aligned(64))) int32_t PMC1_xm2[lagh];
    __attribute__((aligned(64))) int32_t PMC1_xpo[lagh];
    __attribute__((aligned(64))) double PMC2_acor[lagh];
    __attribute__((aligned(64))) double PMC2_acov[lagh];
    __attribute__((aligned(64))) double PMC2_xarcoef[lagh];
    __attribute__((aligned(64))) double PMC2_xv[lagh];
    __attribute__((aligned(64))) double PMC2_xaic[lagh];
    __attribute__((aligned(64))) double PMC2_xparcor[lagh];
    __attribute__((aligned(64))) double PMC2_xdicm[lagh];
    __attribute__((aligned(64))) double PMC2_xb[lagh];
    __attribute__((aligned(64))) double PMC2_xa[lagh];
    __attribute__((aligned(64))) int32_t PMC2_xm1[lagh];
    __attribute__((aligned(64))) int32_t PMC2_xm2[lagh];
    __attribute__((aligned(64))) int32_t PMC2_xpo[lagh];
    __attribute__((aligned(64))) double PMC3_acor[lagh];
    __attribute__((aligned(64))) double PMC3_acov[lagh];
    __attribute__((aligned(64))) double PMC3_xarcoef[lagh];
    __attribute__((aligned(64))) double PMC3_xv[lagh];
    __attribute__((aligned(64))) double PMC3_xaic[lagh];
    __attribute__((aligned(64))) double PMC3_xparcor[lagh];
    __attribute__((aligned(64))) double PMC3_xdicm[lagh];
    __attribute__((aligned(64))) double PMC3_xb[lagh];
    __attribute__((aligned(64))) double PMC3_xa[lagh];
    __attribute__((aligned(64))) int32_t PMC3_xm1[lagh];
    __attribute__((aligned(64))) int32_t PMC3_xm2[lagh];
    __attribute__((aligned(64))) int32_t PMC3_xpo[lagh];
    __attribute__((aligned(64))) double PMC4_acor[lagh];
    __attribute__((aligned(64))) double PMC4_acov[lagh];
    __attribute__((aligned(64))) double PMC4_xarcoef[lagh];
    __attribute__((aligned(64))) double PMC4_xv[lagh];
    __attribute__((aligned(64))) double PMC4_xaic[lagh];
    __attribute__((aligned(64))) double PMC4_xparcor[lagh];
    __attribute__((aligned(64))) double PMC4_xdicm[lagh];
    __attribute__((aligned(64))) double PMC4_xb[lagh];
    __attribute__((aligned(64))) double PMC4_xa[lagh];
    __attribute__((aligned(64))) int32_t PMC4_xm1[lagh];
    __attribute__((aligned(64))) int32_t PMC4_xm2[lagh];
    __attribute__((aligned(64))) int32_t PMC4_xpo[lagh];
    __attribute__((aligned(64))) double PMC5_acor[lagh];
    __attribute__((aligned(64))) double PMC5_acov[lagh];
    __attribute__((aligned(64))) double PMC5_xarcoef[lagh];
    __attribute__((aligned(64))) double PMC5_xv[lagh];
    __attribute__((aligned(64))) double PMC5_xaic[lagh];
    __attribute__((aligned(64))) double PMC5_xparcor[lagh];
    __attribute__((aligned(64))) double PMC5_xdicm[lagh];
    __attribute__((aligned(64))) double PMC5_xb[lagh];
    __attribute__((aligned(64))) double PMC5_xa[lagh];
    __attribute__((aligned(64))) int32_t PMC5_xm1[lagh];
    __attribute__((aligned(64))) int32_t PMC5_xm2[lagh];
    __attribute__((aligned(64))) int32_t PMC5_xpo[lagh];
    __attribute__((aligned(64))) double PMC6_acor[lagh];
    __attribute__((aligned(64))) double PMC6_acov[lagh];
    __attribute__((aligned(64))) double PMC6_xarcoef[lagh];
    __attribute__((aligned(64))) double PMC6_xv[lagh];
    __attribute__((aligned(64))) double PMC6_xaic[lagh];
    __attribute__((aligned(64))) double PMC6_xparcor[lagh];
    __attribute__((aligned(64))) double PMC6_xdicm[lagh];
    __attribute__((aligned(64))) double PMC6_xb[lagh];
    __attribute__((aligned(64))) double PMC6_xa[lagh];
    __attribute__((aligned(64))) int32_t PMC6_xm1[lagh];
    __attribute__((aligned(64))) int32_t PMC6_xm2[lagh];
    __attribute__((aligned(64))) int32_t PMC6_xpo[lagh];
   
    double  __restrict * PMC0_xw   = NULL;
    double  __restrict * PMC0_xz   = NULL;
    double  __restrict * PMC0_xRs  = NULL;
    double  __restrict * PMC0_xchi = NULL;
    int32_t __restrict * PMC0_xndt = NULL;
    double __restrict  * PMC0_xdic = NULL;
    double  __restrict * PMC1_xw   = NULL;
    double  __restrict * PMC1_xz   = NULL;
    double  __restrict * PMC1_xRs  = NULL;
    double  __restrict * PMC1_xchi = NULL;
    int32_t __restrict * PMC1_xndt = NULL;
    double __restrict  * PMC1_xdic = NULL;
    double  __restrict * PMC2_xw   = NULL;
    double  __restrict * PMC2_xz   = NULL;
    double  __restrict * PMC2_xRs  = NULL;
    double  __restrict * PMC2_xchi = NULL;
    int32_t __restrict * PMC2_xndt = NULL;
    double __restrict  * PMC2_xdic = NULL;
    double  __restrict * PMC3_xw   = NULL;
    double  __restrict * PMC3_xz   = NULL;
    double  __restrict * PMC3_xRs  = NULL;
    double  __restrict * PMC3_xchi = NULL;
    int32_t __restrict * PMC3_xndt = NULL;
    double __restrict  * PMC3_xdic = NULL;
    double  __restrict * PMC4_xw   = NULL;
    double  __restrict * PMC4_xz   = NULL;
    double  __restrict * PMC4_xRs  = NULL;
    double  __restrict * PMC4_xchi = NULL;
    int32_t __restrict * PMC4_xndt = NULL;
    double __restrict  * PMC4_xdic = NULL;
    double  __restrict * PMC5_xw   = NULL;
    double  __restrict * PMC5_xz   = NULL;
    double  __restrict * PMC5_xRs  = NULL;
    double  __restrict * PMC5_xchi = NULL;
    int32_t __restrict * PMC5_xndt = NULL;
    double __restrict  * PMC5_xdic = NULL;
    double  __restrict * PMC6_xw   = NULL;
    double  __restrict * PMC6_xz   = NULL;
    double  __restrict * PMC6_xRs  = NULL;
    double  __restrict * PMC6_xchi = NULL;
    int32_t __restrict * PMC6_xndt = NULL;
    double __restrict  * PMC6_xdic = NULL;
    double  PMC0_xoaic = 0.0;
    double  PMC0_xmean = 0.0;
    int32_t PMC0_xmo  = 0;
    int32_t PMC0_xnc  = 0;
    int32_t PMC0_xk   = 0;
    int32_t PMC0_xl   = 0;
    double PMC1_xoaic = 0.0;
    double PMC1_xmean = 0.0;
    int32_t PMC1_xmo  = 0;
    int32_t PMC1_xnc  = 0;
    int32_t PMC1_xk   = 0;
    int32_t PMC1_xl   = 0;
    double PMC2_xoaic = 0.0;
    double PMC2_xmean = 0.0;
    int32_t PMC2_xmo  = 0;
    int32_t PMC2_xnc  = 0;
    int32_t PMC2_xk   = 0;
    int32_t PMC2_xl   = 0;
    double PMC3_xoaic = 0.0;
    double PMC3_xmean = 0.0;
    int32_t PMC3_xmo  = 0;
    int32_t PMC3_xnc  = 0;
    int32_t PMC3_xk   = 0;
    int32_t PMC3_xl   = 0;
    double PMC4_xoaic = 0.0;
    double PMC4_xmean = 0.0;
    int32_t PMC4_xmo  = 0;
    int32_t PMC4_xnc  = 0;
    int32_t PMC4_xk   = 0;
    int32_t PMC4_xl   = 0;
    double PMC5_xoaic = 0.0;
    double PMC5_xmean = 0.0;
    int32_t PMC5_xmo  = 0;
    int32_t PMC5_xnc  = 0;
    int32_t PMC5_xk   = 0;
    int32_t PMC5_xl   = 0;
    double PMC6_xoaic = 0.0;
    double PMC6_xmean = 0.0;
    int32_t PMC6_xmo  = 0;
    int32_t PMC6_xnc  = 0;
    int32_t PMC6_xk   = 0;
    int32_t PMC6_xl   = 0;
    // OpenMP multithreaded calls to _mm_malloc (using parallel sections)
    // Multithreaded allocation for large dynamic arrays.
#pragma omp parallel sections
    {

         #pragma section
	   {
               PMC0_xw   = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	       PMC0_xz   = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC0_xRs  = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC0_xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC0_xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	       PMC0_xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   }
	 #pragma section
	   {
               PMC1_xw   = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	       PMC1_xz   = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC1_xRs  = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC1_xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC1_xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	       PMC1_xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   }
	 #pragma section
	   {
               PMC2_xw   = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	       PMC2_xz   = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC2_xRs  = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC2_xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC2_xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	       PMC2_xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   }
	 #pragma section
	   {
               PMC3_xw   = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	       PMC3_xz   = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC3_xRs  = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC3_xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC3_xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	       PMC3_xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   }
	 #pragma section
	   {
               PMC4_xw   = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	       PMC4_xz   = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC4_xRs  = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC4_xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC4_xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	       PMC4_xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   }
	 #pragma section
	   {
               PMC5_xw   = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	       PMC5_xz   = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC5_xRs  = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC5_xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC5_xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	       PMC5_xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   }
	 #pragma section
	   {
               PMC6_xw   = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	       PMC6_xz   = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC6_xRs  = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC6_xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	       PMC6_xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	       PMC6_xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   }
        
    }  // omp parallel sections

    const bool isPMC0_NULL = (NULL==PMC0_xw)   || (NULL==PMC0_xz)   || (NULL==PMC0_xRs) ||
                             (NULL==PMC0_xchi) || (NULL==PMC0_xndt) || (NULL==PMC0_xdic);
    if(__builtin_expect(isPMC0_NULL,0)) { MALLOC_FAILED}
    const bool isPMC1_NULL = (NULL==PMC1_xw)   || (NULL==PMC1_xz)   || (NULL==PMC1_xRs) ||
                             (NULL==PMC1_xchi) || (NULL==PMC1_xndt) || (NULL==PMC1_xdic);
    if(__builtin_except(isPMC1_NULL,0)) { MALLOC_FAILED}
    const bool isPMC2_NULL = (NULL==PMC2_xw)   || (NULL==PMC2_xz)   || (NULL==PMC2_xRs) ||
                             (NULL==PMC2_xchi) || (NULL==PMC2_xndt) || (NULL==PMC2_xdic);
    if(__builtin_except(isPMC2_NULL,0)) { MALLOC_FAILED}
    const bool isPMC3_NULL = (NULL==PMC3_xw)   || (NULL==PMC3_xz)   || (NULL==PMC3_xRs) ||
                             (NULL==PMC3_xchi) || (NULL==PMC3_xndt) || (NULL==PMC3_xdic);
    if(__builtin_except(isPMC3_NULL,0)) { MALLOC_FAILED}
    const bool isPMC4_NULL = (NULL==PMC4_xw)   || (NULL==PMC4_xz)   || (NULL==PMC4_xRs) ||
                             (NULL==PMC4_xchi) || (NULL==PMC4_xndt) || (NULL==PMC4_xdic);
    if(__builtin_except(isPMC4_NULL,0)) { MALLOC_FAILED}
    const bool isPMC5_NULL = (NULL==PMC5_xw)   || (NULL==PMC5_xz)   || (NULL==PMC5_xRs) ||
                             (NULL==PMC5_xchi) || (NULL==PMC5_xndt) || (NULL==PMC5_xdic);
    if(__builtin_except(isPMC5_NULL,0)) { MALLOC_FAILED}
    const bool isPMC6_NULL = (NULL==PMC6_xw)   || (NULL==PMC6_xz)   || (NULL==PMC6_xRs) ||
                             (NULL==PMC6_xchi) || (NULL==PMC6_xndt) || (NULL==PMC6_xdic);
    if(__builtin_except(isPMC6_NULL,0)) { MALLOC_FAILED}

#pragma omp parallel sections
    {
         #pragma section
	   {
              autcorf_(&in_data->pmc0[0],&len,&PMC0_acov[0],&PMC0_acor[0],&lagh,&PMC0_xmean);
	      canarmf_(&len,&lagh,&PMC0_acov[0],&PMC0_xarcoef[0],&lagh,&PMC0_xv[0],
	               &PMC0_xaic[0],&PMC0_xoaic,&PMC0_xmo,&PMC0_xparcor[0],&PMC0_xnc,
		       &PMC0_xm1[0],&PMC0_xm2[0],&PMC0_xw[0],&PMC0_xz[0],&PMC0_xRs[0],
		       &PMC0_xchi[0],&PMC0_xndt[0],&PMC0_xdic[0],&PMC0_xdicm[0],
		       &PMC0_xpo[0],&PMC0_xk,&PMC0_xb[0],&PMC0_xl,&PMC0_xa[0],&lagh,&lagh);
	   }
	 #pragma section
	   {
              autcorf_(&in_data->pmc1[0],&len,&PMC1_acov[0],&PMC1_acor[0],&lagh,&PMC1_xmean);
	      canarmf_(&len,&lagh,&PMC1_acov[0],&PMC1_xarcoef[0],&lagh,&PMC1_xv[0],
	               &PMC1_xaic[0],&PMC1_xoaic,&PMC1_xmo,&PMC1_xparcor[0],&PMC1_xnc,
		       &PMC1_xm1[0],&PMC1_xm2[0],&PMC1_xw[0],&PMC1_xz[0],&PMC1_xRs[0],
		       &PMC1_xchi[0],&PMC1_xndt[0],&PMC1_xdic[0],&PMC1_xdicm[0],
		       &PMC1_xpo[0],&PMC1_xk,&PMC1_xb[0],&PMC1_xl,&PMC1_xa[0],&lagh,&lagh);
	   }
	 #pragma section
	   {
               autcorf_(&in_data->pmc2[0],&len,&PMC2_acov[0],&PMC2_acor[0],&lagh,&PMC2_xmean);
	       canarmf_(&len,&lagh,&PMC2_acov[0],&PMC2_xarcoef[0],&lagh,&PMC2_xv[0],
	               &PMC2_xaic[0],&PMC2_xoaic,&PMC2_xmo,&PMC2_xparcor[0],&PMC2_xnc,
		       &PMC2_xm1[0],&PMC2_xm2[0],&PMC2_xw[0],&PMC2_xz[0],&PMC2_xRs[0],
		       &PMC2_xchi[0],&PMC2_xndt[0],&PMC2_xdic[0],&PMC2_xdicm[0],
		       &PMC2_xpo[0],&PMC2_xk,&PMC2_xb[0],&PMC2_xl,&PMC2_xa[0],&lagh,&lagh);
	   }
	 #pragma section
	   {
               autcorf_(&in_data->pmc3[0],&len,&PMC3_acov[0],&PMC3_acor[0],&lagh,&PMC3_xmean);
	       canarmf_(&len,&lagh,&PMC3_acov[0],&PMC3_xarcoef[0],&lagh,&PMC3_xv[0],
	               &PMC3_xaic[0],&PMC3_xoaic,&PMC3_xmo,&PMC3_xparcor[0],&PMC3_xnc,
		       &PMC3_xm1[0],&PMC3_xm2[0],&PMC3_xw[0],&PMC3_xz[0],&PMC3_xRs[0],
		       &PMC3_xchi[0],&PMC3_xndt[0],&PMC3_xdic[0],&PMC3_xdicm[0],
		       &PMC3_xpo[0],&PMC3_xk,&PMC3_xb[0],&PMC3_xl,&PMC3_xa[0],&lagh,&lagh);
	   }
	 #pragma section
	   {
               autcorf_(&in_data->pmc4[0],&len,&PMC4_acov[0],&PMC4_acor[0],&lagh,&PMC4_xmean);
	       canarmf_(&len,&lagh,&PMC4_acov[0],&PMC4_xarcoef[0],&lagh,&PMC4_xv[0],
	               &PMC4_xaic[0],&PMC4_xoaic,&PMC4_xmo,&PMC4_xparcor[0],&PMC4_xnc,
		       &PMC4_xm1[0],&PMC4_xm2[0],&PMC4_xw[0],&PMC4_xz[0],&PMC4_xRs[0],
		       &PMC4_xchi[0],&PMC4_xndt[0],&PMC4_xdic[0],&PMC4_xdicm[0],
		       &PMC4_xpo[0],&PMC4_xk,&PMC4_xb[0],&PMC4_xl,&PMC4_xa[0],&lagh,&lagh);
	   }
	 #pragma section
	   {
               autcorf_(&in_data->pmc5[0],&len,&PMC5_acov[0],&PMC5_acor[0],&lagh,&PMC5_xmean);
	       canarmf_(&len,&lagh,&PMC5_acov[0],&PMC5_xarcoef[0],&lagh,&PMC5_xv[0],
	               &PMC5_xaic[0],&PMC5_xoaic,&PM5_xmo,&PMC5_xparcor[0],&PMC5_xnc,
		       &PMC55_xm1[0],&PMC5_xm2[0],&PM5_xw[0],&PMC5_xz[0],&PMC5_xRs[0],
		       &PMC5_xchi[0],&PMC5_xndt[0],&PMC5_xdic[0],&PMC5_xdicm[0],
		       &PMC5_xpo[0],&PMC5_xk,&PMC5_xb[0],&PMC5_xl,&PMC5_xa[0],&lagh,&lagh);
	   }
	 #pragma section
	   {
             autcorf_(&in_data->pmc6[0],&len,&PMC6_acov[0],&PMC6_acor[0],&lagh,&PMC6_xmean);
	     canarmf_(&len,&lagh,&PMC6_acov[0],&PMC6_xarcoef[0],&lagh,&PMC6_xv[0],
	               &PMC6_xaic[0],&PMC6_xoaic,&PMC6_xmo,&PMC6_xparcor[0],&PMC6_xnc,
		       &PMC6_xm1[0],&PMC6_xm2[0],&PMC6_xw[0],&PMC6_xz[0],&PMC6_xRs[0],
		       &PMC6_xchi[0],&PMC6_xndt[0],&PMC6_xdic[0],&PMC6_xdicm[0],
		       &PMC6_xpo[0],&PMC6_xk,&PMC6_xb[0],&PMC6_xl,&PMC6_xa[0],&lagh,&lagh);
	   }
       } // omp parallel sections

       	fprintf(fp,"Performance Event name: %s\n",event_names[0].c_str());
	fprintf(fp,"mo=%.16f, _oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",
	        PMC0_xmo,PMC0_xoaic,PMC0_xnc,PMC0_xk,PMC0_xl);
	fprintf(fp, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         PMC0_xarcoef[i],PMC0_xv[i],PMC0_xaic[i],PMC0_xparcor[i],PMC0_xdicm[i],PMC0_xb[i],PMC0_xa[i],PMC0_xm1[i],PMC0_xm2[i],PMC0_xpo[i]);}
	fprintf(fp,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fp,"%.16f\n",PMC0_xw[i]);}
        fprintf(fp, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fp, " %.16f %.16f %.16f %d %.16f\n",
           PMC0_xz[i],PMC0_xRs[i],PMC0_xchi[i],PMC0_xndt[i],PMC0_xdic[i]);}

	fprintf(fp,"Performance Event name: %s\n",event_names[1].c_str());
	fprintf(fp,"mo=%.16f, _oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",
	        PMC1_xmo,PMC1_xoaic,PMC1_xnc,PMC1_xk,PMC1_xl);
	fprintf(fp, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         PMC1_xarcoef[i],PMC1_xv[i],PMC1_xaic[i],PMC1_xparcor[i],PMC1_xdicm[i],PMC1_xb[i],PMC1_xa[i],PMC1_xm1[i],PMC1_xm2[i],PMC1_xpo[i]);}
	fprintf(fp,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fp,"%.16f\n",PMC1_xw[i]);}
        fprintf(fp, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fp, " %.16f %.16f %.16f %d %.16f\n",
           PMC1_xz[i],PMC1_xRs[i],PMC1_xchi[i],PMC1_xndt[i],PMC1_xdic[i]);}

	fprintf(fp,"Performance Event name: %s\n",event_names[2].c_str());
	fprintf(fp,"mo=%.16f, _oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",
	        PMC2_xmo,PMC2_xoaic,PMC2_xnc,PMC2_xk,PMC2_xl);
	fprintf(fp, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         PMC2_xarcoef[i],PMC2_xv[i],PMC2_xaic[i],PMC2_xparcor[i],PMC2_xdicm[i],PMC2_xb[i],PMC2_xa[i],PMC2_xm1[i],PMC2_xm2[i],PMC2_xpo[i]);}
	fprintf(fp,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fp,"%.16f\n",PMC2_xw[i]);}
        fprintf(fp, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fp, " %.16f %.16f %.16f %d %.16f\n",
           PMC2_xz[i],PMC2_xRs[i],PMC2_xchi[i],PMC2_xndt[i],PMC2_xdic[i]);}

	fprintf(fp,"Performance Event name: %s\n",event_names[3].c_str());
	fprintf(fp,"mo=%.16f, _oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",
	        PM31_xmo,PMC3_xoaic,PMC3_xnc,PMC3_xk,PMC3_xl);
	fprintf(fp, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         PMC3_xarcoef[i],PMC3_xv[i],PMC3_xaic[i],PMC3_xparcor[i],PMC3_xdicm[i],PMC3_xb[i],PMC3_xa[i],PMC3_xm1[i],PMC3_xm2[i],PMC3_xpo[i]);}
	fprintf(fp,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fp,"%.16f\n",PMC3_xw[i]);}
        fprintf(fp, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fp, " %.16f %.16f %.16f %d %.16f\n",
           PMC3_xz[i],PMC3_xRs[i],PMC3_xchi[i],PMC3_xndt[i],PMC3_xdic[i]);}

	fprintf(fp,"Performance Event name: %s\n",event_names[4].c_str());
	fprintf(fp,"mo=%.16f, _oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",
	        PMC4_xmo,PMC4_xoaic,PMC4_xnc,PMC4_xk,PMC4_xl);
	fprintf(fp, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         PMC4_xarcoef[i],PMC4_xv[i],PMC4_xaic[i],PMC4_xparcor[i],PMC4_xdicm[i],PMC4_xb[i],PMC4_xa[i],PMC4_xm1[i],PMC4_xm2[i],PMC4_xpo[i]);}
	fprintf(fp,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fp,"%.16f\n",PMC4_xw[i]);}
        fprintf(fp, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fp, " %.16f %.16f %.16f %d %.16f\n",
           PMC4_xz[i],PMC4_xRs[i],PMC4_xchi[i],PMC4_xndt[i],PMC4_xdic[i]);}


	fprintf(fp,"Performance Event name: %s\n",event_names[5].c_str());
	fprintf(fp,"mo=%.16f, _oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",
	        PMC5_xmo,PMC5_xoaic,PMC5_xnc,PMC5_xk,PMC5_xl);
	fprintf(fp, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         PMC5_xarcoef[i],PMC5_xv[i],PMC5_xaic[i],PMC5_xparcor[i],PMC5_xdicm[i],PMC5_xb[i],PMC5_xa[i],PMC5_xm1[i],PMC5_xm2[i],PMC5_xpo[i]);}
	fprintf(fp,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fp,"%.16f\n",PMC5_xw[i]);}
        fprintf(fp, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fp, " %.16f %.16f %.16f %d %.16f\n",
           PMC5_xz[i],PMC5_xRs[i],PMC5_xchi[i],PMC5_xndt[i],PMC5_xdic[i]);}

	fprintf(fp,"Performance Event name: %s\n",event_names[6].c_str());
	fprintf(fp,"mo=%.16f, _oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",
	        PMC6_xmo,PMC6_xoaic,PMC6_xnc,PMC6_xk,PMC6_xl);
	fprintf(fp, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         PMC6_xarcoef[i],PMC6_xv[i],PMC6_xaic[i],PMC6_xparcor[i],PMC6_xdicm[i],PMC6_xb[i],PMC6_xa[i],PMC6_xm1[i],PMC6_xm2[i],PMC6_xpo[i]);}
	fprintf(fp,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fp,"%.16f\n",PMC6_xw[i]);}
        fprintf(fp, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fp, " %.16f %.16f %.16f %d %.16f\n",
           PMC6_xz[i],PMC6_xRs[i],PMC6_xchi[i],PMC6_xndt[i],PMC6_xdic[i]);}
	fclose(fp);
	_mm_free(PMC0_xw);   _mm_free(PMC0_xz); _mm_free(PMC0_xRs);
	_mm_free(PMC0_xchi); _mm_free(PMC0_xndt); _mm_free(PMC0_xdic);
	_mm_free(PMC1_xw);   _mm_free(PMC1_xz); _mm_free(PMC1_xRs);
	_mm_free(PMC1_xchi); _mm_free(PMC1_xndt); _mm_free(PMC1_xdic);
	_mm_free(PMC2_xw);   _mm_free(PMC2_xz); _mm_free(PMC2_xRs);
	_mm_free(PMC2_xchi); _mm_free(PMC2_xndt); _mm_free(PMC2_xdic);
	_mm_free(PMC3_xw);   _mm_free(PMC3_xz); _mm_free(PMC3_xRs);
	_mm_free(PMC3_xchi); _mm_free(PMC3_xndt); _mm_free(PMC3_xdic);
	_mm_free(PMC4_xw);   _mm_free(PMC4_xz); _mm_free(PMC4_xRs);
	_mm_free(PMC4_xchi); _mm_free(PMC4_xndt); _mm_free(PMC4_xdic);
	_mm_free(PMC5_xw);   _mm_free(PMC5_xz); _mm_free(PMC5_xRs);
	_mm_free(PMC5_xchi); _mm_free(PMC5_xndt); _mm_free(PMC5_xdic);
	_mm_free(PMC6_xw);   _mm_free(PMC6_xz); _mm_free(PMC6_xRs);
	_mm_free(PMC6_xchi); _mm_free(PMC6_xndt); _mm_free(PMC6_xdic);

}
			     


















































#endif /*__GMS_PMC_EVENTS_TIME_SERIES_ANALYSIS_HPP__*/
