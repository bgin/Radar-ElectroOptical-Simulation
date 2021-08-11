

#ifndef __GMS_HW_METRICS_TIME_SERIES_ANALYSIS_H__
#define __GMS_HW_METRICS_TIME_SERIES_ANALYSIS_H__

#include <cstdint>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "Timsac_iface.h"
#include "descriptive_statistics.hpp"
#include "convert_numeric_data_types.hpp"







/*
    Apply Time-Series analysis (Timsac) subroutine "CANARM".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "CANARM".
    Attempt to calculate the descritpive statistics if result
    of Wilk-Shapiro normality test allows it. 
*/
__attribute__((hot))
__attribute__((aligned(32)))
template<int32_t len>
bool
hw_perf_metrics_canarm(const double * __restrict __attribute__((aligned(64))) data,
                       const char   * __restrict fname,
		       const char   * __restrict metric_name) {
		       
    static_assert(len > 100000, "Input data length can not exceed -- **100000** elements!!");
    const int32_t lagh = (int32_t)(std::sqrtf((int32_t)len));
    const int32_t n2   = len/2; // shapiro-wilk 'a' array length.
    const std::size_t  lag2len = (std::size_t)(lagh*lagh);
    const std::size_t  lag3len = (std::size_t)lag2len*len;
    constexpr float w_limit = 0.05f;
    __attribute__((aligned(64))) double acor[lagh]    = {};
    __attribute__((aligned(64))) double acov[lagh]    = {};
    __attribute__((aligned(64))) double xarcoef[lagh] = {};
    __attribute__((aligned(64))) double xv[lagh]      = {};
    __attribute__((aligned(64))) double xaic[lagh]    = {};
    __attribute__((aligned(64))) double xparcor[lagh] = {};
    __attribute__((aligned(64))) double xdicm[lagh]   = {};
    __attribute__((aligned(64))) double xb[lagh]      = {};
    __attribute__((aligned(64))) double xa[lagh]      = {};
    __attribute__((aligned(64))) int32_t xm1[lagh]    = {};
    __attribute__((aligned(64))) int32_t xm2[lagh]    = {};
    __attribute__((aligned(64))) int32_t xpo[lagh]    = {};
    double  __restrict * xw   = NULL;
    double  __restrict * xz   = NULL;
    double  __restrict * xRs  = NULL;
    double  __restrict * xchi = NULL;
    int32_t __restrict * xndt = NULL;
    double __restrict  * xdic = NULL;
    float  __restrict  * a    = NULL; // shapiro-wilk calculated weights
    float  __restrict  * df32 = NULL;
    FILE * fptr = NULL;
    double xoaic = 0.0;
    double xmean = 0.0;
    int32_t xmo  = 0;
    int32_t xnc  = 0;
    int32_t xk   = 0;
    int32_t xl   = 0;
    float   w    = 0.0f; // swilk W-statistics
    float   pw   = 0.0f; // swilk P-value of W
    int32_t ifault = -1; // swilk error indicator
    float srsd = 0.0f;
    float svar = 0.0f;
    float skew = 0.0f;
    float kurt = 0.0f;
    float autocor = 0.0f;
   
    float xmid = 0.0f;
    float xmean = 0.0f;
    float xmidm = 0.0f;
    float xmed = 0.0f;
    float smin = 0.0f;
    float smax = 0.0f;
 
    float xrange = 0.0f;
    float xsd = 0.0f;
    float xrelsd = 0.0f;
    float xvar = 0.0f;
    
    const bool init = false; // swilk init argument.
    // OpenMP multithreaded calls to _mm_malloc (using parallel sections)
#pragma omp parallel sections
    {
            #pragma section
	    {
               xw = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	      
	    }
	    #pragma section
	    {
               xz = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	    
	    }
	    #pragma section
	    {
               xRs = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	     
	    }
	    #pragma section
	    {
              xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	     
	    }
	    #pragma section
	    {
               xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	      
	    }
	    #pragma section
	    {
               xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	    }
	    #pragma section
	    {
               a    = reinterpret_cast<float*>(_mm_malloc((std::size_t)n2*sizeof(float),64));
	    }
	    #pragma section
	    {
               df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)n*sizeof(float),64));
	    }
    }

       // Single thread checks the returned pointers!!
        const bool is_ptr_null = NULL==a    || NULL==xdic || NULL==xndt ||
	                         NULL==xchi || NULL==xRs  || NULL==xz   ||
				 NULL==xw   || NULL==df32;
	if(is_ptr_null) {MALLOC_FAILED}
	
	printf("Calling: AUTCORF TIMSAC subroutine!!\n");
	autcorf_(&data[0],&len,&acov[0],&acor[0],&lagh,&xmean);
	printf("Calling: CANARMF TIMSAC subroutine!!\n");
	canarmf_(&len,&lagh,&acov[0],&xarcoef[0],&lagh,&xv[0],&xaic[0],&xoaic[0],
	         &xmo,&xparcor[0],&xnc,&xm1[0],&xm2[0],&xw[0],&xz[0],&xRs[0],
		 &xchi[0],&xndt[0],&xdic[0],&xdicm[0],&xpo[0],&xk,&xb[0],&xl,
		 &xa[0],&lagh,&lagh);
	printf("Starting File I/O operations!!\n");
	fptr = fopen(fname,"a+");
	if(NULL==fptr) {
            printf("File open error: %s\n",fname1);
	    std::exit(EXIT_FAILURE);
	}
	fprintf(fptr,"HW Metric name: %s\n",metric_name);
	fprintf(fptr,"mo=%.16f, oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",xmo,xoaic,xnc,xk,xl);
	fprintf(fptr, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fptr,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         xarcoef[i],xv[i],xaic[i],xparcor[i],xdicm[i],xb[i],xa[i],xm1[i],xm2[i],xpo[i]);}
	fprintf(fptr,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fptr,"%.16f\n",xw[i]);}
        fprintf(fptr, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fptr, " %.16f %.16f %.16f %d %.16f\n",
           xz[i],xRs[i],xchi[i],xndt[i],xdic[i]);}
	fprintf(fptr, "End of CANARMF results dump\n");
         // Sort a samples arrays in ascending order
	//std::sort(data,data+len);
	cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
	printf("Calling  Shapiro-Wilk normality test subroutine!!\n");
        std::sort(a,a+len);
        swilk(init,&df32[0],len,len,n2,&a[0],w,pw,ifault);
	if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
	fprintf(fptr,"Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",w,pw);
	if(pw<w_limit) fprintf(fptr,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
        if(pw>w_limit) {
	  
	   fprintf(fptr,"Descriptive Statistics calculations!!\n");
	   fprintf(fptr,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fptr,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fptr,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fptr,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fptr,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fptr,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fptr,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fptr,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fptr,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
	fclose(fptr);
	_mm_free(df32);
	_mm_free(a);
	_mm_free(xdic);
	_mm_free(xndt);
	_mm_free(xchi);
	_mm_free(xRs);
	_mm_free(xz);
	_mm_free(xw);
}


/*
    Apply Time-Series analysis (Timsac) subroutine "MULCOR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "MULCOR".
    No descriptive statistics computations for this function.
*/
bool
hw_perf_metrics_mulcor(const double * __restrict, //multidimensional data
                       const char   * __restrict,
		       const int32_t, // Number of dimensions
		       const int32_t) __attribute__((hot))
				                __attribute__((aligned(32)));


/*
    Apply Time-Series analysis (Timsac) subroutine "MULSPE".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "MULSPE".
    No descriptive statistics computations for this function.
*/
bool
hw_perf_metrics_mulspe(const double * __restrict, // Multidimensional data
                                 const char   * __restrict,
				 const int32_t,  // Number of dimensions
				 const int32_t) __attribute__((hot))
				                __attribute__((aligned(32)));


/*
    Apply Time-Series analysis (Timsac) subroutine "UNIMAR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "UNIMAR".
    Attempt to calculate the descritpive statistics if result
    of Wilk-Shapiro normality test allows it. 
*/
bool
hw_perf_metrics_unimar(const double * __restrict,
                                 const char   * __restrict,
				 const int32_t,
				 const bool)    __attribute__((hot))
				                __attribute__((aligned(32)));



/*
    Apply Time-Series analysis (Timsac) subroutine "UNIBAR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "UNIBAR".
   
*/					  
bool
hw_perf_metrics_unibar(const double * __restrict,
                                 const char   * __restrict,
				 const int32_t,
				 const bool)    __attribute__((hot))
				                __attribute__((aligned(32)));


/*
    Apply Time-Series analysis (Timsac) subroutine "EXSAR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "EXSAR".
   
*/
bool
hw_perf_metrics_exsar( const double * __restrict,
                                 const char   * __restrict,
				 const int32_t,
				 const bool)    __attribute__((hot))
				                __attribute__((aligned(32)));

/*
    Apply Time-Series analysis (Timsac) subroutine "BISPEC".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "BISPEC".
  
*/
bool
hw_perf_metrics_bispec(const double * __restrict,
                                 const char   * __restrict,
				 const int32_t,
				 const bool)    __attribute__((hot))
				                __attribute__((aligned(32)));


/*
    Apply Time-Series analysis (Timsac) subroutine "THIRMO".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "THIRMO".
   
*/
bool
hw_perf_metrics_thirmo(const double * __restrict,
                                 const char   * __restrict,
				 const int32_t,
				 const bool)    __attribute__((hot))
				                __attribute__((aligned(32)));


/*
    Apply Time-Series analysis (Timsac) subroutine "AUTOCOR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "AUTOCOR".
   
*/
bool
hw_perf_metrics_autocor(const double * __restrict,
                                  const char   * __restrict,
				  const int32_t,
				  const bool)    __attribute__((hot))
				                __attribute__((aligned(32)));

#endif /*__GMS_HW_METRICS_TIME_SERIES_ANALYSIS_H__*/
