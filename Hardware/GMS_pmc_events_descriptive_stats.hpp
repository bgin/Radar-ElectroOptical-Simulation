

#ifndef __GMS_PMC_EVENTS_DESCRIPTIVE_STATS_HPP__
#define __GMS_PMC_EVENTS_DESCRIPTIVE_STATS_HPP__


#include <cstdint>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "GMS_config.h"
#include "GMS_descriptive_statistics.hpp"
#include "GMS_convert_numeric_data_types.hpp"

__attribute__((hot))
__attribute__((aligned(32)))
void
pmc_events_descriptive_stats(const double * __restrict  __attribute__((aligned(64))) data,
			     const int32_t len,
			     const char   * __restrict fname,
			     const char   * __restrict event_name) {
     FILE * fp = NULL;
     if(__builtin_except(fopen(&fp,fname,"a+"),0) != 0) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
      }
      const int32_t len2   = len/2;
      constexpr float w_limit = 0.05f;
      float  __restrict  * a    = NULL;  
      float  __restrict  * df32 = NULL;  
      float   w    = 0.0f;               
      float   pw   = 0.0f;               
      int32_t ifault = -1;               
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
      a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
      if(__builtin_except(NULL==a,0)) { MALLOC_FAILED}
      df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
      if(__builtin_except(NULL==df32,0)) { MALLOC_FAILED}
      
      cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
      printf("Calling  Shapiro-Wilk normality test subroutine!!\n");
      std::sort(df32,df32+len);
      swilk(init,&df32[0],len,len,n2,&a[0],w,pw,ifault);
      if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
      fprintf(fptr,"Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",w,pw);
      if(pw<w_limit) fprintf(fp,"Warning!! -- [pw < 0.05] i.e. normality limit -- Data is not normally distributed!!\n");
      if(pw>w_limit) {
	  
	   fprintf(fp,"Descriptive Statistics for HW event: %s\n",event_name);
	   fprintf(fp,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fp,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fp,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fp,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fp,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fptr,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fp,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fp,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fp,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
	fclose(fp);
	_mm_free(df32);
	_mm_free(a);
}
			     








#endif /*__GMS_PMC_EVENTS_DESCRIPTIVE_STATS_HPP__*/
