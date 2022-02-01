

#include <cmath>
#include <omp.h>
#include "GMS_laplace_inv2D.h"
#include "GMS_cephes.h"
#include "GMS_indices.h"


bool
gms::math::
gaver_stehfest_r8_1(double(*Fx)(double s1, double s2),
                    double * __restrict __ATTR_ALIGN__(64) y2d,
	            double * __restrict __ATTR_ALIGN__(64) al,
	            double * __restrict __ATTR_ALIGN__(64) om,
	            double * __restrict __ATTR_ALIGN__(64) t,
	            const double c1,
	            const double c2,
	            const int32_t nappr,
	            const int32_t ncoeff,
	            int32_t & niters) {

       if(__builtin_expect(ncoeff<1,0) &&
          __builtin_expect(nappr<1,0)) {
          return (false);
       }
       double ti;
       double tj;
       double sum1;
       double sum2;
       double t0;
       int32_t i,j,ii,jj;
       int32_t n;
       n = ncoeff+ncoeff;
#if defined(__INTEL_COMPILER) || defined(__ICC)
       __assume_aligned(y2d,64);
       __assume_aligned(al,64);
       __assume_aligned(om,64);
       __assume_aligned(t,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
      y2d = (double*)__builtin_assume_aligned(y2d,64);
      al  = (double*)__builtin_assume_aligned(al,64);
      om  = (double*)__builtin_assume_aligned(om,64);
      t   = (double*)__builtin_assume_aligned(t,64);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
      for(i = 1; i < nappr; ++i) {
             _mm_prefetch((const char*)&t[i+16],_MM_HINT_T0);
          for(j = 1; j < nappr; ++j) {
	      _mm_prefetch((const char*)&t[j+16],_MM_HINT_T0);
              sum1 = 0.0;
	      ti = t[i];
	      const double ei = std::exp(ti*c1);
	      tj = t[j];
	      const double ej = std::exp(tj*c2);
	      for(ii = 1; ii < n; ++ii) {
                  sum2 = 0.0;
		  const double omii = om[ii]
		  t0   = c1+al[ii]/ti;
#pragma omp simd reduction(+:sum2)  aligned(om:64) \
                  aligned(al:64) linear(jj:1) 
		  for(jj = 1; jj < n; ++jj) {
                      sum2 += om[jj]*Fx(t0,c2+al[jj]/tj);
		  }
		  sum1 += omii*sum2;
	      }
	      sum1 = sum1/ti/tj;
	      y2d[Ix2D(i,nappr,j)] = ei*ej*sum1;
	  }
      }
      niters = nappr*nappr*n*n;
      return (true);
}



bool
gms::math::
gaver_stehfest_r4_1(float(*Fx)(float s1, float s2),
                    float * __restrict __ATTR_ALIGN__(64) y2d,
	            float * __restrict __ATTR_ALIGN__(64) al,
	            float * __restrict __ATTR_ALIGN__(64) om,
	            float * __restrict __ATTR_ALIGN__(64) t,
	            const float c1,
	            const float c2,
	            const int32_t nappr,
	            const int32_t ncoeff,
	            int32_t & niters) {

       if(__builtin_expect(ncoeff<1,0) &&
          __builtin_expect(nappr<1,0)) {
          return (false);
       }
       float ti;
       float tj;
       float sum1;
       float sum2;
       float t0;
       int32_t i,j,ii,jj;
       int32_t n;
       n = ncoeff+ncoeff;
#if defined(__INTEL_COMPILER) || defined(__ICC)
       __assume_aligned(y2d,64);
       __assume_aligned(al,64);
       __assume_aligned(om,64);
       __assume_aligned(t,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
      y2d = (float*)__builtin_assume_aligned(y2d,64);
      al  = (float*)__builtin_assume_aligned(al,64);
      om  = (float*)__builtin_assume_aligned(om,64);
      t   = (float*)__builtin_assume_aligned(t,64);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
      for(i = 1; i < nappr; ++i) {
             _mm_prefetch((const char*)&t[i+16],_MM_HINT_T0);
          for(j = 1; j < nappr; ++j) {
	      _mm_prefetch((const char*)&t[j+16],_MM_HINT_T0);
              sum1 = 0.0;
	      ti = t[i];
	      const float ei = cephes_expf(ti*c1);
	      tj = t[j];
	      const float ej = cephes_expf(tj*c2);
	      for(ii = 1; ii < n; ++ii) {
                  sum2 = 0.0;
		  const float omii = om[ii]
		  t0   = c1+al[ii]/ti;
#pragma omp simd reduction(+:sum2)  aligned(om:64) \
                  aligned(al:64) linear(jj:1) 
		  for(jj = 1; jj < n; ++jj) {
                      sum2 += om[jj]*Fx(t0,c2+al[jj]/tj);
		  }
		  sum1 += omii*sum2;
	      }
	      sum1 = sum1/ti/tj;
	      y2d[Ix2D(i,nappr,j)] = ei*ej*sum1;
	  }
      }
      niters = nappr*nappr*n*n;
      return (true);
}

/*
    Adapted from H. Weber implementation.
*/
bool
gms::math::
gaver_stehfest_coff_r8_1(int32_t ncoeff, 
			 double * __restrict __ATTR_ALIGN__(64) al, 
			 double * __restrict __ATTR_ALIGN__(64) om) {

     if(__builtin_expect(ncoeff<1,0)) {
        return (false);
     }
     double * __restrict x = NULL;
     double * __restrict y = NULL;
     constexpr double ln2 = 0.6931471805599453094172;
     int32_t i,j,k;
     int32_t nx,sx,n;
     n = coeff+coeff;
     nh = n/2;
     const std::size_t len1 = (std::size_t)(n+2);
     const std::size_t len2 = (std::size_t)(n/2-1+2);
     x = (double*)_mm_malloc(len1,64);
     if(__builtin_expect(NULL==x,0) &&
        __builtin_expec(len1!=0,0)) {
        return (false);
      }
     y = (double*)_mm_malloc(len2,64);
     if(__builtin_expect(NULL==y,0) &&
        __builtin_expec(len2!=0,0)) {
        return (false);
      }
     x[0] = 1.0;
#if defined(__INTEL_COMPILER) || defined(__ICC)
     __assume_aligned(x,64);
     __assume_aligned(y,64);
     __assume_aligned(al,64);
     __assume_aligned9om,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     x = (double*)__builtin_assume_aligned(x,64);
     y = (double*)__builtin_assume_aligned(y,64);
     om= (double*)__builtin_assume_aligned(om,64);
     al= (double*)__builtin_assume_aligned(al,64);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
#pragma omp simd linear(i:1) aligned(x:64)
     for(i = 1; i < n; ++i) { x[i] = x[i-1]*(double)i;}
     y[i] = 2.0/x[nh-1];
     for(i = 2; i < nh; ++i) {
         const double xi = x[i];
	 const double xi1= x[i-1];
	 const double xnhi=x[nh-i];
	 const double den = xi*xi1*xnhi;
         y[i] = std::pow(i,nh)*x[2*i]/den;
     }
     sx = 2*sign<int32_t>(nh-(nh/2)*2)-1);
     for(i = 1; i < n; ++i) {
         int32_t ix;
         om[i] = 0.0;
	 ix = (i<nh)?i:nh;
	 for(k = (i+1)/2; k < ix; ++k) {
             om[i] += y[k]/x[i-1]*x[2*k-1];
	 }
	 om[i] *= (double)sx;
	 sx = -sx;
     }
#pragma omp simd linear(i:1) aligned(al:64,om)      
     for(i = 1; i < n; ++i) {
         al[i] = ln2*(double)i;
	 om[i] *= om[i]*ln2;
     }
     _mm_free(y);
     _mm_free(x);
     return (true);
}
