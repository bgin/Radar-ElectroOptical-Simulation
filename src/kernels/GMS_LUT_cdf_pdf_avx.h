
#ifndef __GMS_LUT_CDF_PDF_AVX_H__
#define __GMS_LUT_CDF_PDF_AVX_H__



#include <immintrin.h>


#if defined(__GNUC__) && !defined(__INTEL_COMPILER)

/////////////////////////////////////////////////////
 /* Used by the function: gamma_log_ymm4r8*/  
 
extern const  double gamma_log_ymm4r8_c[28];

extern const  double gamma_log_ymm4r8_p1[32];

extern  const  double gamma_log_ymm4r8_p2[32];

extern const  double gamma_log_ymm4r8_p4[32];

extern const  double gamma_log_ymm4r8_q1[32];

extern const  double gamma_log_ymm4r8_q2[32];

extern const  double gamma_log_ymm4r8_q4[32];

////////////////////////////////////////////////////
/* Used by the function: gamma_log_ymm8r4*/    

extern const  float gamma_log_ymm8r4_c[56];

extern const  float gamma_log_ymm8r4_p1[64];

extern const  float gamma_log_ymm8r4_p2[64];

extern const  float gamma_log_ymm8r4_p4[64];

extern const  float gamma_log_ymm8r4_q1[64];

extern const  float gamma_log_ymm8r4_q2[64];

extern  const  float gamma_log_ymm8r4_q4[64];

////////////////////////////////////////////////////
/* Used by the function: bessel_i0_ymm4r8*/  

extern const  double bessel_i0_ymm4r8_p[60];

extern const  double bessel_i0_ymm4r8_pp[32];

extern const  double bessel_i0_ymm4r8_q[20];

extern const  double bessel_i0_ymm4r8_qq[28];

///////////////////////////////////////////////////
/* Used by the function: bessel_i0_ymm8r4*/  

extern const  float bessel_i0_ymm8r4_p[120];

extern const  float bessel_i0_ymm8r4_pp[64];

extern const  float bessel_i0_ymm8r4_q[40];

extern  const  float bessel_i0_ymm8r4_qq[56];

///////////////////////////////////////////////////
/* Used by the function: bessel_i1_ymm4r8*/   

extern const  double bessel_i1_ymm4r8_p[60];

extern const  double bessel_i1_ymm4r8_pp[32];

extern const  double bessel_i1_ymm4r8_q[20];

extern const  double bessel_i1_ymm4r8_qq[24];

///////////////////////////////////////////////////
/* Used by the function: bessel_i1_ymm8r4*/    

extern const  float bessel_i1_ymm8r4_p[120];

extern const  float bessel_i1_ymm8r4_pp[64];

extern const  float bessel_i1_ymm8r4_q[40];

extern const  float bessel_i1_ymm8r4_qq[48];

//////////////////////////////////////////////////
/* Used by the function: normal_01_cdf_inv_ymm4r8*/

extern const  double normal_01_cdf_inv_ymm4r8_a[32];

extern const  double normal_01_cdf_inv_ymm4r8_b[32];

extern const  double normal_01_cdf_inv_ymm4r8_c[32];

extern const  double normal_01_cdf_inv_ymm4r8_d[32];

extern const  double normal_01_cdf_inv_ymm4r8_e[32];

extern const  double normal_01_cdf_inv_ymm4r8_f[32];

/////////////////////////////////////////////////////
/* Used by the function: normal_01_cdf_inv_ymm8r4*/  

extern const  float normal_01_cdf_inv_ymm8r4_a[64];

extern const  float normal_01_cdf_inv_ymm8r4_b[64];

extern const  float normal_01_cdf_inv_ymm8r4_c[64];

extern const  float normal_01_cdf_inv_ymm8r4_d[64];

extern const  float normal_01_cdf_inv_ymm8r4_e[64];

extern const  float normal_01_cdf_inv_ymm8r4_f[64];

////////////////////////////////////////////////////
/* Used by the function: gamma_ymm4r8*/ 

extern const  double gamma_ymm4r8_c[28];

extern const  double gamma_ymm4r8_p[32];

extern const  double gamma_ymm4r8_q[32];

///////////////////////////////////////////////////
/* Used by the function: gamma_ymm8r4*/      

extern const  float gamma_ymm8r4_c[56];

extern const  float gamma_ymm8r4_p[64];

extern const  float gamma_ymm8r4_q[64];



#elif !defined(__GNUC__) && defined(__INTEL_COMPILER)

extern const __m256d gamma_log_ymm4r8_c[7];

extern const __m256d gamma_log_ymm4r8_p1[8];

extern  const __m256d gamma_log_ymm4r8_p2[8];

extern const __m256d gamma_log_ymm4r8_p4[8];

extern const __m256d gamma_log_ymm4r8_q1[8];

extern const __m256d gamma_log_ymm4r8_q2[8];

extern const __m256d gamma_log_ymm4r8_q4[8];

extern const __m256 gamma_log_ymm8r4_c[7];

extern const __m256 gamma_log_ymm8r4_p1[8];

extern const __m256 gamma_log_ymm8r4_p2[8];

extern const __m256 gamma_log_ymm8r4_p4[8];

extern const __m256 gamma_log_ymm8r4_q1[8];

extern const __m256 gamma_log_ymm8r4_q2[8];

extern const __m256 gamma_log_ymm8r4_q4[8];

extern const __m256d bessel_i0_ymm4r8_p[15];

extern const __m256d bessel_i0_ymm4r8_pp[8];

extern const __m256d bessel_i0_ymm4r8_q[5];

extern const __m256d bessel_i0_ymm4r8_qq[7];

extern const __m256 bessel_i0_ymm8r4_p[15];

extern const __m256 bessel_i0_ymm8r4_pp[8];

extern const __m256 bessel_i0_ymm8r4_q[5];

extern const __m256 bessel_i0_ymm8r4_qq[7];

extern const __m256d  bessel_i1_ymm4r8_p[15];

extern const __m256d bessel_i1_ymm4r8_pp[8];

extern const __m256d bessel_i1_ymm4r8_q[5];

extern const __m256d bessel_i1_ymm4r8_qq[6];

extern const __m256  bessel_i1_ymm8r4_p[15];

extern const __m256 bessel_i1_ymm8r4_pp[8];

extern const __m256 bessel_i1_ymm8r4_q[5];

extern const __m256 bessel_i1_ymm8r4_qq[6];

extern const __m256d  normal_01_cdf_inv_ymm4r8_a[8];

extern const __m256d  normal_01_cdf_inv_ymm4r8_b[8];

extern const __m256d  normal_01_cdf_inv_ymm4r8_c[8];

extern const __m256d  normal_01_cdf_inv_ymm4r8_d[8];

extern const __m256d  normal_01_cdf_inv_ymm4r8_e[8];

extern const __m256d  normal_01_cdf_inv_ymm4r8_f[8];

extern const __m256  normal_01_cdf_inv_ymm8r4_a[8];

extern  const __m256  normal_01_cdf_inv_ymm8r4_b[8];

extern  const __m256  normal_01_cdf_inv_ymm8r4_c[8];

extern const __m256d  normal_01_cdf_inv_ymm8r4_d[8];

extern const __m256  normal_01_cdf_inv_ymm8r4_e[8];

extern const __m256  normal_01_cdf_inv_ymm8r4_f[8];

extern const __m256d  gamma_ymm4r8_c[7];

extern const __m256d  gamma_ymm4r8_p[8];

extern const __m256d  gamma_ymm4r8_q[8];

extern const __m256  gamma_ymm8r4_c[7];

extern const  __m256 gamma_ymm8r4_p[8];

extern  const __m256  gamma_ymm8r4_q[8];

#endif

















































#endif /*__GMS_LUT_CDF_PDF_AVX_H__*/
