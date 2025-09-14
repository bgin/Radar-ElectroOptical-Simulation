
#ifndef __GMS_LUT_CDF_PDF_SSE_H__
#define __GMS_LUT_CDF_PDF_SSE_H__


#include <immintrin.h>


#if defined(__GNUC__) && !defined(__INTEL_COMPILER)

/* Used by the function: gamma_log_xmm2r8*/  
 
extern const  double gamma_log_xmm2r8_c[14];

extern const  double gamma_log_xmm2r8_p1[16];

extern const  double gamma_log_xmm2r8_p2[16];

extern const  double gamma_log_xmm2r8_p4[16];

extern const  double gamma_log_xmm2r8_q1[16]; 

extern const  double gamma_log_xmm2r8_q2[16];

extern const  double gamma_log_xmm2r8_q4[16];

/////////////////////////////////////////////
/* Used by the function: gamma_log_xmm4r4*/ 

extern const  float gamma_log_xmm4r4_c[28];

extern const  float gamma_log_xmm4r4_p1[32];

extern const  float gamma_log_xmm4r4_p2[32];

extern const  float gamma_log_xmm4r4_p4[32];

extern const  float gamma_log_xmm4r4_q1[32];

extern const  float gamma_log_xmm4r4_q2[32];

extern const  float gamma_log_xmm4r4_q4[32];

/////////////////////////////////////////////     
/* Used by the function: bessel_i0_xmm2r8*/   

extern const  double bessel_i0_xmm2r8_p[30];

extern const  double bessel_i0_xmm2r8_pp[16];

extern const  double bessel_i0_xmm2r8_q[10];

extern const  double bessel_i0_xmm2r8_qq[14];

/////////////////////////////////////////////
/* Used by the function: bessel_i0_xmm4r4*/    

extern const  float bessel_i0_xmm4r4_p[60];

extern  const  float bessel_i0_xmm4r4_pp[32];

extern const  float bessel_i0_xmm4r4_q[20];

extern const  float bessel_i0_xmm4r4_qq[28];

/////////////////////////////////////////////
/* Used by the function: bessel_i1_xmm2r8*/  

extern  const  double bessel_i1_xmm2r8_p[30];

extern const  double bessel_i1_xmm2r8_pp[16];

extern const  double bessel_i1_xmm2r8_q[10];

extern const  double bessel_i1_xmm2r8_qq[12];

/////////////////////////////////////////////
/* Used by the function: bessel_i1_xmm4r4*/   

extern const  float bessel_i1_xmm4r4_p[60];

extern const  float bessel_i1_xmm4r4_pp[32];

extern const  float bessel_i1_xmm4r4_q[20];

extern const  float bessel_i1_xmm4r4_qq[24];

/////////////////////////////////////////////
/* Used by the function: normal_01_cdf_inv_xmm2r8*/   

extern const  double normal_01_cdf_inv_xmm2r8_a[16];

extern const  double normal_01_cdf_inv_xmm2r8_b[16];

extern const  double normal_01_cdf_inv_xmm2r8_c[16];

extern const  double normal_01_cdf_inv_xmm2r8_d[16];

extern const  double normal_01_cdf_inv_xmm2r8_e[16];

extern  const  double normal_01_cdf_inv_xmm2r8_f[16];

//////////////////////////////////////////////
/* Used by the function: normal_01_cdf_inv_xmm4r4*/  

extern const  float normal_01_cdf_inv_xmm4r4_a[32];

extern const  float normal_01_cdf_inv_xmm4r4_b[32];

extern const  float normal_01_cdf_inv_xmm4r4_c[32];

extern const  float normal_01_cdf_inv_xmm4r4_d[32];

extern const  float normal_01_cdf_inv_xmm4r4_e[32];

extern const  float normal_01_cdf_inv_xmm4r4_f[32];

////////////////////////////////////////////
/* Used by the function: gamma_xmm2r8*/ 

extern const  double gamma_xmm2r8_c[14];

extern const  double gamma_xmm2r8_p[16];

extern const  double gamma_xmm2r8_q[16];

////////////////////////////////////////////
/* Used by the function: gamma_xmm4r4*/ 

extern const  float gamma_xmm4r4_c[29];

extern const  float gamma_xmm4r4_p[32];

extern const  float gamma_xmm4r4_q[32];

 

#elif !defined(__GNUC__) && defined(__INTEL_COMPILER)

extern const __m128d gamma_log_xmm2r8_c[7];

extern const __m128d gamma_log_xmm2r8_p1[8];

extern const __m128d gamma_log_xmm2r8_p2[8];

extern  const __m128d gamma_log_xmm2r8_p4[8];

extern const  __m128d gamma_log_xmm2r8_q1[8];

extern const __m128d gamma_log_xmm2r8_q2[8];

extern const __m128d gamma_log_xmm2r8_q4[8];

extern const __m128 gamma_log_xmm4r4_c[7];

extern const __m128 gamma_log_xmm4r4_p1[8];

extern const __m128 gamma_log_xmm4r4_p2[8];

extern  const __m128 gamma_log_xmm4r4_p4[8];

extern const  __m128 gamma_log_xmm4r4_q1[8];

extern const __m128 gamma_log_xmm4r4_q2[8];

extern const __m128 gamma_log_xmm4r4_q4[8];

extern const  __m128d bessel_i0_xmm2r8_p[15];

extern const __m128d bessel_i0_xmm2r8_pp[8];

extern const __m128d bessel_i0_xmm2r8_q[5];

extern const __m128d bessel_i0_xmm2r8_qq[7];

extern const __m128 bessel_i0_xmm4r4_p[15];

extern  const __m128 bessel_i0_xmm4r4_pp[8];

extern const __m128 bessel_i0_xmm4r4_q[5];

extern  const __m128 bessel_i0_xmm4r4_qq[7];

extern const __m128d  bessel_i1_xmm2r8_p[15];

extern const __m128d bessel_i1_xmm2r8_pp[8];

extern const __m128d bessel_i1_xmm2r8_q[5];

extern const __m128d bessel_i1_xmm2r8_qq[6];

extern const  __m128  bessel_i1_xmm4r4_p[15];

extern const __m128 bessel_i1_xmm4r4_pp[8];

extern const __m128 bessel_i1_xmm4r4_q[5];

extern const __m128 bessel_i1_xmm4r4_qq[6];

extern  const __m128d  normal_01_cdf_inv_xmm2r8_a[8];

extern const __m128d  normal_01_cdf_inv_xmm2r8_b[8];

extern const __m128d  normal_01_cdf_inv_xmm2r8_c[8];

extern const __m128d  normal_01_cdf_inv_xmm2r8_d[8]; 

extern const __m128d  normal_01_cdf_inv_xmm2r8_e[8];

extern const __m128d  normal_01_cdf_inv_xmm2r8_f[8];

extern const __m128  normal_01_cdf_inv_xmm4r4_a[8];

extern const __m128  normal_01_cdf_inv_xmm4r4_b[8];

extern const __m128  normal_01_cdf_inv_xmm4r4_c[8];

extern const __m128d  normal_01_cdf_inv_xmm4r4_d[8];

extern const __m128  normal_01_cdf_inv_xmm4r4_e[8];

extern const __m128  normal_01_cdf_inv_xmm4r4_f[8];

extern const __m128d  gamma_xmm2r8_c[7];

extern const __m128d  gamma_xmm2r8_p[8];

extern const  __m128d  gamma_xmm2r8_q[8];

extern const __m128  gamma_xmm4r4_c[7];

extern const  __m128 gamma_xmm4r4_p[8];

extern const __m128  gamma_xmm4r4_q[8];

#endif

























#endif /*__GMS_LUT_CDF_PDF_SSE_H__*/
