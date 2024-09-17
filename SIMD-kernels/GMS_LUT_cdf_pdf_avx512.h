
/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#ifndef __GMS_LUT_CDF_PDF_AVX512_H__
#define __GMS_LUT_CDF_PDF_AVX512_H__



#include <immintrin.h>



#if defined(__GNUC__) && !defined(__INTEL_COMPILER)

 /* Used by the function: gamma_log_zmm8r8*/               

extern const  double gamma_log_zmm8r8_c[56];
               
 /* Used by the function: gamma_log_zmm8r8*/               

extern const  double gamma_log_zmm8r8_p1[64];
        
/* Used by the function: gamma_log_zmm8r8*/               

extern const  double gamma_log_zmm8r8_p2[64]; 
        
/* Used by the function: gamma_log_zmm8r8*/               
extern const  double gamma_log_zmm8r8_p4[64];
        
/* Used by the function: gamma_log_zmm8r8*/               
extern const  double gamma_log_zmm8r8_q1[64];
        
/* Used by the function: gamma_log_zmm8r8*/               
extern const  double gamma_log_zmm8r8_q2[64];
        
/* Used by the function: gamma_log_zmm8r8*/               
extern  const  double gamma_log_zmm8r8_q4[64];

 /* Used by the function: gamma_log_zmm16r4*/               
extern const  float gamma_log_zmm16r4_c[112];
        
 /* Used by the function: gamma_log_zmm16r4*/               
extern const  float gamma_log_zmm16r4_p1[128];

/* Used by the function: gamma_log_zmm16r4*/               
extern const  float gamma_log_zmm16r4_p2[128];
        
/* Used by the function: gamma_log_zmm16r4*/               
extern const  float gamma_log_zmm16r4_p4[128];

/* Used by the function: gamma_log_zmm16r4*/               
extern const  float gamma_log_zmm16r4_q1[128];
        
/* Used by the function: gamma_log_zmm16r4*/               
extern const  float gamma_log_zmm16r4_q2[128];
        
/* Used by the function: gamma_log_zmm16r4*/               
extern const  float gamma_log_zmm16r4_q4[128];

////////////////////////////////////////////////////////////////////////////////////       
/* Used by the function: bessel_i0_zmm8r8*/      
         
extern const  double bessel_i0_zmm8r8_p[120];
       
/* Used by the function: bessel_i0_zmm8r8*/               
extern const  double bessel_i0_zmm8r8_pp[64];
       
/* Used by the function: bessel_i0_zmm8r8*/               
extern const  double bessel_i0_zmm8r8_q[40];
       
/* Used by the function: bessel_i0_zmm8r8*/               
extern  const  double bessel_i0_zmm8r8_qq[56];
       
/* Used by the function: bessel_i0_zmm16r4*/               
extern const  float bessel_i0_zmm16r4_p[240];
       
/* Used by the function: bessel_i0_zmm16r4*/               
extern const  float bessel_i0_zmm16r4_pp[128];

/* Used by the function: bessel_i0_zmm16r4*/               
extern const  float bessel_i0_zmm16r4_q[80];
/* Used by the function: bessel_i0_zmm8r8*/               

extern const  float bessel_i0_zmm16r4_qq[112];
       
/* Used by the function: bessel_i1_zmm8r8*/               

extern const  double bessel_i1_zmm8r8_p[120]; 
/* Used by the function: bessel_i1_zmm8r8*/               

extern const  double bessel_i1_zmm8r8_pp[64];
       
/* Used by the function: bessel_i1_zmm8r8*/               

extern const  double bessel_i1_zmm8r8_q[40];
       
/* Used by the function: bessel_i1_zmm8r8*/               

extern const  double bessel_i1_zmm8r8_qq[48];
         
/* Used by the function: bessel_i1_zmm16r4*/               

extern const  float bessel_i1_zmm16r4_p[240];
       
/* Used by the function: bessel_i1_zmm16r4*/               

extern const  float bessel_i1_zmm16r4_pp[128];
/* Used by the function: bessel_i1_zmm16r4*/               

extern const  float bessel_i1_zmm16r4_q[80];
       
/* Used by the function: bessel_i1_zmm16r4*/               

extern const  float bessel_i1_zmm16r4_qq[96];
/* Used by the function: normal_01_cdf_inv_zmm8r8*/               

extern const  double normal_01_cdf_inv_zmm8r8_a[64];
/* Used by the function: normal_01_cdf_inv_zmm8r8*/               

extern const  double normal_01_cdf_inv_zmm8r8_b[64];
         
/* Used by the function: normal_01_cdf_inv_zmm8r8*/               

extern const  double normal_01_cdf_inv_zmm8r8_c[64];
         
/* Used by the function: normal_01_cdf_inv_zmm8r8*/               

extern const  double normal_01_cdf_inv_zmm8r8_d[64];
/* Used by the function: normal_01_cdf_inv_zmm8r8*/               

extern const  double normal_01_cdf_inv_zmm8r8_e[64];
         
/* Used by the function: normal_01_cdf_inv_zmm8r8*/               

extern const  double normal_01_cdf_inv_zmm8r8_f[64];

//////////////////////////////////////////////////////////////////////////////////////////////////
/* Used by the function: normal_01_cdf_inv_zmm16r4*/               

extern const  float normal_01_cdf_inv_zmm16r4_a[128];
         
/* Used by the function: normal_01_cdf_inv_zmm16r4*/               

extern const  float normal_01_cdf_inv_zmm16r4_b[128];
         
/* Used by the function: normal_01_cdf_inv_zmm16r4*/               

extern const  float normal_01_cdf_inv_zmm16r4_c[128];
         
/* Used by the function: normal_01_cdf_inv_zmm16r4*/               

extern const  float normal_01_cdf_inv_zmm16r4_d[128];
         
/* Used by the function: normal_01_cdf_inv_zmm16r4*/               

extern const  float normal_01_cdf_inv_zmm16r4_e[128];
/* Used by the function: normal_01_cdf_inv_zmm16r4*/               

extern const  float normal_01_cdf_inv_zmm16r4_f[128];
///////////////////////////////////////////////////////////////////////////
/* Used by the function: gamma_zmm8r8*/               

extern const  double gamma_zmm8r8_c[58];
         

extern const  double gamma_zmm8r8_p[64];
         

extern const  double gamma_zmm8r8_q[64];
         
/* Used by the function: gamma_zmm16r4*/               
extern const  float gamma_zmm16r4_c[116]
         
extern const  float gamma_zmm16r4_p[128];
         
extern const  float gamma_zmm16r4_q[128];






         
#elif !defined(__GNUC__) && defined(__INTEL_COMPILER)


extern const  __m512d gamma_log_zmm8r8_c[7];
        

extern const   __m512d gamma_log_zmm8r8_p1[8]; 
                                    
extern const   __m512d gamma_log_zmm8r8_p2[8]; 
       

extern const   __m512d gamma_log_zmm8r8_p4[8];
       

extern const   __m512d gamma_log_zmm8r8_q1[8];
       

extern const   __m512d gamma_log_zmm8r8_q2[8]; 
       

extern const   __m512d gamma_log_zmm8r8_q4[8]; 
       

extern const   __m512 gamma_log_zmm16r4_c[7]; 
        

extern const   __m512 gamma_log_zmm16r4_p1[8]; 
        
                                    
extern const   __m512 gamma_log_zmm16r4_p2[8]; 
       

extern const   __m512 gamma_log_zmm16r4_p4[8]; 
       

extern const   __m512 gamma_log_zmm16r4_q1[8];
       

extern const   __m512 gamma_log_zmm16r4_q2[8];
       

extern const  __m512 gamma_log_zmm16r4_q4[8];
/////////////////////////////////////////////////////////////////////////////////////////////////       

extern const   __m512d bessel_i0_zmm8r8_p[15]; 
                                     
extern const   __m512d bessel_i0_zmm8r8_pp[8]; 

extern const   __m512d bessel_i0_zmm8r8_q[5]; 
                                           
extern const   __m512d bessel_i0_zmm8r8_qq[7]; 

extern const   __m512 bessel_i0_zmm16r4_p[15];
 
                                   
extern const  __m512 bessel_i0_zmm16r4_pp[8]; 

extern const   __m512 bessel_i0_zmm16r4_q[5];  
                                           
extern const   __m512 bessel_i0_zmm16r4_qq[7]; 
       

extern const   __m512d  bessel_i1_zmm8r8_p[15]; 
       

extern const   __m512d bessel_i1_zmm8r8_pp[8];  
       

extern const   __m512d bessel_i1_zmm8r8_q[5];  
       

extern const   __m512d bessel_i1_zmm8r8_qq[6]; 

extern const   __m512  bessel_i1_zmm16r4_p[15]; 

extern const   __m512 bessel_i1_zmm16r4_pp[8]; 
       

extern const   __m512 bessel_i1_zmm16r4_q[5];  

extern const   __m512 bessel_i1_zmm16r4_qq[6];  
///////////////////////////////////////////////////////////////////////////

extern const   __m512d  normal_01_cdf_inv_zmm8r8_a[8]; 
        

extern const   __m512d  normal_01_cdf_inv_zmm8r8_b[8];
        

extern const   __m512d  normal_01_cdf_inv_zmm8r8_c[8];
        

extern const   __m512d  normal_01_cdf_inv_zmm8r8_d[8]; 

extern const   __m512d  normal_01_cdf_inv_zmm8r8_e[8];
        

extern const   __m512d  normal_01_cdf_inv_zmm8r8_f[8]; 
      

extern const   __m512  normal_01_cdf_inv_zmm16r4_a[8]; 
        

extern const   __m512  normal_01_cdf_inv_zmm16r4_b[8]; 
        

extern const   __m512  normal_01_cdf_inv_zmm16r4_c[8]; 

extern const   __m512d  normal_01_cdf_inv_zmm16r4_d[8]; 
        

extern const   __m512  normal_01_cdf_inv_zmm16r4_e[8];
        

extern const   __m512  normal_01_cdf_inv_zmm16r4_f[8];  
/////////////////////////////////////////////////////////////////////////////////


extern const  __m512d  gamma_zmm8r8_c[7];
        

extern const   __m512d  gamma_zmm8r8_p[8]; 
        

extern const   __m512d  gamma_zmm8r8_q[8]; 
        

extern const   __m512  gamma_zmm16r4_c[7]; 
        

extern const    __m512 gamma_zmm16r4_p[8];  
        
         
extern const   __m512  gamma_zmm16r4_q[8];  




#endif




































#endif /*__GMS_LUT_CDF_PDF_AVX512_H__*/
