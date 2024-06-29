

#ifndef __GMS_FAST_APPROX_ZMM16R4_HPP__
#define __GMS_FAST_APPROX_ZMM16R4_HPP__ 290620240949


namespace file_info {

 const unsigned int GMS_FAST_APPROX_ZMM16R4_MAJOR = 1U;
 const unsigned int GMS_FAST_APPROX_ZMM16R4_MINOR = 0U;
 const unsigned int GMS_FAST_APPROX_ZMM16R4_MICRO = 0U;
 const unsigned int GMS_FAST_APPROX_ZMM16R4_FULLVER =
  1000U*GMS_FAST_APPROX_ZMM16R4_MAJOR+100U*GMS_FAST_APPROX_ZMM16R4_MINOR+10U*GMS_FAST_APPROX_ZMM16R4_MICRO;
 const char * const GMS_FAST_APPROX_ZMM16R4_CREATION_DATE = "29-06-2024 09:49 +00200 (SAT 29 JUN 2024 09:49 GMT+2)";
 const char * const GMS_FAST_APPROX_ZMM16R4_BUILD_DATE    = __DATE__ " " __TIME__;
 const char * const GMS_FAST_APPROX_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const GMS_FAST_APPROX_ZMM16R4_SYNOPSIS      = "Fast approximation of some functions, based on the book: Approximation for Digital Computers, by Cecil Hastings.";


}

#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

         namespace math {
           
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 fast_log10_zmm16r4(const __m512 x) {
		            
		            const __m512 C1 = _mm512_set1_ps(0.86304f);
		            const __m512 C3 = _mm512_set1_ps(0.36415f);
		            const __m512 ONE= _mm512_set1_ps(1.0f);
		            register __m512 rat,cubed;
		            register __m512 t0,t1;
		            register __m512 log10;
		            t0    = _mm512_sub_ps(x,one);
		            t1    = _mm512_add_ps(x,one);
		            rat   = _mm512_div_ps(t0,t1);
		            cubed = _mm512_mul_ps(rat,_mm512_mul_ps(rat,rat));
		            log10  = _mm512_fmadd_ps(C1,rat,_mm512_mul_ps(C3,cubed));
		            return (log10); 
		      }
       
     }

}



#endif /*__GMS_FAST_APPROX_ZMM16R4_HPP__*/

