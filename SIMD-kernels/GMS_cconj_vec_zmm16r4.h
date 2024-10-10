
#ifndef __GMS_CCONJ_VEC_ZMM16R4_H__
#define __GMS_CCONJ_VEC_ZMM16R4_H__ 281220221036


namespace file_version {

    const unsigned int GMS_CCONJ_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CCONJ_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CCONJ_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CCONJ_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CCONJ_VEC_ZMM16R4_MAJOR+
      100U*GMS_CCONJ_VEC_ZMM16R4_MINOR+
      10U*GMS_CCONJ_VEC_ZMM16R4_MICRO;
    const char * const GMS_CCONJ_VEC_ZMM16R4_CREATION_DATE = "28-12-2022 10:36  +00200 (WED 28 DEC 2022 GMT+2)";
    const char * const GMS_CCONJ_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CCONJ_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CCONJ_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector conjugate operations."

}

#include <cstdint>




namespace gms {

      
          namespace math {

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cconjv_zmm16r4_unroll_16x_u(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n); 

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	            void cconjv_zmm16r4_unroll_16x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n); 

                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cconjv_zmm16r4_unroll_10x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n); 

                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cconjv_zmm16r4_unroll_10x_u(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n); 

                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cconjv_zmm16r4_unroll_6x_u(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n); 

                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cconjv_zmm16r4_unroll_6x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n); 


                  
                   

        } // math


} // gms















#endif /*__GMS_CCONJ_VEC_ZMM16R4_H__*/
