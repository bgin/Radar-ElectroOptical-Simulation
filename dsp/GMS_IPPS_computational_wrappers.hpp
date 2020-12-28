

#ifndef __GMS_IPPS_COMPUTATIONAL_WRAPPERS_HPP__
#define __GMS_IPPS_COMPUTATIONAL_WRAPPERS_HPP__

/*
    Work functions for IPPS computational functions.
*/




#include <cstdint>
#include <ipps.h>
#include "GMS_config.h"

















  




 










#include <immintrin.h>

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
bool      vecf32_fill_ippsRandUniform_32f(Ipp32f * __restrict pDst, int32_t len,
                                          Ipp32f low, Ipp32f high) {

      IppsRandUniState_32f * __restrict pRandStateObj = NULL;
      IppStatus stat;
      int32_t sizeRndObj;
      uint32_t seed; 
      stat = ippsRandUniformGetSize_32f(&sizeRndObj);
      if(stat != ippStsNoErr) { goto Failed;}
      pRandStateObj = (IppsRandUniState_32f*)ippsMalloc_32f(sizeRndObj);
      if(NULL==pRandStateObj && sizeRndObj != 0) {goto Failed;}
      _rdseed32_step(&seed);
      if(0==seed) { seed = (uint32_t)__rdtsc();}
      stat = ippsRandUniformInit_32f(pRandStateObj,low,high,seed);
      if(stat != ippStsNoErr) { goto Failed;}
      stat = ippsRandUniform_32f(pDst,len,pRandStateObj);
      if(stat != ippStsNoErr) { goto Failed;}
      ippsFree(pRandStateObj);
      return (true); // Success
Failed:
     {
        if(NULL!=pRandStateObj) {ippsFree(pRandStateObj);}
        return (false);
   }
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
bool      vecf64_fill_ippsRandUniform_64f(Ipp64f * __restrict pDst, int32_t len,
                                          Ipp64f low, Ipp64f high) {

      IppsRandUniState_64f * __restrict pRandStateObj = NULL;
      IppStatus stat;
      int32_t sizeRndObj;
      uint32_t seed;
      stat = ippsRandUniformGetSize_64f(&sizeRndObj);
      if(stat != ippStsNoErr) { goto Failed;}
      pRandStateObj = (IppsRandUniState_64f*)ippsMalloc_64f(sizeRndObj);
      if(NULL==pRandStateObj && sizeRndObj != 0) { goto Failed;}
      _rdseed32_step(&seed);
      if(0==seed) {seed = (uint32_t)__rdtsc();}
      stat = ippsRandUniformInit_64f(pRandStateObj,low,high,seed);
      if(stat != ippStsNoErr) { goto Failed;}
      stat = ippsRandUniform_64f(pDst,len,pRandStateObj);
      if(stat != ippStsNoErr) { goto Failed;}
      ippsFree(pRandStatObj);
      return (true);
Failed:
     {
        if(NULL!=pRandStateObj) {ippsFree(pRandStateObj);}
        return (false);
   }      
}







__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
bool gms_vecf32_fill_ippsRandGauss_f32(Ipp32f * pDst, int32_len,
                                       Ipp32f mean,   Ipp32f stdDev) {

      IppsRandGaussState_32f * __restrict pRandGaussState = NULL;
      IppStatus stat;
      int32_t sizeRndObj;
      uint32_t seed;
      stat = ippsRandGaussGetSize_32f(&sizeRndObj);
      if(stat != ippStsNoErr) { goto Failed;}
      pRandGaussState = (IppsRandGaussState_32f*)ippMalloc_32f(sizeRndObj);
      if(NULL==pRandGaussState && 0 != sizeRndObj) { goto Failed;}
      _rdseed32_step(&seed);
      if(0==seed) { seed = (uint32_t)__rdtsc();}
      stat = ippsRandGaussInit_32f(pRandGaussState,mean,stdDev,seed);
      if(stat != ippStsNoErr) { goto Failed;}
      stat = ippsRandGauss_32f(pDst,len,pRandGaussState);
      if(stat != ippStsNoErr) { goto Failed;}
      ippsFree(pRandGaussState);
      return (true);
Failed:
         {
           if(NULL!=pRandGaussState) { ippsFree(pRandGaussState); }
	   return (false);
   }
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
bool gms_vecf64_fill_ippsRandGauss_f64(Ipp64f * pDst, int32_len,
                                       Ipp64f mean,   Ipp64f stdDev) {

      IppsRandGaussState_64f * __restrict pRandGaussState = NULL;
      IppStatus stat;
      int32_t sizeRndObj;
      uint32_t seed;
      stat = ippsRandGaussGetSize_64f(&sizeRndObj);
      if(stat != ippStsNoErr) { goto Failed;}
      pRandGaussState = (IppsRandGaussState_64f*)ippMalloc_64f(sizeRndObj);
      if(NULL==pRandGaussState && 0 != sizeRndObj) { goto Failed;}
      _rdseed32_step(&seed);
      if(0==seed) { seed = (uint32_t)__rdtsc();}
      stat = ippsRandGaussInit_64f(pRandGaussState,mean,stdDev,seed);
      if(stat != ippStsNoErr) { goto Failed;}
      stat = ippsRandGauss_64f(pDst,len,pRandGaussState);
      if(stat != ippStsNoErr) { goto Failed;}
      ippsFree(pRandGaussState);
      return (true);
Failed:
         {
           if(NULL!=pRandGaussState) { ippsFree(pRandGaussState); }
	   return (false);
   }
}







#endif /*__GMS_IPPS_COMPUTATIONAL_WRAPPERS_HPP__*/
