

#ifndef __GMS_IPPS_COMPUTATIONAL_WRAPPERS_HPP__
#define __GMS_IPPS_COMPUTATIONAL_WRAPPERS_HPP__

/*
    Wrappers for IPPS computational functions.
*/




#include <cstdint>
#include <ipps.h>
#include "GMS_config.h"

/*
        Vector Initialization Functions
*/


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__ ((assume_aligned(64)))
static inline
IppStatus gms_ippsCopy_32f(const Ipp32f * pSrc, Ipp32f * pDst, int32_t len) {

     return (ippsCopy_32f(pSrc,pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsCopy_64f(const Ipp64f * pSrc, Ipp64f * pDst, int32_t len) {

      return (ippsCopy_64f(pSrc,pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsCopy_32fc(const Ipp32fc * pSrc, Ipp32fc * pDst, int32_t len) {

      return (ippsCopy_32fc(pSrc,pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsCopy_64fc(const Ipp64fc * pSrc, Ipp64fc * pDst, int32_t len) {

      return (ippsCopy_64fc(pSrc,pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsCopy_8u(const Ipp8u * pSrc, Ipp8u * pDst, int32_t len) {

      return (ippsCopy_8u(pSrc,pDst,len));
}

/*
    Parameters
pSrc Pointer to the source vector.
pDst Pointer to the destination vector.
len  Number of elements to copy.
srcBitOffset Offset, in bits, from the first byte of the source vector.
dstBitOffset Offset, in bits, from the first byte of the destination vector.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsCopyLE_1u(const Ipp8u * pSrc, int32_t srcBitOffset,
                            Ipp8u * pDst, int32_t dstBitOffset, int32_t len) {

      return (ippsCopyLE_1u(pSrc,srcBitOffset,pDst,dstBitOffset,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsCopyBE_1u(const Ipp8u * pSrc, int32_t srcBitOffset,
                            Ipp8u * pDst, int32_t dstBitOffset, int32_t len) {

      return (ippsCopyBE_1u(pSrc,srcBitOffset,pDst,dstBitOffset,len));
}


/*
    Parameters
pSrc Pointer to the source vector used to initialize pDst .
pDst Pointer to the destination vector to be initialized.
len Number of elements to move.
    Description
This function moves the first len elements from a source vector pSrc into the destination vector pDst . If
some parts of the source and destination vectors are overlapping, then the function ensures that the original
source bytes in the overlapping parts are moved (it means that they are copied before being overwritten) to
the appropriate parts of the destination vector.
    Return Values
ippStsNoErr Indicates no error.
ippStsNullPtrErr Indicates an error when the pSrc or pDst pointer is NULL .
ippStsSizeErr Indicates an error when len is less than or equal to zero.


     Example
The example below shows how to use the function ippsMove .
Ipp8u pSrc[10] = { "123456789" };
Ipp8u pDst[6];
int len = 6;
IppStatus status;
status = ippsMove_8u ( pSrc, pDst, len );
if(ippStsNoErr != status)
printf("Intel(R) IPP Error: %s",ippGetStatusString(status));
Result:
pSrc = 123456789
pDst = 123456
*/

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsMove_8u(const Ipp8u * pSrc, Ipp8u * pDst, int32_t len) {

      return (ippsMove_8u(pSrc,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsMove_32f(const Ipp32f * pSrc, Ipp32f * pDst, int32_t len) {

       return (ippsMove_32f(pSrc,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsMove_64f(const Ipp64f * pSrc, Ipp64f * pDst, int32_t len) {

       return (ippsMove_64f(pSrc,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsMove_32fc(const Ipp32fc * pSrc, Ipp32fc * pDst, int32_t len) {

       return (ippsMove_32fc(pSrc,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsMove_64fc(const Ipp64fc * pSrc, Ipp64fc * pDst, int32_t len) {

       return (ippsMove_64fc(pSrc,pDst,len));
}

/*
   Set functions
   
   Initializes vector elements to a specified common
   value.
*/

/*
     Parameters
pDst Pointer to the vector to be initialized.
len Number of elements to initialize.
val Value used to initialize the vector pDst .
Description
This function initializes the first len elements of the real or complex vector pDst to contain the same value
val .
Return Values
ippStsNoErr Indicates no error.
ippStsNullPtrErr Indicates an error when the pDst pointer is NULL .
ippStsSizeErr Indicates an error when len is less than or equal to zero.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippSet_8u(Ipp8u val, Ipp8u * pDst, int32_t len) {

       return (ippsSet_8u(val,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippSet_32f(Ipp32f val, Ipp32f * pDst, int32_t len) {

       return (ippsSet_32f(val,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippSet_64f(Ipp64f val, Ipp64f * pDst, int32_t len) {

       return (ippsSet_64f(val,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippSet_32fc(Ipp32fc val, Ipp32fc * pDst, int32_t len) {

       return (ippsSet_32fc(val,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippSet_64fc(Ipp64fc val, Ipp64fc * pDst, int32_t len) {

       return (ippsSet_64fc(val,pDst,len));
}


/*
    Zero
    Initializes a vector to zero.
    
*/

/*
    Parameters
pDst Pointer to the vector to be initialized to zero.
len Number of elements to initialize.
Description
This function initializes the first len elements of the vector pDst to zero. If pDst is a complex vector, both
real and imaginary parts are zeroed.
Return Values
ippStsNoErr Indicates no error.
ippStsNullPtrErr Indicates an error when the pDst pointer is NULL .
ippStsSizeErr Indicates an error when len is less than or equal to zero.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsZero_8u(Ipp8u * pDst, int32_t len) {

        return (ippsZero_8u(pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsZero_32f(Ipp32f * pDst, int32_t len) {

        return (ippsZero_32f(pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsZero_64f(Ipp64f * pDst, int32_t len) {

        return (ippsZero_64f(pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsZero_32fc(Ipp32fc * pDst, int32_t len) {

        return (ippsZero_32fc(pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsZero_64fc(Ipp64fc * pDst, int32_t len) {

        return (ippsZero_64fc(pDst,len));
}

/*
    Sample-Generating Functions
*/

/*

    Parameters
pDst Magnitude of the tone, that is, the maximum value attained by
the wave.
Pointer to the phase of the tone relative to a cosine wave. It
must be in range [0.0, 2π). You can use the returned value to
compute the next continuous data block.
Frequency of the tone relative to the sampling frequency. It
must be in the interval [0.0, 0.5) for real tone and in [0.0, 1.0)
for complex tone.
Pointer to the array that stores the samples.
len Number of samples to be computed.


hint
Suggests using specific code. The possible values for the hint
argument are described in Hint Arguments.

     Description
This function generates the tone with the specified frequency rFreq , phase pPhase , and magnitude magn .
The function computes len samples of the tone, and stores them in the array pDst . For real tones, each
generated value x[n] is defined as:
x[n] = magn * cos(2πn*rFreq + phase)
For complex tones, x[n] is defined as:
x[n] = magn * (cos(2πn*rFreq + phase)+j* sin(2πn*rFreq + phase))
The parameter hint suggests using specific code, which provides for either fast but less accurate calculation,
or more accurate but slower execution.
Return Values
ippStsNoErr Indicates no error.
ippStsNullPtrErr Indicates an error when the pDst or pPhase pointer is NULL .
ippStsSizeErr Indicates an error when len is less than, or equal to zero.
ippStsToneMagnErr Indicates an error when magn is less than, or equal to zero.
ippStsToneFreqErr
ippStsTonePhaseErr
Indicates an error when rFreq is negative, or greater than, or
equal to 0.5 for real tone and to 1.0 for complex tone.
Indicates an error when the pPhase value is negative, or
greater than or equal to IPP_2PI .
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTone_32f(Ipp32f * pDst, int32_t len, Ipp32f mag, Ipp32f freq,
                           Ipp32f * pPhase, IppHintAlgorithm hint) {

	 return (ippsTone_32f(pDst,len,mag,freq,pPhase,hint));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTone_64f(Ipp64f * pDst, int32_t len, Ipp64f mag, Ipp64f freq,
                           Ipp64f * pPhase, IppHintAlgorithm hint) {

	 return (ippsTone_64f(pDst,len,mag,freq,pPhase,hint));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTone_32fc(Ipp32fc * pDst, int32_t len, Ipp32f mag, Ipp32f freq,
                           Ipp32f * pPhase, IppHintAlgorithm hint) {

	 return (ippsTone_32fc(pDst,len,mag,freq,pPhase,hint));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTone_64fc(Ipp64fc * pDst, int32_t len, Ipp64f mag, Ipp64f freq,
                           Ipp64f * pPhase, IppHintAlgorithm hint) {

	 return (ippsTone_64fc(pDst,len,mag,freq,pPhase,hint));
}


/*
    Triangle
    Generates a triangle with a given frequency, phase,
    and magnitude.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTriangle_32f(Ipp32f * pDst, int32_t len, Ipp32f mag, Ipp32f freq,
                               Ipp32f asym, Ipp32f * pPhase) {

	 return (ippsTriangle_32f(pDst,len,mag,freq,asym,pPhase));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTriangle_64f(Ipp64f * pDst, int32_t len, Ipp64f mag, Ipp64f freq,
                               Ipp64f asym, Ipp64f * pPhase) {

	 return (ippsTriangle_64f(pDst,len,mag,freq,asym,pPhase));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTriangle_32fc(Ipp32fc * pDst, int32_t len, Ipp32f mag, Ipp32f freq,
                               Ipp32f asym, Ipp32f * pPhase) {

	 return (ippsTriangle_32fc(pDst,len,mag,freq,asym,pPhase));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1, 2)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsTriangle_64fc(Ipp64fc * pDst, int32_t len, Ipp64f mag, Ipp64f freq,
                               Ipp64f asym, Ipp64f * pPhase) {

	 return (ippsTriangle_64fc(pDst,len,mag,freq,asym,pPhase));
}


/*
   RandUniformInit
   Initializes a noise generator with uniform distribution.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandUniformInit_32f(IppsRandUniState_32f * pRandUniState, Ipp32f low,
                                  Ipp32f high, uint32_t seed) {

	  return (ippsRandUniformInit_32f(&pRandUniState,low,high,seed));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandUniformInit_64f(IppsRandUniState_64f * pRandUniState, Ipp64f low,
                                  Ipp64f high, uint32_t seed) {

	  return (ippsRandUniformInit_64f(&pRandUniState,low,high,seed));
}


/*
    RandUniformGetSize
    Computes the length of the uniform distribution
    generator structure
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandUniformGetSize_32f(int * pRandUniformStateSize) {

          return (ippsRandUniformGetSize_32f(&pRandUniformStateSize));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandUniformGetSize_64f(int * pRandUniformStateSize) {

          return (ippsRandUniformGetSize_64f(&pRandUniformStateSize));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandUniform_32f(Ipp32f * pDst, int32_t len, IppsRandUniState_32f *
				  pRandUniState) {

	  return (ippsRandUniform_32f(pDst,len,pRandUniState));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandUniform_64f(Ipp64f * pDst, int32_t len, IppsRandUniState_64f *
				  pRandUniState) {

	  return (ippsRandUniform_64f(pDst,len,pRandUniState));
}

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

/*
        RandGaussInit
Initializes a noise generator with Gaussian
distribution.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandGaussInit_32f(IppsRandGaussState_32f * pRandGaussState, Ipp32f mean,
                                    Ipp32f stdDev, uint32_t seed) {

	return (ippsRandGaussinit_32f(pRandGaussState,mean,stdDev,seed));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandGaussInit_64f(IppsRandGaussState_64f * pRandGaussState, Ipp64f mean,
                                    Ipp64f stdDev, uint32_t seed) {

	return (ippsRandGaussinit_64f(pRandGaussState,mean,stdDev,seed));
}

/*
    RandGaussGetSize
Computes the length of the Gaussian distribution
generator structure.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandGaussGetSize_32f(int32_t * pRandGaussStateSize) {

        return (ippsRandGaussGetSize_32f(&pRandGaussStateSize));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1)))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandGaussGetSize_64f(int32_t * pRandGaussStateSize) {

        return (ippsRandGaussGetSize_64f(&pRandGaussStateSize));
}

/*
   RandGauss
Generates the pseudo-random samples with a
Gaussian distribution.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandGauss_32f(Ipp32f * pDst, int32_t len,
                                IppsRandGaussState_32f * pRandGaussState) {

	  return (ippsRandGauss_32f(pDst,len,pRandGaussState));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsRandGauss_64f(Ipp64f * pDst, int32_t len,
                                IppsRandGaussState_64f * pRandGaussState) {

	  return (ippsRandGauss_64f(pDst,len,pRandGaussState));
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

/*
    VectorJaehne
    Creates a Jaehne vector.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsVectorJaehne_32f(Ipp32f * pDst, int32_t len, Ipp32f magn) {

         return (ippsVectorJaehne_32f(pDst,len,magn));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsVectorJaehne_64f(Ipp64f * pDst, int32_t len, Ipp64f magn) {

         return (ippsVectorJaehne_64f(pDst,len,magn));
}

/*
    VectorSlope
Creates a slope vector.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsVectorSlope_32f(Ippf32 * pDst, int32_t len, Ipp32f offset, Ipp32f slope) {

         return (ippsVectorSlope_32f(pDst,len,offset,slope));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsVectorSlope_64f(Ippf64 * pDst, int32_t len, Ipp64f offset, Ipp64f slope) {

         return (ippsVectorSlope_64f(pDst,len,offset,slope));
}


/*
    AddC
    Adds a constant value to each element of a vector.
*/

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_32f(const Ipp32f * pSrc, Ipp32f val, Ipp32f * pDst, int32_t len) {

         return (ippsAddC_32f(pSrc,val,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_64f(const Ipp64f * pSrc, Ipp64f val, Ipp64f * pDst, int32_t len) {

         return (ippsAddC_64f(pSrc,val,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_32fc(const Ipp32fc * pSrc, Ipp32fc val, Ipp32fc * pDst, int32_t len) {

         return (ippsAddC_32fc(pSrc,val,pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_64fc(const Ipp64fc * pSrc, Ipp64fc val, Ipp64fc * pDst, int32_t len) {

         return (ippsAddC_64fc(pSrc,val,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_32f_I(Ipp32f val, Ipp32f * pSrcDst, int32_t len) {

         return (ippsAddC_32f_I(val,pSrcDst,len));
}



__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_64f_I(Ipp64f val, Ipp64f * pSrcDst, int32_t len) {

         return (ippsAddC_64f_I(val,pSrcDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_32fc_I(Ipp32fc val, Ipp32fc * pSrcDst, int32_t len) {

         return (ippsAddC_32fc_I(val,pSrcDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAddC_64fc_I(Ipp64fc val, Ipp64fc * pSrcDst, int32_t len) {

         return (ippsAddC_64fc_I(val,pSrcDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2,3)
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_32f(const Ipp32f * pSrc1, const Ipp32f * pSrc2, Ipp32f * pDst, int32_t len) {

         return (ippsAdd_32f(pSrc1,pSrc2,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2,3))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_64f(const Ipp64f * pSrc1, const Ipp64f * pSrc2, Ipp64f * pDst, int32_t len) {

         return (ippsAdd_64f(pSrc1,pSrc2,pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2,3))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_32fc(const Ipp32fc * pSrc1, const Ipp32fc * pSrc2, Ipp32fc * pDst, int32_t len) {

         return (ippsAdd_32fc(pSrc1,pSrc2,pDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2,3))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_64fc(const Ipp64fc * pSrc1, const Ipp64fc * pSrc2, Ipp64fc * pDst, int32_t len) {

         return (ippsAdd_64fc(pSrc1,pSrc2,pDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_32f_I(const Ipp32f * pSrc, Ipp32f * pSrcDst, int32_t len) {

          return (ippsAdd_32f_I(pSrc,pSrcDst,len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_64f_I(const Ipp64f * pSrc, Ipp64f * pSrcDst, int32_t len) {

          return (ippsAdd_64f_I(pSrc,pSrcDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_32fc_I(const Ipp32fc * pSrc, Ipp32fc * pSrcDst, int32_t len) {

          return (ippsAdd_32fc_I(pSrc,pSrcDst,len));
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__attribute__((nonnull (1,2))
__attribute__((assume_aligned(64)))
static inline
IppStatus gms_ippsAdd_64fc_I(const Ipp64fc * pSrc, Ipp64fc * pSrcDst, int32_t len) {

          return (ippsAdd_64fc_I(pSrc,pSrcDst,len));
}




#endif /*__GMS_IPPS_COMPUTATIONAL_WRAPPERS_HPP__*/
