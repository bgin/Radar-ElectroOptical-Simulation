
#ifndef __GMS_HOST_DEV_COMPLEX_H__
#define __GMS_HOST_DEV_COMPLEX_H__


/*

     Simple host device complex arithmetic auxiliary routines.
*/

#include <cuComplex.h>



__host__ __device__
__forceinline__ cuFloatComplex
create_complex_c4_1(const float re,
		  const float im) {
  return (make_cuFloatComplex(re,im));
}


__host__ __device__
__forceinline__ cuFloatDouble
create_complex_c8_1(const double re,
		  const double im) {
  return (make_cuDoubleComplex(re,im));
}


__host__ __device__
__forceinline__ float
real_c4_1(const cuFloatComplex c) {
    return (cuCrealf(c));
}


__host__ __device__
__forceinline__ double
real_c8_1(const cuDoubleComplex c) {
    return (cuCreal(c));
}


__host__ __device__
__forceinline__ float
imag_c4_1(const cuFloatComplex c) {
     return (cuCimagf(c));
}


__host__ __device__
__forceinline__ double
imag_c8_1(const cuDoubleComplex c) {
     return (cuCimag(c));
}


__host__ __device__
__forceinline__ cuFloatComplex
conjugate_c4_1(const cuFloatComplex c) {
      return (cuConjf(c));
}


__host__ __device__
__forceinline__ cuDoubleComplex
conjugate_c8_1(const cuDoubleComplex c) {
      return (cuConj(c));
}


__host__ __device__
__forceinline__ float
cabs_c4_1(const cuFloatComplex c) {
       return (cuCabsf(c));
}


__host__ __device__
__forceinline__ double
cabs_c8_1(const cuDoubleComplex c) {
       return (cuCabs(c));
}


__host__ __device__
__forceinline__ cuFloatComplex
cdiv_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
       return (cuCdivf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
cdiv_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
       return (cuCdiv(c1,c2));
}


__host__ __device__
__forceinline__ cuFloatComplex
cmul_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
       return (cuCmulf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
cmul_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
       return (cuCmul(c1,c2));
}


__host__ __device__
__forceinline__ cuFloatComplex
csub_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
      return (cuCsubf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
csub_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
      return (cuCsub(c1,c2));
}


__host__ __device__
__forceinline__ cuFloatComplex
cadd_c4_1(const cuFloatComplex c1,
          const cuFloatComplex c2) {
      return (cuCaddf(c1,c2));
}


__host__ __device__
__forceinline__ cuDoubleComplex
cadd_c8_1(const cuDoubleComplex c1,
          const cuDoubleComplex c2) {
      return (cuCadd(c1,c2));
}


__host__ __device__
__forceinline__ float
cmag_c4_1(const cuFloatComplex c) {
    const float re=cuCrealf(c);
    const float im=cuCimagf(c);
    return (re*re+im*im);
}


__host__ __device__
__forceinline__ double
cmag_c8_1(const cuDoubleComplex c) {
    const double re=cuCreal(c);
    const double im=cuCimag(c);
    return (re*re+im*im);
}















#endif /*__GMS_HOST_DEV_COMPLEX_H__*/
