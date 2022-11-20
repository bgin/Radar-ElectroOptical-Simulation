

#ifndef __GMS_ANTENNA_TYPES_H__
#define __GMS_ANTENNA_TYPES_H__ 201120221117


namespace file_info {

     const unsigned int GMS_ANTENNA_TYPES_MAJOR = 1;
     const unsigned int GMS_ANTENNA_TYPES_MINOR = 1;
     const unsigned int GMS_ANTENNA_TYPES_MICRO = 0;
     const unsigned int GMS_ANTENNA_TYPES_FULLVER =
       1000U*GMS_ANTENNA_TYPES_MAJOR+100U*GMS_ANTENNA_TYPES_MINOR+
       10U*GMS_ANTENNA_TYPES_MICRO;
     const char * const GMS_ANTENNA_TYPES_CREATION_DATE = "20-11-2022 11:17 +00200 (SUN 20 NOV 2022 GMT+2)";
     const char * const GMS_ANTENNA_TYPES_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_ANTENNA_TYPES_SYNOPSIS      = "Antenna common data types."

}


/*
 Purpose:
 !                        Derived data types for 'antenna_sensor' module implementation.
 !                        Various characteristics of different antenna types  
 !                        Based mainly on book titled (rus):          
 !                        Проектирование антенно фидерных устройств. Жук М.С. Молочков Ю.Б
*/


#include <cstdint>
#include <complex>
#include "GMS_config"


namespace gms {


          namespace  radiolocation

              typedef struct __ATTR_ALIGN__(32) E_c4_t {
                      // Complex electric field.
                      std::complex<float> * __restrict e_x
                      std::complex<float> * __restrict e_y
                      std::complex<float> * __restrict e_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif
              } E_c4_t;

              
               typedef struct __ATTR_ALIGN__(32) H_c4_t {
                      // Complex magnetic field.
                      std::complex<float> * __restrict h_x
                      std::complex<float> * __restrict h_y
                      std::complex<float> * __restrict h_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif
              } H_c4_t;


               typedef struct __ATTR_ALIGN__(32) E_c8_t {
                      // Complex electric field.
                      std::complex<double> * __restrict e_x
                      std::complex<double> * __restrict e_y
                      std::complex<double> * __restrict e_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif
              } E_c8_t;

              
               typedef struct __ATTR_ALIGN__(32) H_c8_t {
                      // Complex magnetic field.
                      std::complex<double> * __restrict h_x
                      std::complex<double> * __restrict h_y
                      std::complex<double> * __restrict h_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif
              } H_c8_t;
   

              typedef struct __ATTR_ALIGN__(64) E_r4_t {
                      //  ! Complex Electric  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                      float * __restrict e_xr;
                      float * __restrict e_xi;
                      float * __restrict e_yr;
                      float * __restrict e_yi;
                      float * __restrict e_zr;
                      float * __restrict e_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif
              } E_r4_t;


              typedef struct __ATTR_ALIGN__(64) H_r4_t {
                      //  ! Complex Magnetic  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                      float * __restrict h_xr;
                      float * __restrict h_xi;
                      float * __restrict h_yr;
                      float * __restrict h_yi;
                      float * __restrict h_zr;
                      float * __restrict h_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif
              } H_r4_t;


              typedef struct __ATTR_ALIGN__(64) E_r8_t {
                      //  ! Complex Electric  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                      double * __restrict e_xr;
                      double * __restrict e_xi;
                      double * __restrict e_yr;
                      double * __restrict e_yi;
                      double * __restrict e_zr;
                      double * __restrict e_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif
              } E_r8_t;


              typedef struct __ATTR_ALIGN__(64) H_r8_t {
                      //  ! Complex Magnetic  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                      double * __restrict h_xr;
                      double * __restrict h_xi;
                      double * __restrict h_yr;
                      double * __restrict h_yi;
                      double * __restrict h_zr;
                      double * __restrict h_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif
              } H_r8_t;
               
       }

}














#endif /*__GMS_ANTENNA_TYPES_H__*/
