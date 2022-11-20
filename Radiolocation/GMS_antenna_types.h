

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


              typedef struct __ATTR_ALIGN__(32) JE_c4_t {
                      // Complex magnetic current
                      std::complex<float> * __restrict je_x;
                      std::complex<float> * __restrict je_y;
                      std::complex<float> * __restrict je_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JE_c4_t;


              typedef struct __ATTR_ALIGN__(32) JM_c4_t {
                      // Complex electric current
                      std::complex<float> * __restrict jm_x;
                      std::complex<float> * __restrict jm_y;
                      std::complex<float> * __restrict jm_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JM_c4_t;


              typedef struct __ATTR_ALIGN__(32) JE_c8_t {
                      // Complex magnetic current
                      std::complex<double> * __restrict je_x;
                      std::complex<double> * __restrict je_y;
                      std::complex<double> * __restrict je_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JE_c8_t;


              typedef struct __ATTR_ALIGN__(32) JM_c8_t {
                      // Complex electric current
                      std::complex<double> * __restrict jm_x;
                      std::complex<double> * __restrict jm_y;
                      std::complex<double> * __restrict jm_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JM_c8_t;


              typedef struct __ATTR_ALIGN__(64) JE_r4_t {
                     // ! Complex Electric Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      float * __restrict je_xr;
                      float * __restrict je_xi;
                      float * __restrict je_yr;
                      float * __restrict je_yi;
                      float * __restrict je_zr;
                      float * __restrict je_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JE_r4_t;


              typedef struct __ATTR_ALIGN__(64) JM_r4_t {
                     // ! Complex Magnetic Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      float * __restrict jm_xr;
                      float * __restrict jm_xi;
                      float * __restrict jm_yr;
                      float * __restrict jm_yi;
                      float * __restrict jm_zr;
                      float * __restrict jm_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JM_r4_t;


              typedef struct __ATTR_ALIGN__(64) JE_r8_t {
                     // ! Complex Electric Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      double * __restrict je_xr;
                      double * __restrict je_xi;
                      double * __restrict je_yr;
                      double * __restrict je_yi;
                      double * __restrict je_zr;
                      double * __restrict je_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JE_r8_t;


              typedef struct __ATTR_ALIGN__(64) JM_r8_t {
                     // ! Complex Magnetic Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      double * __restrict jm_xr;
                      double * __restrict jm_xi;
                      double * __restrict jm_yr;
                      double * __restrict jm_yi;
                      double * __restrict jm_zr;
                      double * __restrict jm_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JM_r8_t;


              typedef struct __ATTR_ALIGN__(32) eikr_c4_t {
                      // Time-Harmonic complex exponential
                      float               * __restrict R;
                      std::complex<float> * __restrict ce;
                      float                            k;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif              
              } eikr_c4_t;


               typedef struct __ATTR_ALIGN__(32) eikr_c8_t {
                      // Time-Harmonic complex exponential
                      double               * __restrict R;
                      std::complex<double> * __restrict ce;
                      double                            k;
                      int32_t                           npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif              
              } eikr_c8_t;


              typedef struct __ATTR_ALIGN__(32) eikr_r4_t {
                      // ! Time-Harmonic complex exponential decomposed into 
                      // ! real and imaginary parts
                      float                 * __restrict R;
                      float                 * __restrict e_re;
                      float                 * __restrict e_im;
                      float                              k;
                      int32_t                            npts;
              } eikr_r4_t;


              typedef struct __ATTR_ALIGN__(32) eikr_r8_t {
                      // ! Time-Harmonic complex exponential decomposed into 
                      // ! real and imaginary parts
                      double                 * __restrict R;
                      double                 * __restrict e_re;
                      double                 * __restrict e_im;
                      double                              k;
                      int32_t                             npts;
              } eikr_r8_t;


              // ! Formula (1-37)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f137_r4_t {
                       
                       float                 * __restrict sinth;
                       float                 * __restrict F;
                       float                              ith;
                       float                              iph;
                       float                              ifac;
                       float                              omega;
                       float                              avsl;
                       int32_t                            nth;
                       int32_t                            nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif  
              } f137_r4_t;


              // ! Formula (1-37)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f137_r8_t {
                       
                       double                 * __restrict sinth;
                       double                 * __restrict F;
                       double                              ith;
                       double                              iph;
                       double                              ifac;
                       double                              omega;
                       double                              avsl;
                       int32_t                             nth;
                       int32_t                             nph;
              } f137_r8_t;


              // ! Formula (1-38)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f138_r4_t {
                       
                       float                 * __restrict sinth;
                       float                 * __restrict Fsqr;
                       float                              ith;
                       float                              iph;
                       float                              ifac;
                       float                              omega;
                       float                              avsl;
                       int32_t                            nth;
                       int32_t                            nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif  
              } f138_r4_t;


              // ! Formula (1-38)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f138_r8_t {
                       
                       double                 * __restrict sinth;
                       double                 * __restrict Fsqr;
                       double                              ith;
                       double                              iph;
                       double                              ifac;
                       double                              omega;
                       double                              avsl;
                       int32_t                             nth;
                       int32_t                             nph;
              } f138_r8_t;

               
       }

}














#endif /*__GMS_ANTENNA_TYPES_H__*/
