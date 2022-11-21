

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


              //! Formula (1-39)
              // ! Dispersion coefficient
              typedef struct __ATTR_ALIGN__(64) f139_r4_t {

                       float                  * __restrict P;
                       float                  * __restrict sinth;
                       float                               ith;
                       float                               iph;
                       float                               omega;
                       int32_t                             nth;
                       int32_t                             nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,26)
#endif 
              } f139_r4_t;


              //! Formula (1-39)
              // ! Dispersion coefficient
              typedef struct __ATTR_ALIGN__(64) f139_r8_t {

                       double                  * __restrict P;
                       double                  * __restrict sinth;
                       double                               ith;
                       double                               iph;
                       double                               omega;
                       int32_t                              nth;
                       int32_t                              nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,16)
#endif 
              } f139_r8_t; 


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_c4_t {

                       struct JE_c4_t                        jec4;
                       struct eikr_c4_t                      ec4;
                       std::complex<float>                   ifac;
                       std::complex<float>      * __restrict he_x;
                       std::complex<float>      * __restrict he_y;
                       std::complex<float>      * __restrict he_z;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,60)
#endif
             } Hve_c4_t;


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_c8_t {

                       struct JE_c8_t                         jec8;
                       struct eikr_c8_t                       ec8;
                       std::complex<double>                   ifac;
                       std::complex<double>      * __restrict he_x;
                       std::complex<double>      * __restrict he_y;
                       std::complex<double>      * __restrict he_z;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,52)
#endif
             } Hve_c8_t;


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_r4_t {

                       struct JE_r4_t                        jer4;
                       struct eikr_r4_t                      er4;
                       float                                 if_re;
                       float                                 if_im;
                       float                    * __restrict he_xr;
                       float                    * __restrict he_xi;
                       float                    * __restrict he_yr;
                       float                    * __restrict he_yi;
                       float                    * __restrict he_zr;
                       float                    * __restrict he_zi;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Hve_r4_t;


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_r8_t {

                       struct JE_r8_t                         jer8;
                       struct eikr_r8_t                       er8;
                       double                                 if_re;
                       double                                 if_im;
                       double                    * __restrict he_xr;
                       double                    * __restrict he_xi;
                       double                    * __restrict he_yr;
                       double                    * __restrict he_yi;
                       double                    * __restrict he_zr;
                       double                    * __restrict he_zi;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,30)
#endif
             } Hve_r8_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hvm_c4_t {

                       struct JM_c4_t                        jmc4;
                       struct eikr_c4_t                      ec4;
                       std::complex<float>                   ifac;
                       std::complex<float>      * __restrict hm_x;
                       std::complex<float>      * __restrict hm_y;
                       std::complex<float>      * __restrict hm_z;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,60)
#endif
             } Hvm_c4_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hvm_c8_t {

                       struct JM_c8_t                         jmc8;
                       struct eikr_c8_t                       ec8;
                       std::complex<double>                   ifac;
                       std::complex<double>      * __restrict hm_x;
                       std::complex<double>      * __restrict hm_y;
                       std::complex<double>      * __restrict hm_z;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,52)
#endif
             } Hvm_c8_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hve_r4_t {

                       struct JM_r4_t                        jmr4;
                       struct eikr_r4_t                      er4;
                       float                                 if_re;
                       float                                 if_im;
                       float                    * __restrict hm_xr;
                       float                    * __restrict hm_xi;
                       float                    * __restrict hm_yr;
                       float                    * __restrict hm_yi;
                       float                    * __restrict hm_zr;
                       float                    * __restrict hm_zi;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Hvm_r4_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hvm_r8_t {

                       struct JM_r8_t                         jmr8;
                       struct eikr_r8_t                       er8;
                       double                                 if_re;
                       double                                 if_im;
                       double                    * __restrict hm_xr;
                       double                    * __restrict hm_xi;
                       double                    * __restrict hm_yr;
                       double                    * __restrict hm_yi;
                       double                    * __restrict hm_zr;
                       double                    * __restrict hm_zi;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,30)
#endif
             } Hvm_r8_t;


             //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_r4_t {

                       struct JE_r4_t                         jer4;
                       struct eikr_r4_t                       er4;
                       float                     * __restrict costh;
                       float                     * __restrict ne_xr;
                       float                     * __restrict ne_xi;
                       float                     * __restrict ne_yr;
                       float                     * __restrict ne_yi;
                       float                     * __restrict ne_zr;
                       float                     * __restrict ne_zi;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_r4_t;

              //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_r8_t {

                       struct JE_r8_t                          jer8;
                       struct eikr_r8_t                         er8;
                       double                     * __restrict costh;
                       double                     * __restrict ne_xr;
                       double                     * __restrict ne_xi;
                       double                     * __restrict ne_yr;
                       double                     * __restrict ne_yi;
                       double                     * __restrict ne_zr;
                       double                     * __restrict ne_zi;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_r8_t;


             //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_c4_t {

                       struct JE_c4_t                         jec4;
                       struct eikr_c4_t                       ec4;
                       float                     * __restrict costh;
                       std::complex<float>       * __restrict ne_x;
                       std::complex<float>       * __restrict ne_y;
                       std::complex<float>       * __restrict ne_z;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_c4_t;


              //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_c8_t {

                       struct JE_c8_t                          jcr8;
                       struct eikr_c8_t                         ec8;
                       double                     * __restrict costh;
                       std::complex<double>       * __restrict ne_x;
                       std::complex<double>       * __restrict ne_y;
                       std::complex<double>       * __restrict ne_z;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_c8_t;


             //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_r4_t {

                       struct JM_r4_t                         jmr4;
                       struct eikr_r4_t                       er4;
                       float                     * __restrict costh;
                       float                     * __restrict nm_xr;
                       float                     * __restrict nm_xi;
                       float                     * __restrict nm_yr;
                       float                     * __restrict nm_yi;
                       float                     * __restrict nm_zr;
                       float                     * __restrict nm_zi;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_r4_t;

              //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_r8_t {

                       struct JM_r8_t                          jmr8;
                       struct eikr_r8_t                         er8;
                       double                     * __restrict costh;
                       double                     * __restrict nm_xr;
                       double                     * __restrict nm_xi;
                       double                     * __restrict nm_yr;
                       double                     * __restrict nm_yi;
                       double                     * __restrict nm_zr;
                       double                     * __restrict nm_zi;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_r8_t;


             //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_c4_t {

                       struct JM_c4_t                         jmc4;
                       struct eikr_c4_t                       ec4;
                       float                     * __restrict costh;
                       std::complex<float>       * __restrict nm_x;
                       std::complex<float>       * __restrict nm_y;
                       std::complex<float>       * __restrict nm_z;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_c4_t;


              //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_c8_t {

                       struct JM_c8_t                          jmc8;
                       struct eikr_c8_t                         ec8;
                       double                     * __restrict costh;
                       std::complex<double>       * __restrict nm_x;
                       std::complex<double>       * __restrict nm_y;
                       std::complex<double>       * __restrict nm_z;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_c8_t;



            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_r4_t {

                       struct Nev_r4_t                         ner4;
                       struct eikr_r4_t                        er4;   // 96-bytes
                       float                      * __restrict hv_xr;
                       float                      * __restrict hv_xi;
                       float                      * __restrict hv_yr;
                       float                      * __restrict hv_yi;
                       float                      * __restrict hv_zr;
                       float                      * __restrict hv_zi;
                       float                                   if_re;
                       float                                   if_im; // 54-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,42)
#endif           
            } FFhev_r4_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_r8_t {

                       struct Nev_r8_t                          ner8;
                       struct eikr_r8_t                         er8;   // 96-bytes
                       double                      * __restrict hv_xr;
                       double                      * __restrict hv_xi;
                       double                      * __restrict hv_yr;
                       double                      * __restrict hv_yi;
                       double                      * __restrict hv_zr;
                       double                      * __restrict hv_zi;
                       double                                   if_re;
                       double                                   if_im; // 54-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,34)
#endif           
            } FFhev_r8_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_c4_t {

                       struct Nev_c4_t                         nec4;
                       struct eikr_c4_t                        ec4;   // 96-bytes
                       std::complex<float>        * __restrict hv_x;
                       std::complex<float>        * __restrict hv_y;
                       std::complex<float>        * __restrict hv_z;
                       std::complex<float>                     ifac; // 32-bytes
         
            } FFhev_r4_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_c8_t {

                       struct Nev_c8_t                          nec8;
                       struct eikr_c8_t                         ec8;   // 96-bytes
                       std::complex<double>        * __restrict hv_x;
                       std::complex<double>        * __restrict hv_y;
                       std::complex<double>        * __restrict hv_z;
                       std::complex<double>                     ifac;
                       
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,56)
#endif           
            } FFhev_c8_t;

               
       } // radiolocation

} // gms














#endif /*__GMS_ANTENNA_TYPES_H__*/
