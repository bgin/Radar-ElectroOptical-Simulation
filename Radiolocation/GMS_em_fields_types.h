

#ifndef __GMS_EM_FIELDS_TYPES_H__
#define __GMS_EM_FIELDS_TYPES_H__ 190120221323

namespace file_info {

     const unsigned int GMS_EM_FIELDS_TYPES_MAJOR = 1;
     const unsigned int GMS_EM_FIELDS_TYPES_MINOR = 1;
     const unsigned int GMS_EM_FIELDS_TYPES_MICRO = 0;
     const unsigned int GMS_EM_FIELDS_TYPES_FULLVER =
       1000U*GMS_EM_FIELDS_TYPES_MAJOR+100U*GMS_EM_FIELDS_TYPES_MINOR+
       10U*GMS_EM_FIELDS_TYPES_MICRO;
     const char * const GMS_EM_FIELDS_TYPES_CREATION_DATE = "19-01-2022 13:23 +00200 (WED 19 JAN 2022 GMT+2)";
     const char * const GMS_EM_FIELDS_TYPES_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_EM_FIELDS_TYPES_SYNOPSIS      = "Electric and Magnetic Fields PAoS and SoA data types"

}

#include <cstdint>
#include "GMS_avx512c16f32.h"
#include "GMS_avx512c8f64.h"
#include "GMS_avx512vecf32.h"
#include "GMS_avx512vecf64.h"
#include "GMS_config.h"


namespace gms {


                /*
                           !==============================================================
                           !                   PAoS data types
                           !==============================================================
                 */

                 //  Complex Electric Field 1D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(16) H_X_C16 {

		       ZMM16c4 * __restrict H_x;
                       int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif
		 } H_X_C16;


		 // Complex Electric Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) H_XY_C16 {

                       ZMM16c4 * __restrict H_x;
		       ZMM16c4 * __restrict H_y;
		       int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       
		 } H_XY_C16;


		 //  Complex Electric Field 3D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) H_XYZ_C16 {

                       ZMM16c4 * __restrict H_x;
		       ZMM16c4 * __restrict H_y;
		       ZMM16c4 * __restrict H_z;
		       int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif			       
		 } H_XYZ_C16;

using ZMM16r4 = AVX512Vec16;


                 // ! Real Electric Field 3D (Packed-AoS) data type decomposed.
                 typedef struct __ATTR_ALIGN__(16) H_X_R16 {

                       ZMM16r4 * __restrict H_x;
		       int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif		       
		 } H_X_R16;


		 // Real Electric Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) H_XY_R16 {

		       ZMM16r4 * __restrict H_x;
		       ZMM16r4 * __restrict H_y;
		       int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       
                 } H_XY_R16;


		 // Real Electric Field 3D (Packed-AoS) data type decomposed.
                 typedef struct __ATTR_ALIGN__(32) H_XYZ_R16 {

                        ZMM16r4 * __restrict H_x;
			ZMM16r4 * __restrict H_y;
			ZMM16r4 * __restrict H_z;
                        int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif					
		 } H_XYZ_R16;


		 // Complex Magnetic Field 1D (Packed-AoS) data type decomposed.
                 typedef struct __ATTR_ALIGN__(16) B_X_C16 {

                        ZMM16c4 * __restrict B_x;
                        int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif		       		     
		 } B_X_C16;


		 // Complex Magnetic  Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XY_C16 {

                         ZMM16c4 * __restrict B_x;
			 ZMM16c4 * __restrict B_y;
                         int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       			 
		 } B_XY_C16;


		 //  Complex Magnetic Field 3D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XYZ_C16 {

                        ZMM16c4 * __restrict B_x;
			ZMM16c4 * __restrict B_y;
			ZMM16c4 * __restrict B_z;
                        int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif						
		 } B_XYZ_C16;


		 //! Real Magnetic Field 1D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(16) B_X_R16 {

                        ZMM16r4 * __restrict B_x;
                        int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif		       				
		 } B_X_R16;


		 // Real Magnetic Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XY_R16 {

                         ZMM16r4 * __restrict B_x;
			 ZMM16r4 * __restrict B_y;
                         int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       				 
		 } B_XY_R16;


		 // Real Magnetic Field 3D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XYZ_R16 {

                         ZMM16r4 * __restrict B_x;
			 ZMM16r4 * __restrict B_y;
			 ZMM16r4 * __restrict B_z;
                         int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif						 
		 } B_XYZ_R16;


		 //  Complex Surface currents 3D (Packed-AoS) data type.
		 typedef struct __ATTR_ALIGN__(32) Js_XYZ_C16 {

                        ZMM16c4 * __restrict Js_x;
			ZMM16c4 * __restrict Js_y;
			ZMM16c4 * __restrict Js_z;
                        int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif				
		 } Js_XYZ_C16;


		 /*
                          !=======================================================================================!
                          !!               Packed AoS double precision data types definitions
                          !=======================================================================================!
                    */


		 // Complex Electric Field 1D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(16) H_X_C8 {

                         ZMM8c8 * __restrict H_x;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif		       	
		 } H_X_C8;


		 // Complex Electric Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) H_XY_C8 {

                         ZMM8c8 * __restrict H_x;
			 ZMM8c8 * __restrict H_y;
                         int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       				 
		 } H_XY_C8;


		 // Complex Electric Field 3D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) H_XYZ_C8 {

                         ZMM8c8 * __restrict H_x;
			 ZMM8c8 * __restrict H_y;
			 ZMM8c8 * __restrict H_z;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif			
		 } H_XYZ_C8;

using ZMM8r8 = AVX512Vec8;

		 // Real Electric Field 1D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(16) H_X_R8 {

                         ZMM8r8 * __restrict H_x;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif		       	
		 } H_X_R8;

		 //  Real Electric Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) H_XY_R8 {

                         ZMM8r8 * __restrict H_x;
			 ZMM8r8 * __restrict H_y;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       	
		 } H_XY_R8;


		 // Real Electric Field 3D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) H_XYZ_R8 {

                         ZMM8r8 * __restrict H_x;
			 ZMM8r8 * __restrict H_y;
			 ZMM8r8 * __restrict H_z;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif			
		 } H_XYZ_R8;


		 // Complex Magnetic Field 1D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(16) B_X_C8 {

                         ZMM8c8 * __restrict B_x;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif		       	
		 } B_X_C8;


		 // Complex Magnetic  Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XY_C8 {

                         ZMM8c8 * __restrict B_x;
			 ZMM8c8 * __restrict B_y;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       	
		 } B_XY_C8;


		 // Complex Magnetic Field 3D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XYZ_C8 {

                         ZMM8c8 * __restrict B_x;
			 ZMM8c8 * __restrict B_y;
			 ZMM8c8 * __restrict B_z;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif			
		 } B_XYZ_C8;


		 // Real Magnetic Field 1D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(16) B_X_R8 {

                         ZMM8r8 * __restrict B_x;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif		       	
		 } B_X_R8;


		 // Real Magnetic  Field 2D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XY_R8 {

                         ZMM8r8 * __restrict B_x;
			 ZMM8r8 * __restrict B_y;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,12)
#endif		       	
		 } B_XY_R8;


		 // Real Magnetic Field 3D (Packed-AoS) data type decomposed.
		 typedef struct __ATTR_ALIGN__(32) B_XYZ_R8 {

                         ZMM8r8 * __restrict B_x;
			 ZMM8r8 * __restrict B_y;
			 ZMM8r8 * __restrict B_z;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif			
		 } B_XYZ_R8;


		 // Complex Surface currents 3D (Packed-AoS) data type.
		 typedef struct __ATTR_ALIGN__(32) Js_XYZ_C8 {

                         ZMM8c8 * __restrict Js_x;
			 ZMM8c8 * __restrict Js_y;
			 ZMM8c8 * __restrict Js_z;
			 int32_t n_pts; // n-points of evaluation, i.e. R-points
#if(USE_STRUCT_PADDING) == 1		       
		       PAD_TO(1,4)
#endif			
		 } Js_XYZ_C8;


		 
		 

} // GMS


#endif /*__GMS_EM_FIELDS_TYPES_H__*/
