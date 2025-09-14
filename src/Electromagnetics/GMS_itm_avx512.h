
#ifndef __GMS_ITM_AVX512_H__
#define __GMS_ITM_AVX512_H__

/*
    

SOFTWARE DISCLAIMER / RELEASE

This software was developed by employees of the National Telecommunications and Information Administration (NTIA), an agency of the Federal Government and is provided to you as a public service. Pursuant to Title 15 United States Code Section 105, works of NTIA employees are not subject to copyright protection within the United States.

The software is provided by NTIA “AS IS.” NTIA MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NTIA does not warrant or make any representations regarding the use of the software or the results thereof, including but not limited to the correctness, accuracy, reliability or usefulness of the software.

To the extent that NTIA holds rights in countries other than the United States, you are hereby granted the non-exclusive irrevocable and unconditional right to print, publish, prepare derivative works and distribute the NTIA software, in any medium, or authorize others to do so on your behalf, on a royalty-free basis throughout the World.

You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change.

You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property.

Please provide appropriate acknowledgments of NTIA’s creation of the software in any copies or derivative works of this software.
   
   
   Manually vectorized by Bernard Gingold on 11-12-2021 18:59
*/
namespace file_info {

     const unsigned int GMS_ITM_AVX512_MAJOR = 1;
     const unsigned int GMS_ITM_AVX512_MINOR = 0;
     const unsigned int GMS_ITM_AVX512_MICRO = 0;
     const unsigned int GMS_ITM_AVX512_FULLVER = 1000*GMS_ITM_AVX512_MAJOR+
                                                 100*GMS_ITM_AVX512_MINOR+
						 10*GMS_ITM_AVX512_MICRO;
     const char * const GMS_ITM_AVX512_CREATION_DATE = "11-12-2021 10:45 AM +00200 (SAT 11 DEC 2021 GMT+2)";
     const char * const GMS_ITM_AVX512_BUILD_DATE    = __DATE__ ":"__TIME__;
     const char * const GMS_ITM_AVX512_PROGRAMMER    = "Programmer: Bernard Gingold, contact: beniekg@gmail.com, adapted from ITM v1.4 model.";
     const char * const GMS_ITM_AVX512_DESCRIPTION   = "Manual AVX512 vectorization of ITM v1.4 suitable functions."
}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


// List of valid polarizations
#define POLARIZATION__HORIZONTAL                0
#define POLARIZATION__VERTICAL                  1

// List of valid siting criteria
#define SITING_CRITERIA__RANDOM                 0
#define SITING_CRITERIA__CAREFUL                1
#define SITING_CRITERIA__VERY_CAREFUL           2

// List of valid radio climates
#define CLIMATE__EQUATORIAL                     1
#define CLIMATE__CONTINENTAL_SUBTROPICAL        2
#define CLIMATE__MARITIME_SUBTROPICAL           3
#define CLIMATE__DESERT                         4
#define CLIMATE__CONTINENTAL_TEMPERATE          5
#define CLIMATE__MARITIME_TEMPERATE_OVER_LAND   6
#define CLIMATE__MARITIME_TEMPERATE_OVER_SEA    7

// List of valid modes of propagation
#define MODE__NOT_SET                           0
#define MODE__LINE_OF_SIGHT                     1
#define MODE__DIFFRACTION                       2
#define MODE__TROPOSCATTER                      3

namespace gms {

          namespace math {


	        typedef struct __attribute__((aligned(64))) CDatav16 {

		        __m512 theta_hzn;     // Terminal horizon angles
			__m512 d_hzn_meter;   // Terminal horizon distances, in meters
                        __m512 h_e_meter;     // Terminal effective heights, in meters
			__m512 N_s;           // Surface refractivity, in N-Units
			__m512 delta_h_meter; // Terrain irregularity parameter, in meters
			__m512 A_ref_db;      // Reference attenuation, in dB
			__m512 A_fs_db;       // Free space basic transmission loss, in dB
			__m512 d_km;          // Path distance, in km
		}CDatav16;

		typedef struct __attribute__((aligned(64))) CDatav8 {

		        __m512d theta_hzn;     // Terminal horizon angles
			__m512d d_hzn_meter;   // Terminal horizon distances, in meters
                        __m512d h_e_meter;     // Terminal effective heights, in meters
			__m512 dN_s;           // Surface refractivity, in N-Units
			__m512d delta_h_meter; // Terrain irregularity parameter, in meters
			__m512d A_ref_db;      // Reference attenuation, in dB
			__m512d A_fs_db;       // Free space basic transmission loss, in dB
			__m512d d_km;          // Path distance, in km
		}CDatav8;


		

		    
                      __ATTR_ALWAYS_INLINE__
		     static inline
                      void
		      data_preload_v16() {

                          const volatile __m512 t0 = cl1[0];
			  const volatile __m512 t1 = cl1[1];
			  const volatile __m512 t2 = cl2[0];
			  const volatile __m512 t3 = cl2[1];
			  const volatile __m512 t4 = cl2[2];
			  const volatile __m512 t5 = cl2[3];
			  const volatile __m512 t6 = cl2[4];
			  const volatile __m512 t7 = cl2[5];
			  const volatile __m512 t8 = cl3[0];
			  const volatile __m512 t9 = cl3[1];
		      }
 /*=============================================================================
 |
 |  Description:  Compute delta_h_d
 |
 |        Input:  d__meter       - Path distance, in meters
 |                delta_h__meter - Terrain irregularity parameter
 |
 |      Outputs:  [None]
 |
 |      Returns:  delta_h_d      - Terrain irregularity of path
 |
 *===========================================================================*/


		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      __m512
		      itm_terrain_roughness_zmm16r4(const __m512 d_meter, // this variable shall be a negative
		                                    const __m512 delta_h_meter) {
                            const register __m512 _1 = _mm512_set1_ps(1.0F);
                            // // [ERL 79-ITS 67, Eqn 3], with distance in meters instead of kilometers
			    const register __m512 t0 = _mm512_div_ps(d_meter,cl1[0]);
			    const register __m512 t1 = _mm512_sub_ps(_1,_mm512_mul_ps(,_mm512_exp_ps(t0)));
			    return (_mm512_mul_ps(delta_h_meter,t1));
		      }

/*=============================================================================
 |
 |  Description:  Approximate to ideal knife edge diffraction loss
 |
 |        Input:  v2             - v^2 parameter
 |
 |      Outputs:  [None]
 |
 |      Returns:  A(v, 0)        - Loss, in dB
 |
 *===========================================================================*/

                    
                     __ATTR_ALWAYS_INLINE__
		     static inline
		     __m512
		     itm_fresnel_integral_zmm16r4(const __m512 v2) {

                          
			  const register __m512 t0    = _mm512_mul_ps(cl2[1],v2);
			  const register __m512 c0    = _mm512_sub_ps(_mm512_fmadd_ps(cl2[2],
			                                              _mm512_sqrt_ps(v2),
								      cl2[3],t0));
			  const register __m512 c1    = _mm512_fmadd_ps(cl2[4],_mm512_log10_ps(v2),cl2[5]);
			  const __m512 fi;
			  const __mmask16 k = 0x0;
			  k  = _mm512_cmp_mask_ps(v2,cl2[0],_CMP_LT_OQ);
			  fi = _mm512_mask_blend_ps(k,c1,c0);
			  return (fi);
		    }

/*=============================================================================
 |
 |  Description:  Compute the knife-edge diffraction loss
 |
 |        Input:  d__meter          - Distance of interest, in meters
 |                f__mhz            - Frequency, in MHz
 |                a_e__meter        - Effective earth radius, in meters
 |                theta_los         - Angular distance of line-of-sight region
 |                d_hzn__meter[2]   - Horizon distances, in meters
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_k__db        - Knife-edge diffraction loss, in dB
 |
 *===========================================================================*/		      

		      
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      __m512
		      itm_knife_edge_diffraction_zmm16r4(const __m512 d_meter,
		                                         const __m512 f_mhz,
						         const __m512 a_e_meter,
						         const __m512 theta_los,
						         const __m512 d_hzn_meter1,
						         const __m512 d_hzn_meter2) {
			    //const register __m512 _0_0795775   = _mm512_set1_ps(0.0795775F);
			    //const register __m512 _47_7        = _mm512_set1_ps(47.7F);
			    const register __m512 t0           = _mm512_mul_ps(cl3[0],_mm512_div_ps(f_mhz,cl3[1]));
			    const register __m512 d_ml_meter   = _mm512_sub_ps(d_hzn_meter1,d_hzn_meter2);
			    const register __m512 d_nlos_meter = _mm512_sub_ps(d_meter,d_ml_meter);
			    const register __m512 theta_nlos   = _mm512_sub_ps(_mm512_div_ps(d_meter,a_e_meter),theta_los);
			    const register __m512 t1           = _mm512_mul_ps(theta_nlos,theta_nlos);
			    const register __m512 t2           = _mm512_div_ps(d_nlos_meter,_mm512_add_ps(d_nlos_meter,d_hzn_meter1));
			    const register __m512 t3           = _mm512_div_ps(d_nlos_meter,_mm512_add_ps(d_nlos_meter,d_hzn_meter2));
			    const register __m512 v_1          = _mm512_mul_ps(t0,_mm512_mul_ps(_mm512_mul_ps(t1,d_hzn_meter1),t2));
			    const register __m512 v_2          = _mm512_mul_ps(t0,_mm512_mul_ps(_mm512_mul_ps(t1,d_hzn_meter2),t3));
			    const register __m512 A_k_db       = itm_fresnel_integral_zmm16r4(v_1);
			                                         itm_fresnel_integral_zmm16r4(v_2);
			    return (A_k_db);
		      }

/*=============================================================================
 |
 |  Description:  Compute the smooth earth diffraction loss using the 
 |                Vogler 3-radii method
 |
 |        Input:  d__meter          - Path distance, in meters
 |                f__mhz            - Frequency, in MHz
 |                a_e__meter        - Effective earth radius, in meters
 |                theta_los         - Angular distance of line-of-sight region
 |                d_hzn__meter[2]   - Horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                Z_g               - Complex ground impedance
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_r__db           - Smooth-earth diffraction loss, in dB
 |
 *===========================================================================*/
		      

		     
                      __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      __m512
		      itm_smooth_earth_diffraction_zmm16r4(const __m512 d_meter,
		                                           const __m512 f_mhz,
							   const __m512 a_e_meter,
							   const __m512 theta_los,
							   const __m512 d_hzn_meter1,
							   const __m512 d_hzn_meter2,
							   const __m512 d_h_e_meter1,
							   const __m512 d_h_e_meter2,
							   const __m512 Z_g_re,
							   const __m512 Z_g_im); 

/*=============================================================================
 |
 |  Description:  Height Function, F(x, K) for smooth earth diffraction
 |
 |        Input:  x__km          - Normalized distance, in meters
 |                K              - K value
 |
 |      Outputs:  [None]
 |
 |      Returns:  F(x, K)        - in dB
 |
 *===========================================================================*/		      

                      
                      __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      __m512
		      itm_height_function_zmm16r4(const __m512 x_km,
		                                  const __m512 K); 

/*=============================================================================
 |
 |  Description:  Initialize parameters for point-to-point mode
 |
 |        Input:  f__mhz            - Frequency, in MHz
 |                h_sys__meter      - Average height of the path above
 |                                    mean sea level, in meters
 |                N_0               - Refractivity, in N-Units
 |                pol               - Polarization
 |                                      + 0 : POLARIZATION__HORIZONTAL
 |                                      + 1 : POLARIZATION__VERTICAL
 |                epsilon           - Relative permittivity
 |                sigma             - Conductivity
 |
 |      Outputs:  Z_g               - Complex ground impedance
 |                gamma_e           - Curvature of the effective earth
 |                N_s               - Surface refractivity, in N-Units
 |
 |      Returns:  [None]
 |
 *===========================================================================*/		     

#include "GMS_complex_zmm16r4.hpp"

                    
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
                      void
                      itm_init_p2p_zmm16r4(const __m512 f_mhz,
		                           const __m512 h_sys_meter,
					   const __m512 N_0,
					   const int32_t pol,
					   const __m512 epsilon,
					   const __m512 sigma,
					   ZMM16c4 &Z_g,
					   __m512 &gamm_e,
					   __m512 &N_s); 

		     //
		     //  AVX512 double-precision implementation
		     //
		     //

		     
                      __ATTR_ALWAYS_INLINE__
		    
		      static inline
                      void
		      data_preload_v8() {

                          const volatile __m512d t0 = dl1[0];
			  const volatile __m512d t1 = dl1[1];
			  const volatile __m512d t2 = dl2[0];
			  const volatile __m512d t3 = dl2[1];
			  const volatile __m512d t4 = dl2[2];
			  const volatile __m512d t5 = dl2[3];
			  const volatile __m512d t6 = dl2[4];
			  const volatile __m512d t7 = dl2[5];
			  const volatile __m512d t8 = dl3[0];
			  const volatile __m512d t9 = dl3[1];
		      }


		    
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline
		      __m512d
		      itm_terrain_roughness_zmm8r8(const __m512d d_meter, // this variable shall be a negative
		                                   const __m512d delta_h_meter) {
                            const register __m512d _1 = _mm512_set1_pd(1.0);
                            // // [ERL 79-ITS 67, Eqn 3], with distance in meters instead of kilometers
			    const register __m512d t0 = _mm512_div_pd(d_meter,cl1[0]);
			    const register __m512d t1 = _mm512_sub_pd(_1,_mm512_mul_pd(,_mm512_exp_pd(t0)));
			    return (_mm512_mul_pd(delta_h_meter,t1));
		      }


		   
                     __ATTR_ALWAYS_INLINE__
		   
		     static inline
		     __m512d
		     itm_fresnel_integral_zmm8r8(const __m512d v2) {

                          
			  const register __m512d t0    = _mm512_mul_pd(dl2[1],v2);
			  const register __m512d c0    = _mm512_sub_pd(_mm512_fmadd_pd(dl2[2],
			                                              _mm512_sqrt_pd(v2),
								      dl2[3],t0));
			  const register __m512d c1    = _mm512_fmadd_pd(dl2[4],_mm512_log10_pd(v2),dl2[5]);
			  const __m512d fi;
			  const __mmask8 k = 0x0;
			  k  = _mm512_cmp_mask_pd(v2,dl2[0],_CMP_LT_OQ);
			  fi = _mm512_mask_blend_pd(k,c1,c0);
			  return (fi);
		    }


		   
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline
		      __m512d
		      itm_knife_edge_diffraction_zmm8r8( const __m512d d_meter,
		                                         const __m512d f_mhz,
						         const __m512d a_e_meter,
						         const __m512d theta_los,
						         const __m512d d_hzn_meter1,
						         const __m512d d_hzn_meter2) {
			    //const register __m512 _0_0795775   = _mm512_set1_ps(0.0795775F);
			    //const register __m512 _47_7        = _mm512_set1_ps(47.7F);
			    const register __m512d t0           = _mm512_mul_pd(dl3[0],_mm512_div_pd(f_mhz,dl3[1]));
			    const register __m512d d_ml_meter   = _mm512_sub_pd(d_hzn_meter1,d_hzn_meter2);
			    const register __m512d d_nlos_meter = _mm512_sub_pd(d_meter,d_ml_meter);
			    const register __m512d theta_nlos   = _mm512_sub_pd(_mm512_div_pd(d_meter,a_e_meter),theta_los);
			    const register __m512d t1           = _mm512_mul_pd(theta_nlos,theta_nlos);
			    const register __m512d t2           = _mm512_div_pd(d_nlos_meter,_mm512_add_pd(d_nlos_meter,d_hzn_meter1));
			    const register __m512d t3           = _mm512_div_pd(d_nlos_meter,_mm512_add_pd(d_nlos_meter,d_hzn_meter2));
			    const register __m512d v_1          = _mm512_mul_pd(t0,_mm512_mul_pd(_mm512_mul_ps(t1,d_hzn_meter1),t2));
			    const register __m512d v_2          = _mm512_mul_pd(t0,_mm512_mul_pd(_mm512_mul_ps(t1,d_hzn_meter2),t3));
			    const register __m512d A_k_db       = itm_fresnel_integral_zmm8r8(v_1)+
			                                          itm_fresnel_integral_zmm8r8(v_2);
			    return (A_k_db);
		      }


		   
                     __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      __m512d
		      itm_smooth_earth_diffraction_zmm8r8( const __m512d d_meter,
		                                           const __m512d f_mhz,
							   const __m512d a_e_meter,
							   const __m512d theta_los,
							   const __m512d d_hzn_meter1,
							   const __m512d d_hzn_meter2,
							   const __m512d d_h_e_meter1,
							   const __m512d d_h_e_meter2,
							   const __m512d Z_g_re,
							   const __m512d Z_g_im); 

		    
                      __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      __m512d
		      itm_height_function_zmm8r8(const __m512d x_km,
		                                  const __m512d K); 


		      
     }

}





















#endif /*__GMS_ITM_AVX512_H__*/
