
#ifndef __GMS_ITM_AVX512_HPP__
#define __GMS_ITM_AVX512_HPP__

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

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
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


		namespace {
		           // Cache-aligned const data.
                           static const __m512 cl1[2] = {_mm512_set1_ps(50.0e3F),
			                               	 _mm512_set1_ps(0.8F)};
			   static const __m512 cl2[6] = {_mm512_set1_ps(5.76),
						         _mm512_set1_ps(1.27F), 
                                                         _mm512_set1_ps(9.11F),
			                                 _mm512_set1_ps(6.02F),
							 _mm512_set1_ps(10.0F),
							 _mm512_set1_ps(12.953)};
			                         
			   static const __m512 cl3[2] = {_mm512_set1_ps(0.0795775F),
			                                  _mm512_set1_ps(47.7F)};

			   static const __m512d dl1[2] = {_mm512_set1_pd(50.0e3F),
			                               	 _mm512_set1_pd(0.8F)};
			   static const __m512d dl2[6] = {_mm512_set1_pd(5.76),
						         _mm512_set1_pd(1.27F), 
                                                         _mm512_set1_pd(9.11F),
			                                 _mm512_set1_pd(6.02F),
							 _mm512_set1_pd(10.0F),
							 _mm512_set1_pd(12.953)};
			                         
			   static const __m512d dl3[2] = {_mm512_set1_pd(0.0795775F),
			                                  _mm512_set1_pd(47.7F)};

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 cabs_zmm16r4(const __m512 re,
		                              const __m512 im) {

                                register const __m512 tre = _mm512_mul_ps(re,re);
				register const __m512 tim = _mm512_mul_ps(im,im);
				return (_mm512_sqrt_ps(_mm512_add_ps(tre,tim)));
		     }

		      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d cabs_zmm8r8(const __m512d re,
		                          const __m512d im) {

                                register const __m512d tre = _mm512_mul_pd(re,re);
				register const __m512d tim = _mm512_mul_pd(im,im);
				return (_mm512_sqrt_pd(_mm512_add_pd(tre,tim)));
		     }
		}


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
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


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
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

                     __ATTR_REGCALL__
                     __ATTR_ALWAYS_INLINE__
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
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

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
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
			    const register __m512 A_k_db       = itm_fresnel_integral_zmm16r4(v_1)+
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
		      

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
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
							   const __m512 Z_g_im) {
                           
                            __m512 a_meter[3];
			    __m512 d_km;
			    __m512 F_x_db[2];
			    __m512 K[3];
			    __m512 B_0[3];
			    __m512 x_km[3];
			    __m512 C_0[3];
			    const  __m512 _0_5       = _mm512_set1_ps(0.5F);
			    const  __m512 _0_75      = _mm512_set1_ps(0.75F);
			    const  __m512 _0_3       = _mm512_set1_ps(0.333333333333333333333333333333F);
			    const  __m512 t5         = _mm512_pow_ps(f_mhz,_0_3);
			    const  __m512 _n0_3      = _mm512_set1_ps(-0.33333333333333333333333333333F);
			    const  __m512 _1_607     = _mm512_set1_ps(1.607F);
			    const  __m512 _0_017778  = _mm512_set1_ps(0.017778F);
			    const  __m512 _0_001     = _mm512_set1_ps(0.001F);
			    const  __m512 _0_05751   = _mm512_set1_ps(0.05751F);
			    const  __m512 _10        = _mm512_set1_ps(10.0F);
			    const  __m512 _20        = _mm512_set1_ps(20.0F);
			    const  __m512 a_0_meter  = _mm512_set1_ps(6370.0e3F);
			    const  __m512 theta_nlos = _mm512_sub_ps(_mm512_div_ps(d_meter,a_e_meter),theta_los);
			    const  __m512 d_ML_meter = _mm512_add_ps(d_hzn_meter1,d_hzn_meter2);
			    const  __m512 t0         = _mm512_mul_ps(d_hzn_meter1,d_hzn_meter1);
			    const  __m512 t1         = _mm512_mul_ps(d_hzn_meter2,d_hzn_meter2);
			    const  __m512 t2         = _mm512_div_ps(_mm512_pow_ps(f_mhz,_n0_3),
								                cabs_zmm16r4(Z_g_re,Z_ge_im));
			    const          __m512 t3         = _mm512_mul_ps(d_hzn_meter1,_0_001);
			    const          __m512 t4         = _mm512_mul_ps(d_hzn_meter2,_0_001);
			    
			    // Compute 3 radii.
			     // C_0 is the ratio of the 4/3 earth to effective earth (technically Vogler 1964 ratio is 4/3 to effective earth k value), all raised to the (1/3) power.
                            // C_0 = (4 / 3k) ^ (1 / 3) [Vogler 1964, Eqn 2]
			    a_meter[0] = _mm512_div_ps(_mm512_sub_ps(d_meter,d_ML_meter),theta_nlos);
			    C_0[0]     = _mm512_pow_ps(_mm512_mul_ps(_0_75,
			                                          _mm512_div_ps(a_0_meter,a_meter[0])),_0_3);
			    K[0]       = _mm512_mul_ps(0_017778,_mm512_mul_ps(C_0[0],t2));
                            B_0[0]     = _mm512_sub_ps(_1_607,K[0]);
			    
			    
			    a_meter[1] = _mm512_mul_ps(_0_5,_mm512_div_ps(t0,d_h_e_meter1));
			    C_0[1]     = _mm512_pow_ps(_mm512_mul_ps(_0_75,
			                                          _mm512_div_ps(a_0_meter,a_meter[1])),_0_3);
			    K[1]       = _mm512_mul_ps(0_017778,_mm512_mul_ps(C_0[1],t2));
                            B_0[1]     = _mm512_sub_ps(_1_607,K[1]);
			    const __m512 c1 = _mm512_mul_ps(C_0[1],C_0[1]);
			    x_km[1]    = _mm512_mul_ps(B_0[1],_mm512_mul_ps(c1,_mm512_mul_ps(t5,t3)));
			    F_x_db[0]  = itm_height_function_zmm16r4(x_km[1],K[1]);
			    a_meter[2] = _mm512_mul_ps(_0_5,_mm512_div_ps(t1,d_h_e_meter2));
			    C_0[2]     = _mm512_pow_ps(_mm512_mul_ps(_0_75,
			                                          _mm512_div_ps(a_0_meter,a_meter[2])),_0_3);
			    K[2]       = _mm512_mul_ps(0_017778,_mm512_mul_ps(C_0[2],t2));
                            B_0[2]     = _mm512_sub_ps(_1_607,K[2]);
			    const __m512 c2 = _mm512_mul_ps(C_0[2],C_0[2]);
			    x_km[2]    = _mm512_mul_ps(B_0[2],_mm512_mul_ps(c2,_mm512_mul_ps(t5,t4)));
			    F_x_db[1]  = itm_height_function_zmm16r4(x_km[2],K[2]);
			    d_km       = _mm512_mul_ps(_mm512_mul_ps(a_meter[0],theta_nlos),_0_001);
			    const __m512 c0 = _mm512_mul_ps(C_0[0],C_0[0]);
			    const __m512 x0 = _mm512_add_ps(x_km[1],x_km[2]);
			    x_km[0]    = _mm512_mul_ps(B_0[0],_mm512_mul_ps(c0,_mm512_add_ps(d_km,x0)));
			     // compute height gain functions
			    const register __m512 G_x_db = _mm512_fmsub_ps(_0_05751,x_km[0],
			                                            _mm512_mul_ps(_10,_mm512_log10_ps(x_km[0])));
			    const __m512 x1 = _mm512_sub_ps(F_x_db[0],F_x_db[1]);
			    return (_mm512_sub_ps(G_x_db,_mm512_sub_ps(x1,_20)));
			 
		      }

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

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      itm_height_function_zmm16r4(const __m512 x_km,
		                                  const __m512 K) {

                            // Branch x_km < 200 was removed
			    // Negative -1.0 will be returned.
			    __mmask16 k = 0x0;
			    k = _mm512_cmp_mask_ps(x_km,_mm512_set1_ps(200.0F),_CMP_LT_OQ);
			    if(1==k) { return(_mm512_set1_ps(-1.0F));}
			    const __m512 _0_05751 = _mm512_set1_ps(0.05751F);
			    const __m512 _4_343   = _mm512_set1_ps(4.343F);
			    const __m512 _0_0134  = _mm512_set1_ps(0.0134F);
			    const __m512 _n0_005   = _mm512_set1_ps(-0.005F);
			    const __m512 _117     = _mm512_set1_ps(117.0F);
			    const __m512 _17_372  = _mm512_set1_ps(17.372F);
			    const __m512 _1       = _mm512_set1_ps(1.0F);
			    const __m512 lxkm     = _mm512_log_ps(x_km);
			    const __m512 c0       = _mm512_fmsub_ps(_17_372,lxkm,_117);
			    __m512 t0,t1,t2,t3,res;
			    __mmask16 k2 = 0x0;
			    k2  = _mm512_cmp_mask_ps(x_km,_mm512_set1_ps(2000.0F),_CMP_LT_OQ);
			    t0  = _mm512_fmsub_ps(_0_5751,x_km,_mm512_mul_ps(_4_343,lxkm));
			    t1  = _mm512_mul_ps(_0_134,_mm512_mul_ps(x_km,
			                                _mm512_exp_ps(_mm512_mul_ps(_n0_005,x_km)))); //w
			    t2  = _mm512_mul_ps(_mm512_fmadd_ps(_mm512_sub_ps(_1,t1),t0,t1),c0);
			    res = _mm512_mask_blend_ps(k2,t0,t2);
			    return (res);
		     }



		     //
		     //  AVX512 double-precision implementation
		     //
		     //

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
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


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
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


		     __ATTR_REGCALL__
                     __ATTR_ALWAYS_INLINE__
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
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


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
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


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
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
							   const __m512d Z_g_im) {
                           
                            __m512d a_meter[3];
			    __m512d d_km;
			    __m512d F_x_db[2];
			    __m512d K[3];
			    __m512d B_0[3];
			    __m512d x_km[3];
			    __m512d C_0[3];
			    const  __m512d _0_5       = _mm512_set1_pd(0.5);
			    const  __m512d _0_75      = _mm512_set1_pd(0.75);
			    const  __m512d _0_3       = _mm512_set1_pd(0.333333333333333333333333333333);
			    const          __m512d t5 = _mm512_pow_pd(f_mhz,_0_3);
			    const  __m512d _n0_3      = _mm512_set1_pd(-0.33333333333333333333333333333);
			    const  __m512d _1_607     = _mm512_set1_pd(1.607);
			    const  __m512d _0_017778  = _mm512_set1_pd(0.017778);
			    const  __m512d _0_001     = _mm512_set1_pd(0.001);
			    const  __m512d _0_05751   = _mm512_set1_pd(0.05751);
			    const  __m512d _10        = _mm512_set1_pd(10.0);
			    const  __m512d _20        = _mm512_set1_pd(20.0);
			    const  __m512d a_0_meter  = _mm512_set1_pd(6370.0e3);
			    const  __m512d theta_nlos = _mm512_sub_pd(_mm512_div_pd(d_meter,a_e_meter),theta_los);
			    const  __m512d d_ML_meter = _mm512_add_pd(d_hzn_meter1,d_hzn_meter2);
			    const  __m512d t0         = _mm512_mul_pd(d_hzn_meter1,d_hzn_meter1);
			    const  __m512d t1         = _mm512_mul_pd(d_hzn_meter2,d_hzn_meter2);
			    const  __m512d t2         = _mm512_div_pd(_mm512_pow_pd(f_mhz,_n0_3),
								                cabs_zmm8r8(Z_g_re,Z_ge_im));
			    const          __m512d t3         = _mm512_mul_pd(d_hzn_meter1,_0_001);
			    const          __m512d t4         = _mm512_mul_pd(d_hzn_meter2,_0_001);
			    
			    // Compute 3 radii.
			     // C_0 is the ratio of the 4/3 earth to effective earth (technically Vogler 1964 ratio is 4/3 to effective earth k value), all raised to the (1/3) power.
                            // C_0 = (4 / 3k) ^ (1 / 3) [Vogler 1964, Eqn 2]
			    a_meter[0] = _mm512_div_pd(_mm512_sub_pd(d_meter,d_ML_meter),theta_nlos);
			    C_0[0]     = _mm512_pow_pd(_mm512_mul_pd(_0_75,
			                                          _mm512_div_pd(a_0_meter,a_meter[0])),_0_3);
			    K[0]       = _mm512_mul_pd(0_017778,_mm512_mul_pd(C_0[0],t2));
                            B_0[0]     = _mm512_sub_pd(_1_607,K[0]);
			    
			    
			    a_meter[1] = _mm512_mul_pd(_0_5,_mm512_div_pd(t0,d_h_e_meter1));
			    C_0[1]     = _mm512_pow_pd(_mm512_mul_ps(_0_75,
			                                          _mm512_div_pd(a_0_meter,a_meter[1])),_0_3);
			    K[1]       = _mm512_mul_pd(0_017778,_mm512_mul_pd(C_0[1],t2));
                            B_0[1]     = _mm512_sub_pd(_1_607,K[1]);
			    const __m512d c1 = _mm512_mul_pd(C_0[1],C_0[1]);
			    x_km[1]    = _mm512_mul_pd(B_0[1],_mm512_mul_pd(c1,_mm512_mul_pd(t5,t3)));
			    F_x_db[0]  = itm_height_function_zmm8r8(x_km[1],K[1]);
			    a_meter[2] = _mm512_mul_pd(_0_5,_mm512_div_pd(t1,d_h_e_meter2));
			    C_0[2]     = _mm512_pow_pd(_mm512_mul_pd(_0_75,
			                                          _mm512_div_pd(a_0_meter,a_meter[2])),_0_3);
			    K[2]       = _mm512_mul_pd(0_017778,_mm512_mul_pd(C_0[2],t2));
                            B_0[2]     = _mm512_sub_pd(_1_607,K[2]);
			    const __m512d c2 = _mm512_mul_pd(C_0[2],C_0[2]);
			    x_km[2]    = _mm512_mul_pd(B_0[2],_mm512_mul_pd(c2,_mm512_mul_pd(t5,t4)));
			    F_x_db[1]  = itm_height_function_zmm8r8(x_km[2],K[2]);
			    d_km       = _mm512_mul_pd(_mm512_mul_pd(a_meter[0],theta_nlos),_0_001);
			    const __m512d c0 = _mm512_mul_pd(C_0[0],C_0[0]);
			    const __m512d x0 = _mm512_add_pd(x_km[1],x_km[2]);
			    x_km[0]    = _mm512_mul_pd(B_0[0],_mm512_mul_pd(c0,_mm512_add_pd(d_km,x0)));
			     // compute height gain functions
			    const register __m512d G_x_db = _mm512_fmsub_pd(_0_05751,x_km[0],
			                                            _mm512_mul_pd(_10,_mm512_log10_pd(x_km[0])));
			    const __m512d x1 = _mm512_sub_pd(F_x_db[0],F_x_db[1]);
			    return (_mm512_sub_pd(G_x_db,_mm512_sub_pd(x1,_20)));
			 
		      }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
		      itm_height_function_zmm8r8(const __m512d x_km,
		                                  const __m512d K) {

                            // Branch x_km < 200 was removed
			    // Negative -1.0 will be returned.
			    __mmask8 k = 0x0;
			    k = _mm512_cmp_mask_pd(x_km,_mm512_set1_pd(200.0),_CMP_LT_OQ);
			    if(1==k) { return(_mm512_set1_pd(-1.0));}
			    const __m512d _0_05751 = _mm512_set1_pd(0.05751);
			    const __m512d _4_343   = _mm512_set1_pd(4.343);
			    const __m512d _0_0134  = _mm512_set1_pd(0.0134);
			    const __m512d _n0_005   = _mm512_set1_pd(-0.005);
			    const __m512d _117     = _mm512_set1_pd(117.0);
			    const __m512d _17_372  = _mm512_set1_pd(17.372);
			    const __m512d _1       = _mm512_set1_pd(1.0);
			    const __m512d lxkm     = _mm512_log_pd(x_km);
			    const __m512d c0       = _mm512_fmsub_pd(_17_372,lxkm,_117);
			    __m512d t0,t1,t2,res;
			    __mmask8 k2 = 0x0;
			    k2  = _mm512_cmp_mask_pd(x_km,_mm512_set1_pd(2000.0),_CMP_LT_OQ);
			    t0  = _mm512_fmsub_pd(_0_5751,x_km,_mm512_mul_pd(_4_343,lxkm));
			    t1  = _mm512_mul_pd(_0_134,_mm512_mul_pd(x_km,
			                                _mm512_exp_pd(_mm512_mul_pd(_n0_005,x_km)))); //w
			    t2  = _mm512_mul_pd(_mm512_fmadd_pd(_mm512_sub_pd(_1,t1),t0,t1),c0);
			    res = _mm512_mask_blend_pd(k2,t0,t2);
			    return (res);
		     }



		      
     }

}





















#endif /*__GMS_ITM_AVX512_HPP__*/
