


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




#include "GMS_itm_avx512.h"







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

		    
                      __ATTR_ALWAYS_INLINE__
		      	static inline
		      __m512 cabs_zmm16r4(const __m512 re,
		                              const __m512 im) {

                                 const __m512 tre = _mm512_mul_ps(re,re);
				 const __m512 tim = _mm512_mul_ps(im,im);
				return (_mm512_sqrt_ps(_mm512_add_ps(tre,tim)));
		     }

		      __ATTR_ALWAYS_INLINE__
		    
		      static inline
		      __m512d cabs_zmm8r8(const __m512d re,
		                          const __m512d im) {

                                 const __m512d tre = _mm512_mul_pd(re,re);
				 const __m512d tim = _mm512_mul_pd(im,im);
				return (_mm512_sqrt_pd(_mm512_add_pd(tre,tim)));
		     }
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
		      

		     
                     
		      __m512
		      gms::math::itm_smooth_earth_diffraction_zmm16r4(const __m512 d_meter,
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
			    const  __m512 G_x_db = _mm512_fmsub_ps(_0_05751,x_km[0],
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

                      
                    
		      __m512
		      gms::math::itm_height_function_zmm16r4(const __m512 x_km,
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

                    
                     
                      void
                      gms::math::itm_init_p2p_zmm16r4(const __m512 f_mhz,
		                           const __m512 h_sys_meter,
					   const __m512 N_0,
					   const int32_t pol,
					   const __m512 epsilon,
					   const __m512 sigma,
					   ZMM16c4 &Z_g,
					   __m512 &gamm_e,
					   __m512 &N_s) {

                            const __m512 gamma_a      = _mm512_set1_ps(0.0001569858712715855573F);
			    const __m512 _1           = _mm512_set1_ps(1.0F);
			    const __m512 _0           = _mm512_set1_ps(0.0F);
			    const __m512 _9460        = _mm512_set1_ps(9460.0);
			    const __m512 _179_3       = _mm512_set1_ps(179.3F);
			    const __m512 _0_4665      = _mm512_set1_ps(0.04665F);
			    const __m512 _18000       = _mm512_set1_ps(18000.0F);
			    const __m512 nh_sys_meter = _mm512_sub_ps(_0,h_sys_meter);
			    __mmask16 k = 0x0;
			    k = _mm512_cmp_mask_ps(h_sys_meter,_0,_CMP_EQ_OQ);
			    const  __m512 t0 = _mm512_mul_ps(N_0,
			                                     _mm512_exp_ps(_mm512_div_ps(nh_sys_meter,_9460)));
			    N_s = _mm512_mask_blend_ps(k,t0,N_0);
			    const  __m512 t1 = _mm512_sub_ps(_1,
			                                   _mm512_mul_ps(_0_4665,
							                 _mm512_exp_ps(_mm512_div_ps(N_s,_179_3))));
			    gamma_e = _mm512_mul_ps(gamma_a,t1);
			    const  ZMM16c4 ep_r = ZMM16c4{epsilon,
			                                       _mm512_mul_ps(_18000,
							                    _mm512_div_ps(sigma,f_mhz))};
			    Z_g = _mm512_sqrt_ps(_mm512_sub_ps(ep_r,_1));
			    if(pol == POLARIZATION_VERICAL)
			       Z_g /= ep_r;
		      }

		     //
		     //  AVX512 double-precision implementation
		     //
		     //

		     
                
   
                  
		      __m512d
		      gms::math::itm_smooth_earth_diffraction_zmm8r8( const __m512d d_meter,
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
			    const  __m512d G_x_db = _mm512_fmsub_pd(_0_05751,x_km[0],
			                                            _mm512_mul_pd(_10,_mm512_log10_pd(x_km[0])));
			    const __m512d x1 = _mm512_sub_pd(F_x_db[0],F_x_db[1]);
			    return (_mm512_sub_pd(G_x_db,_mm512_sub_pd(x1,_20)));
			 
		      }


		   
		      __m512d
		      gms::math::itm_height_function_zmm8r8(const __m512d x_km,
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



		      
 
