
#ifndef __GMS_PLANT_DIELECTRIC_AVX512_HPP__
#define __GMS_PLANT_DIELECTRIC_AVX512_HPP__ 050120220920



namespace file_info {

     const unsigned int GMS_PLANT_DIELECTRIC_AVX512_MAJOR = 1;
     const unsigned int GMS_PLANT_DIELECTRIC_AVX512_MINOR = 1;
     const unsigned int GMS_PLANT_DIELECTRIC_AVX512_MICRO = 0;
     const unsigned int gGMS_PLANT_DIELECTRIC_AVX512_FULLVER =
       1000U*GMS_PLANT_DIELECTRIC_AVX512_MAJOR+100U*gGMS_PLANT_DIELECTRIC_AVX512_MINOR+
       10U*gGMS_PLANT_DIELECTRIC_AVX512_MICRO;
     const char * const GMS_PLANT_DIELECTRIC_AVX512_CREATION_DATE = "05-01-2022 09:20 +00200 (WED 05 JUN 2022 GMT+2)";
     const char * const GMS_PLANT_DIELECTRIC_AVX512_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_PLANT_DIELECTRIC_AVX512_SYNOPSIS      = "Vegetation dielectric states AVX512 vectorized."
}

#include "GMS_avx512c16f32.h"
#include "GMS_avx512c8f64.h"
#include "GMS_avx512vecf32.h"
#include "GMS_avx512vecf64.h"
#include "GMS_config.h"


namespace gms {


         namespace math {





                         
	              /*
                               This kernel operates on 16 leaves.
                          */
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      ZMM16c4
		      static inline
		      zmm16c4_leaf_dielectric(const AVX512Vec16 leaf_mg,
		                              const AVX512Vec16 leaf_rho,
					      const AVX512Vec16 leaf_dens,
					      const AVX512Vec16 leaf_tau,
					      const AVX512Vec16 water_tmp,
					      const AVX512Vec16 veg_tmp,
					      const AVX512Vec16 theta,
					      const AVX512Vec16 freq,  // frequency value is a scalar broadcast to ZMM register
					      const bool dry_dens) { 

                           if(dry_dens) {
                              return (zmm16c4_veg_dielectric_2(leaf_mg,
			                                        leaf_rho,
								veg_tmp,
								theta,
								freq));
			   }
			   else {
                              return (zmm16c4_veg_dielectric_1(leaf_mg,
			                                       veg_tmp,
							       theta,
							       freq));
			  }
		  }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      ZMM8c8
		      static inline
		      zmm8c8_leaf_dielectric( const AVX512Vec8 leaf_mg,
		                              const AVX512Vec8 leaf_rho,
					      const AVX512Vec8 leaf_dens,
					      const AVX512Vec8 leaf_tau,
					      const AVX512Vec8 water_tmp,
					      const AVX512Vec8 veg_tmp,
					      const AVX512Vec8 theta,
					      const AVX512Vec8 freq,  // frequency value is a scalar broadcast to ZMM register
					      const bool dry_dens) { 

                           if(dry_dens) {
                              return (zmm8c8_veg_dielectric_2( leaf_mg,
			                                        leaf_rho,
								veg_tmp,
								theta,
								freq));
			   }
			   else {
                              return (zmm8c8_veg_dielectric_1( leaf_mg,
			                                       veg_tmp,
							       theta,
							       freq));
			  }
		  }


		  

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      ZMM16c4
		      static inline
		      zmm16c4_veg_dielectric_2(const AVX512Vec16 mg,
		                               const AVX512Vec16 veg_rho,
					       const AVX512Vec16 tempC,
					       const AVX512Vec16 theta,
					       const AVX512Vec16 freq) {
                           ZMM16c4 x;
			   ZMM16c4 y;
                           ZMM16c4 e;
			   ZMM16c4 f;
			   ZMM16c4 g;
			   ZMM16c4 w;
			   ZMM16c4 result();
			   AVX512Vec16 mv;
			   AVX512Vec16 a;
			   AVX512Vec16 b;
			   AVX512Vec16 c;
			   AVX512Vec16 d;
			   AVX512Vec16 top;
			   AVX512Vec16 fn;
			   AVX512Vec16 en;
			   AVX512Vec16 ein;
			   AVX512Vec16 t0;
			   const ZMM16c4 ONE(1.0f,1.0f);
			   const AVX512Vec16 _1(1.0F);
			   const AVX512Vec16 _1_7(1.7F);
			   const AVX512Vec16 _3_20(3.20F);
			   const AVX512Vec16 _6_5(6.5F);
			   const AVX512Vec16 _1_1109e10(1.1109e-10f);
			   const AVX512Vec16 _n3_824e12(-3.824e-12f);
			   const AVX512Vec16 _6_938e14(6.938e-14f);
			   const AVX512Vec16 _5_096e16(5.096e-16f);
			   const AVX512Vec16 _0_82(0.82f);
			   const AVX512Vec16 _0_166(0.166f);
			   const AVX512Vec16 _1_09(1.09f);
			   const AVX512Vec16 _22_74(22.74f);
			   const AVX512Vec16 _31_4(31.4f);
			   const AVX512Vec16 _59_5(59.5f);
			   const AVX512Vec16 _0(0.0f);
			   const AVX512Vec16 _88_045(88.045f);
			   const AVX512Vec16 _n0_4147f(-0.4147f);
			   const AVX512Vec16 _6_295e4(6.295e-4f);
			   const AVX512Vec16 _1_075e5(1.075e-5f);
			   const AVX512Vec16 _0_707106781186548(0.707106781186548f);
			   const AVX512Vec16 _4_9(4.9f);
			   const AVX512Vec16 _2_9(2.9f);
			   const AVX512Vec16 _55_0(55.0f);
			   mv  = mg*veg_rho/(_1*(_1-veg_rho));
			   t0  = mv*mv;
			   a   = _1_7+_3_20*mv+_6_5*t0;
			   top = _1_1109e10+tempC*(_n3_824e12+tempC*
			         (_6_938e14-tempC*_5_096e16));
			   b   = mv*(_0_82*mv+_0_166);
			   fn  = _1/(top*_1_09);
			   AVX512Vec16 c0 = rad_freq/fn;
			   e   = ZMM16c4(_1.m_v16,c0.m_v16);
			   d   = _22_74;
			   c   = _31_4*t0/(_59_5*t0+_1);
			   c0  = d/freq;
			   f   = ZMM16c4(_0.m_v16,c0.m_v16);
			   en  = _88_045+tempC*(_n0_4147+tempC*(_6_295e4+
                                    tempC*_1_075e5));
			   ein = _4_9;
			   c0  = sqrt(freq/AVX512Vec16(0.18f));
			   w   = _0_707106781186548.m_v16*ONE*c0.m_v16;
			   g   = _1.m_v16*w;
			   c0  = a+b*(_4_9+(en-ein);
			   x   = c0.m_v16/(e-f);
			   y   = c.m_v16*(_2_9.m_v16+_55_0.m_v16/g);
			   result = x+y;
			   return (result);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      ZMM8c8
		      static inline
		      zmm8c8_veg_dielectric_2( const AVX512Vec8 mg,
		                               const AVX512Vec8 veg_rho,
					       const AVX512Vec8 tempC,
					       const AVX512Vec8 theta,
					       const AVX512Vec8 freq) {
                           ZMM8c8 x;
			   ZMM8c8 y;
                           ZMM8c8 e;
			   ZMM8c8 f;
			   ZMM8c8 g;
			   ZMM8c8 w;
			   ZMM8c8 result();
			   AVX512Vec8 mv;
			   AVX512Vec8 a;
			   AVX512Vec8 b;
			   AVX512Vec8 c;
			   AVX512Vec8 d;
			   AVX512Vec8 top;
			   AVX512Vec8 fn;
			   AVX512Vec8 en;
			   AVX512Vec8 ein;
			   AVX512Vec8 t0;
			   const ZMM8c8 ONE(1.0,1.0);
			   const AVX512Vec8 _1(1.0);
			   const AVX512Vec8 _1_7(1.7);
			   const AVX512Vec8 _3_20(3.20);
			   const AVX512Vec8 _6_5(6.5);
			   const AVX512Vec8 _1_1109e10(1.1109e-10);
			   const AVX512Vec8 _n3_824e12(-3.824e-12);
			   const AVX512Vec8 _6_938e14(6.938e-14);
			   const AVX512Vec8 _5_096e16(5.096e-16);
			   const AVX512Vec8 _0_82(0.82);
			   const AVX512Vec8 _0_166(0.166);
			   const AVX512Vec8 _1_09(1.09);
			   const AVX512Vec8 _22_74(22.74);
			   const AVX512Vec8 _31_4(31.4);
			   const AVX512Vec8 _59_5(59.5);
			   const AVX512Vec8 _0(0.0);
			   const AVX512Vec8 _88_045(88.045);
			   const AVX512Vec8 _n0_4147f(-0.4147);
			   const AVX512Vec8 _6_295e4(6.295e-4);
			   const AVX512Vec8 _1_075e5(1.075e-5);
			   const AVX512Vec8 _0_707106781186548(0.707106781186548);
			   const AVX512Vec8 _4_9(4.9);
			   const AVX512Vec8 _2_9(2.9);
			   const AVX512Vec8 _55_0(55.0);
			   mv  = mg*veg_rho/(_1*(_1-veg_rho));
			   t0  = mv*mv;
			   a   = _1_7+_3_20*mv+_6_5*t0;
			   top = _1_1109e10+tempC*(_n3_824e12+tempC*
			         (_6_938e14-tempC*_5_096e16));
			   b   = mv*(_0_82*mv+_0_166);
			   fn  = _1/(top*_1_09);
			   AVX512Vec8 c0 = rad_freq/fn;
			   e   = ZMM8c8(_1.m_v8,c0.m_v8);
			   d   = _22_74;
			   c   = _31_4*t0/(_59_5*t0+_1);
			   c0  = d/freq;
			   f   = ZMM8c8(_0.m_v8,c0.m_v8);
			   en  = _88_045+tempC*(_n0_4147+tempC*(_6_295e4+
                                    tempC*_1_075e5));
			   ein = _4_9;
			   c0  = sqrt(freq/AVX512Vec8(0.18));
			   w   = _0_707106781186548.m_v8*ONE*c0.m_v8;
			   g   = _1.m_v8*w;
			   c0  = a+b*(_4_9+(en-ein);
			   x   = c0.m_v8/(e-f);
			   y   = c.m_v8*(_2_9.m_v8+_55_0.m_v8/g);
			   result = x+y;
			   return (result);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      ZMM16c4
		      static inline
		      zmm16c4_veg_dielectric_1(const AVX512Vec16 mg,
		                               const AVX512Vec16 tempC,
					       const AVX512Vec16 theta,
					       const AVX512Vec16 freq) {

                           ZMM16c4 x;
			   ZMM16c4 y;
                           ZMM16c4 e;
			   ZMM16c4 f;
			   ZMM16c4 g;
			   ZMM16c4 w;
			   ZMM16c4 result();
			   AVX512Vec16 mv;
			   AVX512Vec16 a;
			   AVX512Vec16 b;
			   AVX512Vec16 c;
			   AVX512Vec16 d;
			   AVX512Vec16 top;
			   AVX512Vec16 fn;
			   AVX512Vec16 en;
			   AVX512Vec16 ein;
			   AVX512Vec16 t0;
			   const ZMM16c4 ONE(1.0f,1.0f);
			   const AVX512Vec16 _1_7(1.7f);
			   const AVX512Vec16 _0_74(0.74f);
			   const AVX512Vec16 _6_16(6.16f);
			   const AVX512Vec16 _1_1109e10(1.1109e-10f);
			   const AVX512Vec16 _n3_824e12(-3.824e-12f);
			   const AVX512Vec16 _6_938e14(6.938e-14f);
			   const AVX512Vec16 _5_096e16(5.096e-16f);
			   const AVX512Vec16 _0_55(0.55f);
			   const AVX512Vec16 _0_076(0.076f);
			   const AVX512Vec16 _1_09(1.09f);
			   const AVX512Vec16 _4_64(4.44f);
			   const AVX512Vec16 _7_36(7.36f);
			   const AVX512Vec16 _22_74(22.74f);
			   const AVX512Vec16 _0(0.0f);
			   const AVX512Vec16 _88_045(88.045f);
			   const AVX512Vec16 _n0_4147f(-0.4147f);
			   const AVX512Vec16 _6_295e4(6.295e-4f);
			   const AVX512Vec16 _1_075e5(1.075e-5f);
			   const AVX512Vec16 _0_707106781186548(0.707106781186548f);
			   const AVX512Vec16 _4_9(4.9f);
			   const AVX512Vec16 _2_9(2.9f);
			   const AVX512Vec16 _55_0(55.0f);
			   //!//mv  = mg*veg_rho/(_1*(_1-veg_rho));
			   t0  = mg*mg;
			   a   = _1_7-_0_74*mg+_6_16*t0;
			   top = _1_1109e10+tempC*(_n3_824e12+tempC*
			         (_6_938e14-tempC*_5_096e16));
			   b   = mg*(_0_55*mg-_0_076);
			   fn  = _1/(top*_1_09);
			   AVX512Vec16 c0 = rad_freq/fn;
			   e   = ZMM16c4(_1.m_v16,c0.m_v16);
			   d   = _22_74;
			   c   = _4_64*t0/(_7_36*t0+_1);
			   c0  = d/freq;
			   f   = ZMM16c4(_0.m_v16,c0.m_v16);
			   en  = _88_045+tempC*(_n0_4147+tempC*(_6_295e4+
                                    tempC*_1_075e5));
			   ein = _4_9;
			   c0  = sqrt(freq/AVX512Vec16(0.18f));
			   w   = _0_707106781186548.m_v16*ONE*c0.m_v16;
			   g   = _1.m_v16*w;
			   c0  = a+b*(_4_9+(en-ein);
			   x   = c0.m_v16/(e-f);
			   y   = c.m_v16*(_2_9.m_v16+_55_0.m_v16/g);
			   result = x+y;
			   return (result);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      ZMM16c4
		      static inline
		      zmm8c8_veg_dielectric_1( const AVX512Vec8 mg,
		                               const AVX512Vec8 tempC,
					       const AVX512Vec8 theta,
					       const AVX512Vec8 freq) {

                           ZMM8c8 x;
			   ZMM8c8 y;
                           ZMM8c8 e;
			   ZMM8c8 f;
			   ZMM8c8 g;
			   ZMM8c8 w;
			   ZMM8c8 result();
			   AVX512Vec8 mv;
			   AVX512Vec8 a;
			   AVX512Vec8 b;
			   AVX512Vec8 c;
			   AVX512Vec8 d;
			   AVX512Vec8 top;
			   AVX512Vec8 fn;
			   AVX512Vec8 en;
			   AVX512Vec8 ein;
			   AVX512Vec8 t0;
			   const ZMM8c8 ONE(1.0,1.0);
			   const AVX512Vec8 _1_7(1.7);
			   const AVX512Vec8 _0_74(0.74);
			   const AVX512Vec8 _6_16(6.16);
			   const AVX512Vec8 _1_1109e10(1.1109e-10);
			   const AVX512Vec8 _n3_824e12(-3.824e-12);
			   const AVX512Vec8 _6_938e14(6.938e-14);
			   const AVX512Vec8 _5_096e16(5.096e-16);
			   const AVX512Vec8 _0_55(0.55);
			   const AVX512Vec8 _0_076(0.076);
			   const AVX512Vec8 _1_09(1.09);
			   const AVX512Vec8 _4_64(4.44);
			   const AVX512Vec8 _7_36(7.36);
			   const AVX512Vec8 _22_74(22.74);
			   const AVX512Vec8 _0(0.0);
			   const AVX512Vec8 _88_045(88.045);
			   const AVX512Vec8 _n0_4147f(-0.4147);
			   const AVX512Vec8 _6_295e4(6.295e-4);
			   const AVX512Vec8 _1_075e5(1.075e-5);
			   const AVX512Vec8 _0_707106781186548(0.707106781186548);
			   const AVX512Vec8 _4_9(4.9);
			   const AVX512Vec8 _2_9(2.9);
			   const AVX512Vec8 _55_0(55.0);
			   //!//mv  = mg*veg_rho/(_1*(_1-veg_rho));
			   t0  = mg*mg;
			   a   = _1_7-_0_74*mg+_6_16*t0;
			   top = _1_1109e10+tempC*(_n3_824e12+tempC*
			         (_6_938e14-tempC*_5_096e16));
			   b   = mg*(_0_55*mg-_0_076);
			   fn  = _1/(top*_1_09);
			   AVX512Vec8 c0 = rad_freq/fn;
			   e   = ZMM16c4(_1.m_v8,c0.m_v8);
			   d   = _22_74;
			   c   = _4_64*t0/(_7_36*t0+_1);
			   c0  = d/freq;
			   f   = ZMM8c8(_0.m_v8,c0.m_v8);
			   en  = _88_045+tempC*(_n0_4147+tempC*(_6_295e4+
                                    tempC*_1_075e5));
			   ein = _4_9;
			   c0  = sqrt(freq/AVX512Vec8(0.18));
			   w   = _0_707106781186548.m_v8*ONE*c0.m_v8;
			   g   = _1.m_v8*w;
			   c0  = a+b*(_4_9+(en-ein);
			   x   = c0.m_v8/(e-f);
			   y   = c.m_v8*(_2_9.m_v8+_55_0.m_v8/g);
			   result = x+y;
			   return (result);
		   }


		   
		   

		   

		  


		  













		    
     } // math


} // gms














#endif /*__GMS_PLANT_DIELECTRIC_AVX512_HPP__*/
