
#ifndef __GMS_RADAR_ATMOS_LOSS_HPP__
#define __GMS_RADAR_ATMOS_LOSS_HPP__ 090420221145


namespace file_info {

     const unsigned int GMS_RADAR_ATMOS_LOSS_MAJOR = 1;
     const unsigned int GMS_RADAR_ATMOS_LOSS_MINOR = 1;
     const unsigned int GMS_RADAR_ATMOS_LOSS_MICRO = 0;
     const unsigned int GMS_RADAR_ATMOS_LOSS_FULLVER =
       1000U*GMS_RADAR_ATMOS_LOSS_MAJOR+100U*GMS_RADAR_ATMOS_LOSS_MINOR+
       10U*GMS_RADAR_ATMOS_LOSS_MICRO;
     const char * const GMS_RADAR_ATMOS_LOSS_CREATION_DATE = "09-04-2022 11:45 +00200 (SAT 04 APR 2022 GMT+2)";
     const char * const GMS_RADAR_ATMOS_LOSS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_RADAR_ATMOS_LOSS_SYNOPSIS      = "Radar atmospheric loss equations."

}


#include <cstdint>
#include <cmath> //for double precision
#include "GMS_cephes.h" // single precision
#include "GMS_config.h"
#include "GMS_radar_types.h"


namespace gms {

         namespace radiolocation {



	             __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float blake_atmos_loss_r4_1(const float R,      //nm, target range in vacuum
		                                 const float h_a,     //m, heigth
						 const float f,      //Mhz, radar frequency
						 const float theta,  //deg, angle
						 const float R_m,    ////m, target range
						 const int32_t K) {   // Number of integration points (steps)

			   if(__builtin_expect(K<=1,0)) { return (0.0f);}
			   //constants
			   static const float __ATTR_ALIGN__(64) u_plus[48] = {2.5f,
			                                                       4.6666665f,
									       6.7500000f,
									       8.8000002f,
			                                                       10.8333330f,
									       12.8571424f,
									       14.8750000f,
									       16.8888893f,
									       18.8999996f,
									       20.9090900f,
									       22.9166660f,
									       24.9230766f,
									       26.9285717f,
									       28.9333324f,
									       30.9375000f,
									       32.9411774f,
									       34.9444427f,
									       36.9473686f,
									       38.9500008f,
									       40.9523811f,
                                                                               42.9545441f,
                                                                               44.9565201f,
                                                                               46.9583321f,
                                                                               48.9599991f,
                                                                               50.9615402f,
                                                                               52.9629631f,
                                                                               54.9642868f,
                                                                               56.9655190f,
                                                                               58.9666672f,
                                                                               60.9677429f,
                                                                               62.9687500f,
                                                                               64.9696960f,
                                                                               66.9705887f,
                                                                               68.9714279f,
                                                                               70.9722214f,
                                                                               72.9729767f,
                                                                               74.9736862f,
                                                                               76.9743576f,
                                                                               78.9749985f,
                                                                               80.9756088f,
                                                                               82.9761887f,
                                                                               84.9767456f,
                                                                               86.9772720f,
                                                                               88.9777756f,
                                                                               90.9782639f,
									       0.0f,0.0f,0.0f};

		           static const float __ATTR_ALIGN__(64) u_minus[48] = {2.0000000f,
                                                                                4.5000000f,
                                                                                6.6666665f,
                                                                                8.7500000f,
                                                                                10.8000002f,
                                                                                12.8333330f,
                                                                                14.8571424f,
                                                                                16.8750000f,
                                                                                18.8888893f,
                                                                                20.8999996f,
                                                                                22.9090900f,
                                                                                24.9166660f,
                                                                                26.9230766f,
                                                                                28.9285717f,
                                                                                30.9333324f,
                                                                                32.9375000f,
                                                                                34.9411774f,
                                                                                36.9444427f,
                                                                                38.9473686f,
                                                                                40.9500008f,
                                                                                42.9523811f,
                                                                                44.9545441f,
                                                                                46.9565201f,
                                                                                48.9583321f,
                                                                                50.9599991f,
                                                                                52.9615402f,
                                                                                54.9629631f,
                                                                                56.9642868f,
                                                                                58.9655190f,
                                                                                60.9666672f,
                                                                                62.9677429f,
                                                                                64.9687500f,
                                                                                66.9696960f,
                                                                                68.9705887f,
                                                                                70.9714279f,
                                                                                72.9722214f,
                                                                                74.9729767f,
                                                                                76.9736862f,
                                                                                78.9743576f,
                                                                                80.9749985f,
                                                                                82.9756088f,
                                                                                84.9761887f,
                                                                                86.9767456f,
                                                                                88.9772720f,
                                                                                90.9777756f,
										0.0f,0.0f,0.0f};
                                
			    static const float __ATTR_ALIGN__(64) u_0[48] = {9.0000000f,
                                                                             11.6666670f,
                                                                             15.1666670f,
                                                                             18.8999996f,
                                                                             22.7333336f,
                                                                             26.6190472f,
                                                                             30.5357151f,
                                                                             34.4722214f,
                                                                             38.4222221f,
                                                                             42.3818169f,
                                                                             46.3484840f,
                                                                             50.3205147f,
                                                                             54.2967033f,
                                                                             58.2761917f,
                                                                             62.2583351f,
                                                                             66.2426453f,
                                                                             70.2287598f,
                                                                             74.2163773f,
                                                                             78.2052612f,
                                                                             82.1952362f,
                                                                             86.1861496f,
                                                                             90.1778641f,
                                                                             94.1702881f,
                                                                             98.1633301f,
                                                                             102.1569214f,
                                                                             106.1509933f,
                                                                             110.1455002f,
                                                                             114.1403961f,
                                                                             118.1356354f,
                                                                             122.1311798f,
                                                                             126.1270142f,
                                                                             130.1231079f,
                                                                             134.1194305f,
                                                                             138.1159668f,
                                                                             142.1127014f,
                                                                             146.1096039f,
                                                                             150.1066895f,
                                                                             154.1039124f,
                                                                             158.1012878f,
                                                                             162.0987854f,
                                                                             166.0964050f,
                                                                             170.0941315f,
                                                                             174.0919647f,
                                                                             178.0899048f,
                                                                             182.0879211f,
									     0.0f,0.0f,0.0f};
                               



									       
			                                                
			   constexpr float N    = 0.000313f; //refractivity
			   constexpr float c_e  = -0.149f;   // 1/km, decay constant
			   constexpr float n0   = 1.000313f; // refractive index
			   constexpr float r_km = 6370.0f;   // Earth radius
			   constexpr float r_m  = 6370000.0f; // Earth radius
			   constexpr float a_e  = 8493333.0f; // Earth effective radius
			   constexpr float p0   = 1013.25f;   // atmospheric pressure constant
			   constexpr float alf1 = 5.2561222f; // tropospheric model constants
			   constexpr float alf2 = 0.034164794f;// as above
			   constexpr float alf3 = 11.388265f;  // as above
			   constexpr float T0   = 300.0f;      // standard temperature, K
			   constexpr float C    = 2.0058f;     // absorption coeff const
			   constexpr float z    = 0.017453292519943295769236907685f; // deg-to-rad (PI/180)
			   constexpr float zth  = z*theta;
			   constexpr float f_ghz= f*1000.0f;
			   constexpr float fghz2= f_ghz*f_ghz;
			   const float czth = R_m*cephes_cosf(zth);
			   const float h_m  = R_m*cephes_sinf(zth)+((czth*czth)/(2.0*a_e));
			   const float h_km = h_m/1000.0f;
			   const float delh = h_km/(float)K;
			   float T_k        = 0.0f; //K, atmos temperature
			   float P_k        = 0.0f; //mb, atmos pressure
			   float g_k        = 0.0f; // line breadth constant parameter
			   for(int32_t i = 0; i < K; ++i) {
			       const float k    = (float)i;
                               const float h_k  = h_a*0.001f+k*delh; //km, current height
			       const float h_gm = r_m*h_k*1000.0f/(r_m+h_k*1000.0f); //m, geopotential altitude
			       const float h_gkm= h_gm*1000.0f; //km, geopotential altitude
			       // Atmosphere temperature
			       if(h_gkm<=11.0f) {
                                  T_k = 288.16f-0.0065f*h_k*1000.0f;
			       }
			       else if(h_gkm>11.0f && h_gkm<25.0f) {
                                  T_k = 216.66f;
			       }
			       else {
                                  T_k = 216.66f+0.003f*(h_k-25.0f)*1000.0f;
			       }

			       if(h_gkm <= 11.0f) {
                                  P_k = p0*cephes_powf(T_k*0.003470294280955024986118822876f,alf1);
			       }
			       else if(h_gkm>11.0f && h_gkm<25.0f) {
                                  const float t0 = 226.32/T_k;
				  P_k = t0*cephes_expf(-alf2*(h_k-11.0f)*1000.0f);
			       }
			       else {
                                  P_k = 24.886f*cephes_powf(216.66f/T_k,alf3);
			       }

			       if(h_k <= 8.0f) {
                                  g_k = 0.640f;
			       }
			       else if(h_k>8 && h_k<=25.0f) {
                                  g_k = 0.640f+0.04218f*(h_k-8.0f);
			       }
			       else {
                                  g_k = 1.357f;
			       }

			       const float delfk = g_k*(P_k/p0)*(T0/T_k);  // line-breadth constant
			       const float F0_k  = delfk/((delfk*delfk)+fghz2); // nonresonant contribution

			       for(int32_t j = 1; j < 45; ++j) {
                                   
			       }
			      
			   }
		   }

    }


}


















#endif /*__GMS_RADAR_ATMOS_LOSS_HPP__*/
