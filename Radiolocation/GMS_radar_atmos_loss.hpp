
#ifndef __GMS_RADAR_ATMOS_LOSS_HPP__
#define __GMS_RADAR_ATMOS_LOSS_HPP__ 090420221145


namespace file_info {

     const unsigned int GMS_RADAR_ATMOS_LOSS_MAJOR = 1;
     const unsigned int GMS_RADAR_ATMOS_LOSS_MINOR = 1;
     const unsigned int GMS_RADAR_ATMOS_LOSS_MICRO = 0;
     const unsigned int GMS_RADAR_ATMOS_LOSS_FULLVER =
       1000U*GMS_RADAR_ATMOS_LOSS_MAJOR+100U*GMS_RADAR_ATMOS_LOSS_MINOR+
       10U*GMS_RADAR_ATMOS_LOSS_MICRO;
     const char * const GMS_RADAR_ATMOS_LOSS_CREATION_DATE = "09-04-2022 11:45 +00200 (SAT 09 APR 2022 GMT+2)";
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
		     float blake_atmos_loss_r4_1(
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
                               
                           static const float __ATTR_ALIGN__(64) ser_n[48] = {2.0000000f,
                                                                              6.0000000f,
                                                                              12.0000000f,
                                                                              20.0000000f,
                                                                              30.0000000f,
                                                                              42.0000000f,
                                                                              56.0000000f,
                                                                              72.0000000f,
                                                                              90.0000000f,
                                                                              110.0000000f,
                                                                              132.0000000f,
                                                                              156.0000000f,
                                                                              182.0000000f,
                                                                              210.0000000f,
                                                                              240.0000000f,
                                                                              272.0000000f,
                                                                              306.0000000f,
                                                                              342.0000000f,
                                                                              380.0000000f,
                                                                              420.0000000f,
                                                                              462.0000000f,
                                                                              506.0000000f,
                                                                              552.0000000f,
                                                                              600.0000000f,
                                                                              650.0000000f,
                                                                              702.0000000f,
                                                                              756.0000000f,
                                                                              812.0000000f,
                                                                              870.0000000f,
                                                                              930.0000000f,
                                                                              992.0000000f,
                                                                              1056.0000000f,
                                                                              1122.0000000f,
                                                                              1190.0000000f,
                                                                              1260.0000000f,
                                                                              1332.0000000f,
                                                                              1406.0000000f,
                                                                              1482.0000000f,
                                                                              1560.0000000f,
                                                                              1640.0000000f,
                                                                              1722.0000000f,
                                                                              1806.0000000f,
                                                                              1892.0000000f,
                                                                              1980.0000000f,
                                                                              2070.0000000f,
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
			   //float Sigma_k    = 0.0f;
			   float gamma_k    = 0.0f; // dB/km, absorption coefficient
			   float S1         = 0.0f; // absorption loss integral coeff
			   float S2         = 0.0f; // absorption loss integral coeff
			   float S3         = 0.0f; // absorption loss integral coeff
			   float m_k        = 0.0f; // refractive index model
			   float gamma0     = 0.0f;
			   float h0         = 0.0f;
			   float m0         = 0.0f;
			   const float volatile u_plus_preload  = u_plus[0];
			   const float volatile u_minus_preload = u_minus[0];
			   const float volatile u_0_preload     = u_0[0];
			   const float volatile ser_n_preload   = ser_n[0];
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
                               const float delfk2 = delfk*delfk;
			       float Sigma_K      = 0.0f;
			       for(int32_t j = 1; j < 45; ++j) {
			           const float t0         = f_N_plus[j];
                                   const float Sig1_plus  = delfk/((t0-f_ghz)*(t0-f_ghz))+delfk2;
				   const float Sig2_plus  = delfk/((t0+f_ghz)*(t0+f_ghz))+delfk2;
				   const float t1         = f_N_minus[j];
				   const float Sig1_minus = delfk/((t1-f_ghz)*(t1-f_ghz))+delfk2;
				   const float Sig2_minus = delfk/((t1+f_ghz)*(t1+f_ghz))+delfk2;
				   const float F_N_plus   = Sig1_plus+Sig2_plus;
				   const float F_N_minus  = Sig1_minus+Sig2_minus;
				   const float t2         = u_0[j];
				   const float A1         = u_plus[j]*F_N_plus;
				   const float A2         = u_minus[j]*F_N_minus;
				   const float En         = 2.06844f*ser_n[j];
				   Sigma_k                += Z[j]*((A1+A2+t2*F0_k)*cephes_expf(-(En/T_k)));
			       }
			       const float T_k3   = 1.0f/(T_k*T_k*T_k);
			       if(i==0) {
                                  gamma0          = C*P_k*T_k3*fghz2*Sigma_k;
				  m0              = 1.0f+N*cephes_expf(c_e*h_k);
				  h0              = h_a*0.001f+k*delh;
			       }
			       gamma_k            = C*P_k*T_k3*fghz2*Sigma_k;
			       m_k                = 1.0f+N*cephes_expf(c_e*h_k);
			       
			       const float n0czth = n0*czth;
			       const float term1  = n0czth/(m_k*(1.0f+h_k/r_km));
			       const float term12 = term1*term1;
			       S2                 = gamma_k/cephes_sqrtf(1.0f-term12);
			       S3                 += S2;
			       
			   }
			   const float cterm = n0*czth;
			   const float term1 = cterm/(m0*(1.0f+h0/r_km));
			   const float term12= term1*term1;
			   S1 = gamma0/cephes_sqrtf(1.0f-term12);
			   return (2.0f*(S1+S2+2.0f*S3)*delh*0.5f);
			   
		   }


		    __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double blake_atmos_loss_r8_1(
		                                  const double h_a,     //m, heigth
						  const double f,      //Mhz, radar frequency
						  const double theta,  //deg, angle
						  const double R_m,    ////m, target range
						  const int32_t K) {   // Number of integration points (steps)

			   if(__builtin_expect(K<=1,0)) { return (0.0);}
			   //constants
			   static const double __ATTR_ALIGN__(64) u_plus[48] = {2.5,
			                                                       4.6666665,
									       6.7500000,
									       8.8000002,
			                                                       10.8333330,
									       12.8571424,
									       14.8750000,
									       16.8888893,
									       18.8999996,
									       20.9090900,
									       22.9166660,
									       24.9230766,
									       26.9285717,
									       28.9333324,
									       30.9375000,
									       32.9411774,
									       34.9444427,
									       36.9473686,
									       38.9500008,
									       40.9523811,
                                                                               42.9545441,
                                                                               44.9565201,
                                                                               46.9583321,
                                                                               48.9599991,
                                                                               50.9615402,
                                                                               52.9629631,
                                                                               54.9642868,
                                                                               56.9655190,
                                                                               58.9666672,
                                                                               60.9677429,
                                                                               62.9687500,
                                                                               64.9696960,
                                                                               66.9705887,
                                                                               68.9714279,
                                                                               70.9722214,
                                                                               72.9729767,
                                                                               74.9736862,
                                                                               76.9743576,
                                                                               78.9749985,
                                                                               80.9756088,
                                                                               82.9761887,
                                                                               84.9767456,
                                                                               86.9772720,
                                                                               88.9777756,
                                                                               90.9782639,
									       0.0,0.0,0.0};

		           static const double __ATTR_ALIGN__(64) u_minus[48] = {2.0000000,
                                                                                4.5000000,
                                                                                6.6666665,
                                                                                8.7500000,
                                                                                10.8000002,
                                                                                12.8333330,
                                                                                14.8571424,
                                                                                16.8750000,
                                                                                18.8888893,
                                                                                20.8999996,
                                                                                22.9090900,
                                                                                24.9166660,
                                                                                26.9230766,
                                                                                28.9285717,
                                                                                30.9333324,
                                                                                32.9375000,
                                                                                34.9411774,
                                                                                36.9444427,
                                                                                38.9473686,
                                                                                40.9500008,
                                                                                42.9523811,
                                                                                44.9545441,
                                                                                46.9565201,
                                                                                48.9583321,
                                                                                50.9599991,
                                                                                52.9615402,
                                                                                54.9629631,
                                                                                56.9642868,
                                                                                58.9655190,
                                                                                60.9666672,
                                                                                62.9677429,
                                                                                64.9687500,
                                                                                66.9696960,
                                                                                68.9705887,
                                                                                70.9714279,
                                                                                72.9722214,
                                                                                74.9729767,
                                                                                76.9736862,
                                                                                78.9743576,
                                                                                80.9749985,
                                                                                82.9756088,
                                                                                84.9761887,
                                                                                86.9767456,
                                                                                88.9772720,
                                                                                90.9777756,
										0.0,0.0,0.0};
                                
			    static const double __ATTR_ALIGN__(64) u_0[48] = {9.0000000,
                                                                             11.6666670,
                                                                             15.1666670,
                                                                             18.8999996,
                                                                             22.7333336,
                                                                             26.6190472,
                                                                             30.5357151,
                                                                             34.4722214,
                                                                             38.4222221,
                                                                             42.3818169,
                                                                             46.3484840,
                                                                             50.3205147,
                                                                             54.2967033,
                                                                             58.2761917,
                                                                             62.2583351,
                                                                             66.2426453,
                                                                             70.2287598,
                                                                             74.2163773,
                                                                             78.2052612,
                                                                             82.1952362,
                                                                             86.1861496,
                                                                             90.1778641,
                                                                             94.1702881,
                                                                             98.1633301,
                                                                             102.1569214,
                                                                             106.1509933,
                                                                             110.1455002,
                                                                             114.1403961,
                                                                             118.1356354,
                                                                             122.1311798,
                                                                             126.1270142,
                                                                             130.1231079,
                                                                             134.1194305,
                                                                             138.1159668,
                                                                             142.1127014,
                                                                             146.1096039,
                                                                             150.1066895,
                                                                             154.1039124,
                                                                             158.1012878,
                                                                             162.0987854,
                                                                             166.0964050,
                                                                             170.0941315,
                                                                             174.0919647,
                                                                             178.0899048,
                                                                             182.0879211,
									     0.0,0.0,0.0};
                               
                           static const double __ATTR_ALIGN__(64) ser_n[48] = {2.0000000,
                                                                              6.0000000,
                                                                              12.0000000,
                                                                              20.0000000,
                                                                              30.0000000,
                                                                              42.0000000,
                                                                              56.0000000,
                                                                              72.0000000,
                                                                              90.0000000,
                                                                              110.0000000,
                                                                              132.0000000,
                                                                              156.0000000,
                                                                              182.0000000,
                                                                              210.0000000,
                                                                              240.0000000,
                                                                              272.0000000,
                                                                              306.0000000,
                                                                              342.0000000,
                                                                              380.0000000,
                                                                              420.0000000,
                                                                              462.0000000,
                                                                              506.0000000,
                                                                              552.0000000,
                                                                              600.0000000,
                                                                              650.0000000,
                                                                              702.0000000,
                                                                              756.0000000,
                                                                              812.0000000,
                                                                              870.0000000,
                                                                              930.0000000,
                                                                              992.0000000,
                                                                              1056.000000,
                                                                              1122.0000000,
                                                                              1190.0000000,
                                                                              1260.0000000,
                                                                              1332.0000000,
                                                                              1406.0000000,
                                                                              1482.0000000,
                                                                              1560.0000000,
                                                                              1640.0000000,
                                                                              1722.0000000,
                                                                              1806.0000000,
                                                                              1892.0000000,
                                                                              1980.0000000,
                                                                              2070.0000000,
									      0.0,0.0,0.0};


									       
			                                                
			   constexpr double N    = 0.000313; //refractivity
			   constexpr double c_e  = -0.149;   // 1/km, decay constant
			   constexpr double n0   = 1.000313; // refractive index
			   constexpr double r_km = 6370.0;   // Earth radius
			   constexpr double r_m  = 6370000.0; // Earth radius
			   constexpr double a_e  = 8493333.0; // Earth effective radius
			   constexpr double p0   = 1013.25;   // atmospheric pressure constant
			   constexpr double alf1 = 5.2561222; // tropospheric model constants
			   constexpr double alf2 = 0.034164794;// as above
			   constexpr double alf3 = 11.388265;  // as above
			   constexpr double T0   = 300.0;      // standard temperature, K
			   constexpr double C    = 2.0058;     // absorption coeff const
			   constexpr double z    = 0.017453292519943295769236907685; // deg-to-rad (PI/180)
			   constexpr double zth  = z*theta;
			   constexpr double f_ghz= f*1000.0;
			   constexpr double fghz2= f_ghz*f_ghz;
			   const double czth = R_m*std::cos(zth);
			   const double h_m  = R_m*std::sin(zth)+((czth*czth)/(2.0*a_e));
			   const double h_km = h_m/1000.0;
			   const double delh = h_km/(double)K;
			   double T_k        = 0.0; //K, atmos temperature
			   double P_k        = 0.0; //mb, atmos pressure
			   double g_k        = 0.0; // line breadth constant parameter
			   //float Sigma_k    = 0.0f;
			   double gamma_k    = 0.0; // dB/km, absorption coefficient
			   double S1         = 0.0; // absorption loss integral coeff
			   double S2         = 0.0; // absorption loss integral coeff
			   double S3         = 0.0; // absorption loss integral coeff
			   double m_k        = 0.0; // refractive index model
			   double gamma0     = 0.0;
			   double h0         = 0.0;
			   double m0         = 0.0;
			   const double volatile u_plus_preload  = u_plus[0];
			   const double volatile u_minus_preload = u_minus[0];
			   const double volatile u_0_preload     = u_0[0];
			   const double volatile ser_n_preload   = ser_n[0];
			   for(int32_t i = 0; i < K; ++i) {
			       const double k    = (double)i;
			       
                               const double h_k  = h_a*0.001+k*delh; //km, current height
			       const double h_gm = r_m*h_k*1000.0/(r_m+h_k*1000.0); //m, geopotential altitude
			       const double h_gkm= h_gm*1000.0; //km, geopotential altitude
			       // Atmosphere temperature
			       if(h_gkm<=11.f) {
                                  T_k = 288.16-0.0065*h_k*1000.0;
			       }
			       else if(h_gkm>11.0 && h_gkm<25.0) {
                                  T_k = 216.66;
			       }
			       else {
                                  T_k = 216.66+0.003*(h_k-25.0)*1000.0;
			       }

			       if(h_gkm <= 11.0) {
                                  P_k = p0*std::pow(T_k*0.003470294280955024986118822876,alf1);
			       }
			       else if(h_gkm>11.0 && h_gkm<25.0) {
                                  const float t0 = 226.32/T_k;
				  P_k = t0*std::exp(-alf2*(h_k-11.0)*1000.0);
			       }
			       else {
                                  P_k = 24.886*std::pow(216.66/T_k,alf3);
			       }

			       if(h_k <= 8.0) {
                                  g_k = 0.640;
			       }
			       else if(h_k>8.0 && h_k<=25.0) {
                                  g_k = 0.640+0.04218*(h_k-8.0);
			       }
			       else {
                                  g_k = 1.357;
			       }

			       const double delfk = g_k*(P_k/p0)*(T0/T_k);  // line-breadth constant
			       const double F0_k  = delfk/((delfk*delfk)+fghz2); // nonresonant contribution
                               const double delfk2 = delfk*delfk;
			       double Sigma_K      = 0.0;
			       for(int32_t j = 1; j < 45; ++j) {
			           const double t0         = f_N_plusr8[j];
                                   const double Sig1_plus  = delfk/((t0-f_ghz)*(t0-f_ghz))+delfk2;
				   const double Sig2_plus  = delfk/((t0+f_ghz)*(t0+f_ghz))+delfk2;
				   const double t1         = f_N_minusr8[j];
				   const double Sig1_minus = delfk/((t1-f_ghz)*(t1-f_ghz))+delfk2;
				   const double Sig2_minus = delfk/((t1+f_ghz)*(t1+f_ghz))+delfk2;
				   const double F_N_plus   = Sig1_plus+Sig2_plus;
				   const double F_N_minus  = Sig1_minus+Sig2_minus;
				   const double t2         = u_0[j];
				   const double A1         = u_plus[j]*F_N_plus;
				   const double A2         = u_minus[j]*F_N_minus;
				   const double En         = 2.06844*ser_n[j];
				   Sigma_k                += Zr8[j]*((A1+A2+t2*F0_k)*std::exp(-(En/T_k)));
			       }
			       const double T_k3   = 1.0/(T_k*T_k*T_k);
			       if(i==0) {
                                  gamma0          = C*P_k*T_k3*fghz2*Sigma_k;
				  m0              = 1.0+N*std::exp(c_e*h_k);
				  h0              = h_a*0.001+k*delh;
			       }
			       gamma_k            = C*P_k*T_k3*fghz2*Sigma_k;
			       m_k                = 1.0+N*std::exp(c_e*h_k);
			       
			       const double n0czth = n0*czth;
			       const double term1  = n0czth/(m_k*(1.0+h_k/r_km));
			       const double term12 = term1*term1;
			       S2                 = gamma_k/std::sqrt(1.0-term12);
			       S3                 += S2;
			       
			   }
			   const double cterm = n0*czth;
			   const double term1 = cterm/(m0*(1.0+h0/r_km));
			   const double term12= term1*term1;
			   S1 = gamma0/std::sqrt(1.0f-term12);
			   return (2.0*(S1+S2+2.0*S3)*delh*0.5);
			   
		   }


    }


}


















#endif /*__GMS_RADAR_ATMOS_LOSS_HPP__*/
