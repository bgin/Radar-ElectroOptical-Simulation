
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

	    /*
                  Band     Freq(Ghz)    K_alf (clear air)           k_alf/r (rain)       k_alf/r (snow)
                  UHF       0.4            0.100                       0                  0
                  L         1.3            0.012                       0.0003             0.0003
                  S         3.0            0.015                       0.0013             0.0013
                  C         5.5            0.017                       0.008              0.008
                  X         10             0.24                        0.037              0.002
                  Ku        15             0.055                       0.083              0.004
                  K         22             0.030                       0.23               0.008
                  Ka        35             0.14                        0.57               0.015
                  V         60             35                          1.3                0.03
             */

	     // Simplified Barton atmos loss model

	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    float barton_atmos_loss_r4_1(const float theta, // rad, elevation angle
	                                 const float R_km ,  // km, target range
					 const float k_alf) { // attenuation coefficient (table above)

                   constexpr float z    = 0.017453292519943295769236907685f; // deg-to-rad (PI/180)
		   const float zth      = z*theta;
		   const float th_eff   = 0.00025f/(zth+0.028f);
		   const float R_eff    = 3.0f/cephes_sinf(th_eff);
		   const float arg      = R_km/R_eff;
		   const float arg2     = k_alf*R_eff;
		   return (arg2*(1.0f-cephes_expf(-arg))); // db, two-way atmos attenaytion loss
	    }


	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    double barton_atmos_loss_r8_1(const double theta, // rad, elevation angle
	                                  const double R_km ,  // km, target range
					  const double k_alf) { // attenuation coefficient (table above)

                   constexpr double z    = 0.017453292519943295769236907685; // deg-to-rad (PI/180)
		   const double zth      = z*theta;
		   const double th_eff   = 0.00025/(zth+0.028);
		   const double R_eff    = 3.0/std::sin(th_eff);
		   const double arg      = R_km/R_eff;
		   const double arg2     = k_alf*R_eff;
		   return (arg2*(1.0-std::exp(-arg))); // db, two-way atmos attenaytion loss
	    }

#include <omp.h>
#include "GMS_cquadpack.h"

// Integrand for beamshape_loss_r4 and beamshape_loss_r8
// functions
double
gauss_v_pattern(double t, void * __restrict data) {
     const double t_3db = *(double*)data;
     const double t0    = t/t_3db;
     const double vt = std::exp(-1.3863*t0*t0);
     return (vt*vt*vt*vt);
}
	   // Beamshape Loss
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void beamshape_loss_r4(const float theta, //deg, 3-dB beamwidth				 
                                   const float om_a,  //rpm, antenna rotation rate
				   const float PRF,   //Hz, pulse repetition frequency
				   const int32_t N,   //number of points to simulate antenna pattern
				   const float n,   //number of pulses required for detection
				   const int32_t level, //1->loss at 3-dB level, 2->loss at level specified by number of pulses required for detection
				   float * __restrict __ATTR_ALIGN__(64) f_n,
				   float & L1p       ) { //dB, beamshape loss

		     if(__builtin_expect(N<=10,0)) { return;}
		     float t0;
		     float S_f; // area covered by actual power pattern (two-way propagation)
		     float S_r; // area covered by ideal rectnagular pattern;
		     // Quadpack dqage arguments
		     float a,b,epsbas,epsrel,abserr;
		     int32_t key,neval,ier,last;
		     
		     const float t_s   = 60.0f/om_a; //s, single scan time/s [rotating antenna]
		     const float t_3db = t_s*theta*0.002777777777777777777777777778f; //s,time-on-target at 3dB beamwidth level
                     const float n_3db = t_3db*PRF; // number of pulses received at 3dB level
		     const float t_n   = n/PRF; //s, time-on-target required to receive n-pulses
		     // Antenna pattern simulation
		     const float Ts    = -5.0f*t_3db*0.5f; // pattern limits
		     const float Te    = 5.0f*t_3db*0.5f;  // as above
		     const float delt  = (Te-Ts)/N; //s, antenna pattern increment step
#if defined(__INTEL_COMPILER) || defined(__ICC)
                     __assume_aligned(f_n,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                     f_n = (float*)__builtin_assume_aligned(f_n,64);
#endif
#pragma omp simd simdlen(4) linear(i:1) private(ti,t0)
		     for(int32_t i; i != N; ++i) {
                         const float ti = (float)i;
			 const float t0 = (Ts+ti*delt)/t_3db;
			 f_n[i]         = cephes_expf(-1.3863f*t0*t0);
		     }
		     if(level==1) {
                        t0 = t_3db
		     }
		     else {
                        t0 = tn;
		     }
		     //Beamshape loss calculation
		     a = -100.0f*t0;
		     b =  100.0f*t0;
		     epsabs = 0.0f;
		     epsrel = 0.0000001f;
		     key    = 5;
		     S_f    = (float)dqage(gauss_v_pattern,(double)a,(double)b,
		                           (double)epsabs,(double)epsrel,key,
					   (double)&abserr,&neval,&ier,&last,&t_3db);
		     if(ier>0) {
                        std::cerr << "dqage failed: ier=" << ier << std::endl;
			return;
		     }
		     L1p = 10.0f*cephes_logf(S_r/S_f);
	     }


	       // Beamshape Loss
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void beamshape_loss_r8(const double theta, //deg, 3-dB beamwidth				 
                                   const double om_a,  //rpm, antenna rotation rate
				   const double PRF,   //Hz, pulse repetition frequency
				   const int32_t N,   //number of points to simulate antenna pattern
				   const double n,   //number of pulses required for detection
				   const int32_t level, //1->loss at 3-dB level, 2->loss at level specified by number of pulses required for detection
				   double * __restrict __ATTR_ALIGN__(64) f_n,
				   double & L1p       ) { //dB, beamshape loss

		     if(__builtin_expect(N<=10,0)) { return;}
		     double t0;
		     double S_f; // area covered by actual power pattern (two-way propagation)
		     double S_r; // area covered by ideal rectnagular pattern;
		     // Quadpack dqage arguments
		     double a,b,epsbas,epsrel,abserr;
		     int32_t key,neval,ier,last;
		     
		     const double t_s   = 60.0/om_a; //s, single scan time/s [rotating antenna]
		     const double t_3db = t_s*theta*0.002777777777777777777777777778f; //s,time-on-target at 3dB beamwidth level
                     const double n_3db = t_3db*PRF; // number of pulses received at 3dB level
		     const double t_n   = n/PRF; //s, time-on-target required to receive n-pulses
		     // Antenna pattern simulation
		     const double Ts    = -5.0*t_3db*0.5; // pattern limits
		     const double Te    = 5.0*t_3db*0.5;  // as above
		     const double delt  = (Te-Ts)/N; //s, antenna pattern increment step
#if defined(__INTEL_COMPILER) || defined(__ICC)
                     __assume_aligned(f_n,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                     f_n = (double*)__builtin_assume_aligned(f_n,64);
#endif
#pragma omp simd simdlen(8) linear(i:1) private(ti,t0)
		     for(int32_t i; i != N; ++i) {
                         const double ti = (double)i;
			 const double t0 = (Ts+ti*delt)/t_3db;
			 f_n[i]         = std::exp(-1.3863*t0*t0);
		     }
		     if(level==1) {
                        t0 = t_3db
		     }
		     else {
                        t0 = tn;
		     }
		     //Beamshape loss calculation
		     a = -100.0*t0;
		     b =  100.0*t0;
		     epsabs = 0.0;
		     epsrel = 0.00000001;
		     key    = 5;
		     S_f    = dqage(gauss_v_pattern,a,b,
		                           epsabs,epsrel,key,
					   &abserr,&neval,&ier,&last,&t_3db);
		     if(ier>0) {
                        std::cerr << "dqage failed: ier=" << ier << std::endl;
			return;
		     }
		     L1p = 10.0*std::log(S_r/S_f);
	     }


             
	    /*
                  Inverse normal distribution helper function
                  Adapted from: https://stackoverflow.com/questions/54944716/inverse-of-cumulative-normal-distribution-function-with-parameters
              */
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_PURE__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    float invnorm(const float  p,
	                  const float  mu,
			  const float  sigma) {

                 if(p<=0.0f || p>=1.0f) return(std::numeric_limits<float>::quiet_NaN());
                 const float q = p-0.5f;
		 float r,val;
		 if(std::abs(q)<=0.425f) {
                    r   = 0.180625f-q*q;
		    val = q * (((((((r * 2509.0809287301226727f +
                                     33430.575583588128105f) * r + 67265.770927008700853f) * r +
                                     45921.953931549871457f) * r + 13731.693765509461125f) * r +
                                     1971.5909503065514427f) * r + 133.14166789178437745f) * r +
                                     3.387132872796366608f) /
                              (((((((r * 5226.495278852854561f +
                                     28729.085735721942674f) * r + 39307.89580009271061f) * r +
                                     21213.794301586595867f) * r + 5394.1960214247511077f) * r +
                                     687.1870074920579083f) * r + 42.313330701600911252f) * r + 1);
		 }
		 else {
                      if(q>0.0f) {
                         r = 1.0f-p;
		      }
		      else {
                         r = p;
		      }
		      r = std::sqrt(-std::log(r));
		      
		      if(r<=5.0f) {
                          r   += -1.6f;
			  val =  (((((((r * 7.7454501427834140764e-4f +
                                      0.0227238449892691845833f) * r + 0.24178072517745061177f) *
                                      r + 1.27045825245236838258f) * r +
                                      3.64784832476320460504f) * r + 5.7694972214606914055f) *
                                      r + 4.6303378461565452959f) * r + 1.42343711074968357734f) / 
                                (((((((r *
                                      1.05075007164441684324e-9f + 5.475938084995344946e-4f) *
                                      r + 0.0151986665636164571966f) * r +
                                      0.14810397642748007459f) * r + 0.68976733498510000455f) *
                                      r + 1.6763848301838038494f) * r +
                                      2.05319162663775882187f) * r + 1);
               
		      }
		      else {
                        /* very close to 0.0 or 1.0 */
			  r   += -5.0f;
			  val = (((((((r * 2.01033439929228813265e-7f +
                                       2.71155556874348757815e-5f) * r +
                                       0.0012426609473880784386f) * r + 0.026532189526576123093f) *
                                       r + 0.29656057182850489123f) * r +
                                       1.7848265399172913358f) * r + 5.4637849111641143699f) *
                                       r + 6.6579046435011037772f) / 
                                 (((((((r *
                                       2.04426310338993978564e-15f + 1.4215117583164458887e-7f) *
                                       r + 1.8463183175100546818e-5f) * r +
                                       7.868691311456132591e-4f) * r + 0.0148753612908506148525f)
                                       * r + 0.13692988092273580531f) * r +
                                       0.59983220655588793769f) * r + 1);
		      }

		      if(q<0.0f) {
                         val = -val;
		      }
		 }

		 return (mu+sigma*val);
		 
	   }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_PURE__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    double invnorm(const double  p,
	                   const double  mu,
			   const double  sigma) {

                 if(p<=0.0 || p>=1.0) return(std::numeric_limits<double>::quiet_NaN());
                 const double q = p-0.5;
		 double r,val;
		 if(std::abs(q)<=0.425) {
                    r   = 0.180625-q*q;
		    val = q * (((((((r * 2509.0809287301226727 +
                                     33430.575583588128105) * r + 67265.770927008700853) * r +
                                     45921.953931549871457) * r + 13731.693765509461125) * r +
                                     1971.5909503065514427) * r + 133.14166789178437745) * r +
                                     3.387132872796366608) /
                              (((((((r * 5226.495278852854561 +
                                     28729.085735721942674) * r + 39307.89580009271061) * r +
                                     21213.794301586595867) * r + 5394.1960214247511077) * r +
                                     687.1870074920579083) * r + 42.313330701600911252) * r + 1);
		 }
		 else {
                      if(q>0.0) {
                         r = 1.0-p;
		      }
		      else {
                         r = p;
		      }
		      r = std::sqrt(-std::log(r));
		      
		      if(r<=5.0) {
                          r   += -1.6;
			  val =  (((((((r * 7.7454501427834140764e-4 +
                                      0.0227238449892691845833) * r + 0.24178072517745061177) *
                                      r + 1.27045825245236838258) * r +
                                      3.64784832476320460504) * r + 5.7694972214606914055) *
                                      r + 4.6303378461565452959) * r + 1.42343711074968357734) / 
                                (((((((r *
                                      1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                                      r + 0.0151986665636164571966) * r +
                                      0.14810397642748007459) * r + 0.68976733498510000455) *
                                      r + 1.6763848301838038494) * r +
                                      2.05319162663775882187) * r + 1);
               
		      }
		      else {
                        /* very close to 0.0 or 1.0 */
			  r   += -5.0;
			  val = (((((((r * 2.01033439929228813265e-7 +
                                       2.71155556874348757815e-5) * r +
                                       0.0012426609473880784386) * r + 0.026532189526576123093) *
                                       r + 0.29656057182850489123) * r +
                                       1.7848265399172913358) * r + 5.4637849111641143699) *
                                       r + 6.6579046435011037772) / 
                                 (((((((r *
                                       2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                                       r + 1.8463183175100546818e-5) * r +
                                       7.868691311456132591e-4) * r + 0.0148753612908506148525)
                                       * r + 0.13692988092273580531) * r +
                                       0.59983220655588793769) * r + 1);
		      }

		      if(q<0.0) {
                         val = -val;
		      }
		 }

		 return (mu+sigma*val);
		 
	   }
	   
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void integration_loss(const float P_fa, // false alarm probability
	                          const float P_d,  // detection probability
				  const float n,  // number of pulses integrated
				  float & D_cdb,    // dB, coherent detectibility factor
				  float & L_idb) {  // dB, integration loss
				  
                 constexpr float c0 = 9.2f;
	         constexpr float c1 = 10.0f;
		 constexpr float c2 = 1.0f;
		 constexpr float c3 = 0.5f;
		 float D_c,L_i,lterm,rterm,diff;
		 float nom,denom;
		 lterm = invnorm(c2-P_fa,0.0f,c2);
		 rterm = invnorm(c2-P_d,0.0f,c2);
		 diff  = lterm-rterm;
		 D_c   = c3*diff*dif;
		 D_cdb = c1*std::log10(D_c);
		 nom   = c2+std::sqrt(c2+((c0*n)/D_c));
		 denom = c2+std::sqrt(c2+(c0/D_c));
		 L_i   = nom/denom;
		 L_idb = c1*std::log10(L_i);
	  }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void integration_loss(const double P_fa, // false alarm probability
	                          const double P_d,  // detection probability
				  const double n,  // number of pulses integrated
				  double & D_cdb,    // dB, coherent detectibility factor
				  double & L_idb) {  // dB, integration loss
				  
                 constexpr double c0 = 9.2;
	         constexpr double c1 = 10.0;
		 constexpr double c2 = 1.0;
		 constexpr double c3 = 0.5;
		 double D_c,L_i,lterm,rterm,diff;
		 double nom,denom;
		 lterm = invnorm(c2-P_fa,0.0,c2);
		 rterm = invnorm(c2-P_d,0.0,c2);
		 diff  = lterm-rterm;
		 D_c   = c3*diff*dif;
		 D_cdb = c1*std::log10(D_c);
		 nom   = c2+std::sqrt(c2+((c0*n)/D_c));
		 denom = c2+std::sqrt(c2+(c0/D_c));
		 L_i   = nom/denom;
		 L_idb = c1*std::log10(L_i);
	  }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void fluctuation_loss(const float P_fa,   // false alarm probability  
	                          const float P_d,    // detection probability
				  const float n,      // number of pulses integrated
				  float & D0_db,      // dB, basic swerling detectibility factor
				  float & D1_db,      // dB, basic swerling 1 detectibility factor
				  float & Lf_1_db,    // dB, fluctuation loss for swerling 1 case target
				  float & Lf_db) {    // dB, fluctuation loss
                    // swerling cases i.e. 0,1,2,3,4, only case '1' is used
		    // Swerling 0 -> nonfluctuation target
		    // Swerling 1 -> fluctuation target
		    constexpr float c0 = 0.707106781186547524400844362105f;
		    constexpr float c1 = 0.5f;
		    constexpr float c2 = 1.0f;
		    constexpr float c3 = 10.0f;
		    constexpr float c4 = 0.03f;
		    float D0,D1,lterm,rterm,mul;
		    float den1,den2,nom1,nom2;
		    float den3,Lf_1;
		    lterm = std::sqrt(std::log(c2/P_fa))-c0;
		    rterm = invnorm(c2-P_d,0.0f,c2);
		    mul   = lterm*rterm;
		    D0    = mul*mul-c1;
		    D0_bd = c3*std::log10(D0);
		    D1    = std::log(P_fa)/std::log(P_d);
		    D1_db = c3*std::log10(D1);
		    Lf_1  = D1/D0;
		    Lf_1_db = c3*std::log10(Lf_1)*(c2+c4*std::log10(n));
		    Lf_db = Lf_1_db;
	    }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void fluctuation_loss(const double P_fa,   // false alarm probability  
	                          const double P_d,    // detection probability
				  const double n,      // number of pulses integrated
				  double & D0_db,      // dB, basic swerling detectibility factor
				  double & D1_db,      // dB, basic swerling 1 detectibility factor
				  double & Lf_1_db,    // dB, fluctuation loss for swerling 1 case target
				  double & Lf_db) {    // dB, fluctuation loss
                    // swerling cases i.e. 0,1,2,3,4, only case '1' is used
		    // Swerling 0 -> nonfluctuation target
		    // Swerling 1 -> fluctuation target
		    constexpr double c0 = 0.707106781186547524400844362105;
		    constexpr double c1 = 0.5;
		    constexpr double c2 = 1.0;
		    constexpr double c3 = 10.0;
		    constexpr double c4 = 0.03;
		    double D0,D1,lterm,rterm,mul;
		    double den1,den2,nom1,nom2;
		    double den3,Lf_1;
		    lterm = std::sqrt(std::log(c2/P_fa))-c0;
		    rterm = invnorm(c2-P_d,0.0,c2);
		    mul   = lterm*rterm;
		    D0    = mul*mul-c1;
		    D0_bd = c3*std::log10(D0);
		    D1    = std::log(P_fa)/std::log(P_d);
		    D1_db = c3*std::log10(D1);
		    Lf_1  = D1/D0;
		    Lf_1_db = c3*std::log10(Lf_1)*(c2+c4*std::log10(n));
		    Lf_db = Lf_1_db;
	    }
	    
#include "GMS_root_finding.hpp"

         

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void MIT_processing_loss(float (*f) (float x), // points to erf(-0.5f)-y function, where y = 2.0*Pd-1.0
	                             const float ax,       // lower interval of above function
				     const float bx,       // upper interval of above function (guess value)
	                             const float N,        // NUmber of MTI filter pulses
	                             const float Pfa,      // false alarm probability 
				     const float Pd,       // detection probability
				     const float n,        // number of pulses integrated
				     const int32_t swerling, // swerling case the case 2 is used
				     const int32_t canceler,
				     const int32_t stagger,
				     float & D0_db,         //dB, basic Swerling 0 detectibility factor
				     float & D1_db,         //dB, basic Awerling 1 detectibility factor
				     float & Lf_1_db,       //dB, fluctuation loss for swrling 1 case target
				     float & Lf_db,         //dB, fluctuation loss
				     float & L1f_db,        //dB, fluctuation loss corrected
				     float & Lf,            //dB, fluctuation loss
				     float & L1f,           //dB, fluctuation loss corrected
				     float & Lmti_a_db,     //dB, noise correlation loss
				     float & Lmti_b_db,     //dB, velocity response loss
				     float & Lmti_db)  {    //dB, MTI processing loss 

		  constexpr float c0 = 10.0f;
		  constexpr float c1 = 1.0f;
		  constexpr float c2 = 0.03f;
		  constexpr float c3 = 0.05f;
		  constexpr float c4 = 0.8f;
		  constexpr float c5 = 0.1f;
                  float u,z,D0,D1,term,invPfa,sqr;
		  float fax,fbx,Lf_1,Dsw,D1sw;
		  float a,num1,den1,num2,den2;
		  float radix,n1e,Lmti_a,Lmti_b;
		  float num3,den3,factor,num4,den4;
		  float xzero,fzero;
		  int32_t iflag;
		  xzero = 0.0f;
		  fzero = 0.0f;
		  fax   = f(ax); // erf(guess value)-2.0*Pd-1.0
		  fbx   = f(bx);
		  muller(f,ax,bx,fax,fbx,xzero,fzero,iflag);
		  if(iflag == -2) return; // failed to find the root
		  u = xzero;
		  z = -u;
		  invPfa = c1/Pfa;
		  term = std::sqrt(std::log(invPfa))-z;
		  sqr  = term*term;
		  D0   = sqr-0.5f;
		  D0_db= c0*std::log10(D0);
		  if(canceler==1)
		     a = 0.667f;
		  else if(canceler==2)
		     a = 0.514285714285714285714285714286f;
		  else if(canceler==3)
		     a = 0.425531914893617021276595744681f;
		  num1 = std::log(Pfa);
		  den1 = std::log(Pd);
		  D1   = num1/den1-c1;
		  D1_db= std::log10(D1);
		  Lf_1 = D1/D0;
		  Lf_1_db = c0*std::log10(Lf_1)*(c1+c2*std::log10(n));
		  n1e = a*n;
		  Lf_db = Lf_1_db/n;
		  L1f_db= Lf_1_db/n1e;
		  radix = Lf_db*c5;
		  Lf    = std::pow(c0,radix);
		  radix = L1f_db*c5;
		  L1f   = std::pow(c0,radix);
		  Dsw   = D0*Lf;
		  D1sw  = D0*L1f;
		  Lmti_a= D1sw/Dsw;
		  Lmti_a_db = c0*std::log10(Lmti_a);
		  if(stagger==1) {
                     float term1 = c3/((c1-Pd)*(c1-Pd));
		     factor      = std::pow(term1,N-c1);
		  }
		  else if(stagger==2) {
                     factor      = c5/(c1-Pd);
		  }
		  Lmti_b = c4+factor;
		  Lmti_b_db = c0*std::log10(Lmti_b);
		  Lmti_db   = Lmti_a_db+Lmti_b_db;
	    }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
	    __ATTR_ALWAYS_INLINE
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static
	    inline
	    void MIT_processing_loss(double (*f) (double x), // points to erf(-0.5f)-y function, where y = 2.0*Pd-1.0
	                             const double ax,       // lower interval of above function i.e. 0.0
				     const double bx,       // upper interval of above function (guess value) i.e. 0.5
	                             const double N,        // NUmber of MTI filter pulses
	                             const double Pfa,      // false alarm probability 
				     const double Pd,       // detection probability
				     const double n,        // number of pulses integrated
				     const int32_t swerling, // swerling case the case 2 is used
				     const int32_t canceler,
				     const int32_t stagger,
				     double & D0_db,         //dB, basic Swerling 0 detectibility factor
				     double & D1_db,         //dB, basic Awerling 1 detectibility factor
				     double & Lf_1_db,       //dB, fluctuation loss for swrling 1 case target
				     double & Lf_db,         //dB, fluctuation loss
				     double & L1f_db,        //dB, fluctuation loss corrected
				     double & Lf,            //dB, fluctuation loss
				     double & L1f,           //dB, fluctuation loss corrected
				     double & Lmti_a_db,     //dB, noise correlation loss
				     double & Lmti_b_db,     //dB, velocity response loss
				     double & Lmti_db)  {    //dB, MTI processing loss 

		  constexpr double c0 = 10.0;
		  constexpr double c1 = 1.0;
		  constexpr double c2 = 0.03;
		  constexpr double c3 = 0.05;
		  constexpr double c4 = 0.8;
		  constexpr double c5 = 0.1;
                  double u,z,D0,D1,term,invPfa,sqr;
		  double fax,fbx,Lf_1,Dsw,D1sw;
		  double a,num1,den1,num2,den2;
		  double radix,n1e,Lmti_a,Lmti_b;
		  double num3,den3,factor,num4,den4;
		  double xzero,fzero;
		  int32_t iflag;
		  xzero = 0.0;
		  fzero = 0.0;
		  fax   = f(ax); // erf(guess value)-2.0*Pd-1.0
		  fbx   = f(bx);
		  muller(f,ax,bx,fax,fbx,xzero,fzero,iflag);
		  if(iflag == -2) return; // failed to find the root
		  u = xzero;
		  z = -u;
		  invPfa = c1/Pfa;
		  term = std::sqrt(std::log(invPfa))-z;
		  sqr  = term*term;
		  D0   = sqr-0.5;
		  D0_db= c0*std::log10(D0);
		  if(canceler==1)
		     a = 0.667;
		  else if(canceler==2)
		     a = 0.514285714285714285714285714286;
		  else if(canceler==3)
		     a = 0.425531914893617021276595744681;
		  num1 = std::log(Pfa);
		  den1 = std::log(Pd);
		  D1   = num1/den1-c1;
		  D1_db= std::log10(D1);
		  Lf_1 = D1/D0;
		  Lf_1_db = c0*std::log10(Lf_1)*(c1+c2*std::log10(n));
		  n1e = a*n;
		  Lf_db = Lf_1_db/n;
		  L1f_db= Lf_1_db/n1e;
		  radix = Lf_db*c5;
		  Lf    = std::pow(c0,radix);
		  radix = L1f_db*c5;
		  L1f   = std::pow(c0,radix);
		  Dsw   = D0*Lf;
		  D1sw  = D0*L1f;
		  Lmti_a= D1sw/Dsw;
		  Lmti_a_db = c0*std::log10(Lmti_a);
		  if(stagger==1) {
                     float term1 = c3/((c1-Pd)*(c1-Pd));
		     factor      = std::pow(term1,N-c1);
		  }
		  else if(stagger==2) {
                     factor      = c5/(c1-Pd);
		  }
		  Lmti_b = c4+factor;
		  Lmti_b_db = c0*std::log10(Lmti_b);
		  Lmti_db   = Lmti_a_db+Lmti_b_db;
	    }


			      
    } 


}


















#endif /*__GMS_RADAR_ATMOS_LOSS_HPP__*/
