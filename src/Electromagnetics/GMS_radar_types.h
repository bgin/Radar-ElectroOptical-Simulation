

#ifndef __GMS_RADAR_TYPES_H__
#define __GMS_RADAR_TYPES_H__


/*
     Aggregartion of various data types representing various  Radar components through its
     simulation processes and algorithms.
*/

#include "GMS_config.h"

namespace gms {

        namespace radiolocation {

	       // For simple jamming scenario simulations
               typedef struct __ATTR_ALIGN__(64) RadarParamAoS_R4_1 {

		           float gamm; //m, wavelength
			   float tf;   //sec, coherent processing time
			   float rho;  //usec,pulsewidth
			   float w;    //m, apperture width
                           float Kth;  //beamwidth constant
			   float Ln;   //pattern constant
			   float Ts;   //K, system noise temperature
			   float Fp;   //polarization factor
			   float La;   //dB, troposepheric attenuation
			   float F;    //radar pattern propagation factor
			   float Pt;   //kW, transmitter power
			   float tr;   //usec, PRI
			   float Lt;   //dB, transmitt line loss
			   float h;    //m, apperture height
			   float ha;   //m, phase centre
			   float Frdr; //range dependent response
			   float Dx;   //dB, detectibility factor
			   float Bt;   //Mhz, tuneable bandwidth
			   float Flen; //dB, tropospheric attenuation
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,52)
#endif
	       } RadarParamAoS_R4_1;


	         typedef struct __ATTR_ALIGN__(64) RadarParamAoS_R8_1 {

		           double gamm; //m, wavelength
			   double tf;   //sec, coherent processing time
			   double rho;  //usec,pulsewidth
			   double w;    //m, apperture width
                           double Kth;  //beanwidth constant
			   double Ln;   //pattern constant
			   double Ts;   //K, system noise temperature
			   double Fp;   //polarization factor
			   double La;   //dB, troposepheric attenuation
			   double F;    //radar pattern propagation factor
			   double Pt;   //kW, transmitter power
			   double tr;   //usec, PRI
			   double Lt;   //dB, transmitt line loss
			   double h;    //m, apperture height
			   double ha;   //m, phase centre
			   double Frdr; //range dependent response
			   double Dx;   //dB, detectibility factor
			   double Bt;   //Mhz, tuneable bandwidth
			   double Flen; //dB, tropospheric attenuation
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,40)
#endif
	       } RadarParamAoS_R8_1;


	       // Jammer and target parameters
	        typedef struct __ATTR_ALIGN__(64) JammerParamAoS_R4_1 {

                            float sig;  //m, RSC of target
			    float Pj;   //W, jammer power
			    float Gj;   //dB, jammer antenna gain
			    float Qj;   //dB, jammer noise quality
			    float Flenj;//dB, jammer lens factor
			    float Rj;   //km, jammer range
			    float Bj;   //Mhz,jammer noise BW
			    float Ltj;  //dB, jammer transmit loss
			    float Fpj;  //dB, jammer polarization
			    float Rmj;  //km, jammer screening range
			    float Fj;   //dB, jammer pattern factor of propagation
			    float Laj;  //dB, jammer troposhperic loss
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,8)
#endif			    
		} JammerParamAoS_R4_1;


		 typedef struct __ATTR_ALIGN__(64) JammerParamAoS_R8_1 {

                            double sig;  //m, RSC of target
			    double Pj;   //W, jammer power
			    double Gj;   //dB, jammer antenna gain
			    double Qj;   //dB, jammer noise quality
			    double Flenj;//dB, jammer lens factor
			    double Rj;   //km, jammer range
			    double Bj;   //Mhz,jammer noise BW
			    double Ltj;  //dB, jammer transmit loss
			    double Fpj;  //dB, jammer polarization
			    double Rmj;  //km, jammer screening range
			    double Fj;   //dB, jammer pattern factor of propagation
			    double Laj;  //dB, jammer troposhperic loss
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,32)
#endif			    
		} JammerParamAoS_R8_1;

#include <immintrin.h>
                 
                // SIMD types
		typedef struct __ATTR_ALIGN__(64) RadarParamSIMD_R4_16 {
                
                           __m512 gamm; //m, wavelength
			   __m512 tf;   //sec, coherent processing time
			   __m512 rho;  //usec,pulsewidth
			   __m512 w;    //m, apperture width
                           __m512 Kth;  //beamwidth constant
			   __m512 Ln;   //pattern constant
			   __m512 Ts;   //K, system noise temperature
			   __m512 Fp;   //polarization factor
			   __m512 La;   //dB, troposepheric attenuation
			   __m512 F;    //radar pattern propagation factor
			   __m512 Pt;   //kW, transmitter power
			   __m512 tr;   //usec, PRI
			   __m512 Lt;   //dB, transmitt line loss
			   __m512 h;    //m, apperture height
			   __m512 ha;   //m, phase centre
			   __m512 Frdr; //range dependent response
			   __m512 Dx;   //dB, detectibility factor
			   __m512 Bt;   //Mhz, tuneable bandwidth
			   __m512 Flen; //dB, tropospheric attenuation  
              } RadarParamSIMD_R4_16;


              typedef struct __ATTR_ALIGN__(64) RadarParamSIMD_R8_8 {
                
                           __m512d gamm; //m, wavelength
			   __m512d tf;   //sec, coherent processing time
			   __m512d rho;  //usec,pulsewidth
			   __m512d w;    //m, apperture width
                           __m512d Kth;  //beamwidth constant
			   __m512d Ln;   //pattern constant
			   __m512d Ts;   //K, system noise temperature
			   __m512d Fp;   //polarization factor
			   __m512d La;   //dB, troposepheric attenuation
			   __m512d F;    //radar pattern propagation factor
			   __m512d Pt;   //kW, transmitter power
			   __m512d tr;   //usec, PRI
			   __m512d Lt;   //dB, transmitt line loss
			   __m512d h;    //m, apperture height
			   __m512d ha;   //m, phase centre
			   __m512d Frdr; //range dependent response
			   __m512d Dx;   //dB, detectibility factor
			   __m512d Bt;   //Mhz, tuneable bandwidth
			   __m51d  Flen; //dB, tropospheric attenuation  
              } RadarParamSIMD_R8_8;


               // Jammer and target parameters

              typedef struct __ATTR_ALIGN__(64) JammerParamSIMD_R4_16 {

                            __m512 sig;  //m, RSC of target
			    __m512 Pj;   //W, jammer power
			    __m512 Gj;   //dB, jammer antenna gain
			    __m512 Qj;   //dB, jammer noise quality
			    __m512 Flenj;//dB, jammer lens factor
			    __m512 Rj;   //km, jammer range
			    __m512 Bj;   //Mhz,jammer noise BW
			    __m512 Ltj;  //dB, jammer transmit loss
			    __m512 Fpj;  //dB, jammer polarization
			    __m512 Rmj;  //km, jammer screening range
			    __m512 Fj;   //dB, jammer pattern factor of propagation
			    __m512 Laj;  //dB, jammer troposhperic loss
             } JammerParamSIMD_R4_16;


              typedef struct __ATTR_ALIGN__(64) JammerParamSIMD_R8_8 {

                            __m512d sig;  //m, RSC of target
			    __m512d Pj;   //W, jammer power
			    __m512d Gj;   //dB, jammer antenna gain
			    __m512d Qj;   //dB, jammer noise quality
			    __m512d Flenj;//dB, jammer lens factor
			    __m512d Rj;   //km, jammer range
			    __m512d Bj;   //Mhz,jammer noise BW
			    __m512d Ltj;  //dB, jammer transmit loss
			    __m512d Fpj;  //dB, jammer polarization
			    __m512d Rmj;  //km, jammer screening range
			    __m512d Fj;   //dB, jammer pattern factor of propagation
			    __m512d Laj;  //dB, jammer troposhperic loss
             } JammerParamSIMD_R8_8;

             // Same data types definitions as above, but for __m256 and __m256d.
             
             typedef struct __ATTR_ALIGN__(32) RadarParamSIMD_R4_8 {
                
                           __m256 gamm; //m, wavelength
			   __m256 tf;   //sec, coherent processing time
			   __m256 rho;  //usec,pulsewidth
			   __m256 w;    //m, apperture width
                           __m256 Kth;  //beamwidth constant
			   __m256 Ln;   //pattern constant
			   __m256 Ts;   //K, system noise temperature
			   __m256 Fp;   //polarization factor
			   __m256 La;   //dB, troposepheric attenuation
			   __m256 F;    //radar pattern propagation factor
			   __m256 Pt;   //kW, transmitter power
			   __m256 tr;   //usec, PRI
			   __m256 Lt;   //dB, transmitt line loss
			   __m256 h;    //m, apperture height
			   __m256 ha;   //m, phase centre
			   __m256 Frdr; //range dependent response
			   __m256 Dx;   //dB, detectibility factor
			   __m256 Bt;   //Mhz, tuneable bandwidth
			   __m256 Flen; //dB, tropospheric attenuation  
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,32)
#endif	
              } RadarParamSIMD_R4_8;


              typedef struct __ATTR_ALIGN__(32) RadarParamSIMD_R8_4 {
                
                           __m256d gamm; //m, wavelength
			   __m256d tf;   //sec, coherent processing time
			   __m256d rho;  //usec,pulsewidth
			   __m256d w;    //m, apperture width
                           __m256d Kth;  //beamwidth constant
			   __m256d Ln;   //pattern constant
			   __m256d Ts;   //K, system noise temperature
			   __m256d Fp;   //polarization factor
			   __m256d La;   //dB, troposepheric attenuation
			   __m256d F;    //radar pattern propagation factor
			   __m256d Pt;   //kW, transmitter power
			   __m256d tr;   //usec, PRI
			   __m256d Lt;   //dB, transmitt line loss
			   __m256d h;    //m, apperture height
			   __m256d ha;   //m, phase centre
			   __m256d Frdr; //range dependent response
			   __m256d Dx;   //dB, detectibility factor
			   __m256d Bt;   //Mhz, tuneable bandwidth
			   __m256d Flen; //dB, tropospheric attenuation   
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,32)
#endif	
              } RadarParamSIMD_R8_4;


               // Jammer and target parameters

              typedef struct __ATTR_ALIGN__(64) JammerParamSIMD_R4_8 {

                            __m256 sig;  //m, RSC of target
			    __m256 Pj;   //W, jammer power
			    __m256 Gj;   //dB, jammer antenna gain
			    __m256 Qj;   //dB, jammer noise quality
			    __m256 Flenj;//dB, jammer lens factor
			    __m256 Rj;   //km, jammer range
			    __m256 Bj;   //Mhz,jammer noise BW
			    __m256 Ltj;  //dB, jammer transmit loss
			    __m256 Fpj;  //dB, jammer polarization
			    __m256 Rmj;  //km, jammer screening range
			    __m256 Fj;   //dB, jammer pattern factor of propagation
			    __m256 Laj;  //dB, jammer troposhperic loss
             } JammerParamSIMD_R4_8;


              typedef struct __ATTR_ALIGN__(64) JammerParamSIMD_R8_4 {

                            __m256d sig;  //m, RSC of target
			    __m256d Pj;   //W, jammer power
			    __m256d Gj;   //dB, jammer antenna gain
			    __m256d Qj;   //dB, jammer noise quality
			    __m256d Flenj;//dB, jammer lens factor
			    __m256d Rj;   //km, jammer range
			    __m256d Bj;   //Mhz,jammer noise BW
			    __m256d Ltj;  //dB, jammer transmit loss
			    __m256d Fpj;  //dB, jammer polarization
			    __m256d Rmj;  //km, jammer screening range
			    __m256d Fj;   //dB, jammer pattern factor of propagation
			    __m256d Laj;  //dB, jammer troposhperic loss  
             } JammerParamSIMD_R8_4;



	    // Platform dependent errors
	    // Calculates the angle and range measurement errors cause
	    // by the errors in position and orientation of radar platform
	    typedef struct __ATTR_ALIGN__(64) PlatformErrAoS_R4_1 {

                    float R;      //nm, target range
		    float psi;    //deg,target azimuth
		    float theta;  //deg,target elevation
		    float da1;    //deg, platform error of yaw angle measurement
		    float da2;    //deg, platform error of pitch angle measurement
		    float da3;    //deg, platform error of roll angle measurement
		    float dx1;    //m,   platform center of gravity error x1-axis
		    float dx2;    //m,   platform center of gravity error x2-axis
		    float dx3;    //m,   platform center of gravity erorr x3-axis
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,28)
#endif			    
	    } PlatformErrAoS_R4_1;


	    typedef struct __ATTR_ALIGN__(64) PlatformErrAoS_R8_1 {

                    double R;      //nm, target range
		    double psi;    //deg,target azimuth
		    double theta;  //deg,target elevation
		    double da1;    //deg, platform error of yaw angle measurement
		    double da2;    //deg, platform error of pitch angle measurement
		    double da3;    //deg, platform error of roll angle measurement
		    double dx1;    //m,   platform center of gravity error x1-axis
		    double dx2;    //m,   platform center of gravity error x2-axis
		    double dx3;    //m,   platform center of gravity erorr x3-axis
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,56)
#endif			    
	    } PlatformErrAoS_R8_1;


	    // Propagation errors aggregating data type
	    typedef struct __ATTR_ALIGN__(64) PropagationErrAoS_R4_1 {

                    float R;             //nm, target range
		    float ht;            //m,  target height above surface
		    float ha;            //m,  antenna height above surface
		    float psi0;          //    Fresnel reflection coefficient
		    float psis;          //    specular scattering coefficient
		    float psiv;          //    vegetation coefficient
		    float beta;          //deg,surface slope
		    float th3;           //deg,elevation beamwidth (half-power, i.e. 3db)
		    float thmax;         //deg,beam axis elevation angle
		    float kme;           //    tracking error slope in elevation channel
		    float Ns;            //    refractivity at the radar site
		    int32_t flucts;      //    fluctuations switch, 1==low fluctuations,2==average fluctuations,3==high fluctuations in the -
		                         //    - refractive index
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,16)
#endif						 
	    } PropagationErrAoS_R4_1;


	      // Propagation errors aggregating data type
	    typedef struct __ATTR_ALIGN__(64) PropagationErrAoS_R8_1 {

                    double R;             //nm, target range
		    double ht;            //m,  target height above surface
		    double ha;            //m,  antenna height above surface
		    double psi0;          //    Fresnel reflection coefficient
		    double psis;          //    specular scattering coefficient
		    double psiv;          //    vegetation coefficient
		    double beta;          //deg,surface slope
		    double th3;           //deg,elevation beamwidth (half-power, i.e. 3db)
		    double thmax;         //deg,beam axis elevation angle
		    double kme;           //    tracking error slope in elevation channel
		    double Ns;            //    refractivity at the radar site
		    int32_t flucts;      //    fluctuations switch, 1==low fluctuations,2==average fluctuations,3==high fluctuations in the -
		                         //    - refractive index
#if (USE_STRUCT_PADDING) == 1			   
			   PAD_TO(1,36)
#endif						 
	    } PropagationErrAoS_R8_1;


	    // SIMD versions
	    typedef struct __ATTR_ALIGN__(64) PlatformErrSIMD_R4_16 {

                    __m512 R;      //nm, target range
		    __m512 psi;    //deg,target azimuth
		    __m512 theta;  //deg,target elevation
		    __m512 da1;    //deg, platform error of yaw angle measurement
		    __m512 da2;    //deg, platform error of pitch angle measurement
		    __m512 da3;    //deg, platform error of roll angle measurement
		    __m512 dx1;    //m,   platform center of gravity error x1-axis
		    __m512 dx2;    //m,   platform center of gravity error x2-axis
		    __m512 dx3;    //m,   platform center of gravity erorr x3-axis
    
	    } PlatformErrSIMD_R4_16;


	    typedef struct __ATTR_ALIGN__(64) PlatformErrSIMD_R8_8 {

                    __m512d R;      //nm, target range
		    __m512d psi;    //deg,target azimuth
		    __m512d theta;  //deg,target elevation
		    __m512d da1;    //deg, platform error of yaw angle measurement
		    __m512d da2;    //deg, platform error of pitch angle measurement
		    __m512d da3;    //deg, platform error of roll angle measurement
		    __m512d dx1;    //m,   platform center of gravity error x1-axis
		    __m512d dx2;    //m,   platform center of gravity error x2-axis
		    __m512d dx3;    //m,   platform center of gravity erorr x3-axis
    
	    } PlatformErrSIMD_R8_8;


	        // Propagation errors aggregating data type
	    typedef struct __ATTR_ALIGN__(64) PropagationErrSIMD_R4_16 {

	             // Part of SIMD variables will vector-expanded.
                    __m512 R;             //nm, target range
		    __m512 ht;            //m,  target height above surface
		    __m512 ha;            //m,  antenna height above surface
		    __m512 psi0;          //    Fresnel reflection coefficient
		    __m512 psis;          //    specular scattering coefficient
		    __m512 psiv;          //    vegetation coefficient
		    __m512 beta;          //deg,surface slope
		    __m512 th3;           //deg,elevation beamwidth (half-power, i.e. 3db)
		    __m512 thmax;         //deg,beam axis elevation angle
		    __m512 kme;           //    tracking error slope in elevation channel
		    __m512 Ns;            //    refractivity at the radar site
		    // int32_t flucts;      //    fluctuations switch, 1==low fluctuations,2==average fluctuations,3==high fluctuations in the -
		                          //    - refractive index
		    // Removed from the member list shall be passed as an function argument.			 
	    } PropagationErrSIMD_R4_16;


	          // Propagation errors aggregating data type
	    typedef struct __ATTR_ALIGN__(64) PropagationErrSIMD_R8_8 {

	            // Part of SIMD variables will vector-expanded.
                    __m512d R;             //nm, target range
		    __m512d ht;            //m,  target height above surface
		    __m512d ha;            //m,  antenna height above surface
		    __m512d psi0;          //    Fresnel reflection coefficient
		    __m512d psis;          //    specular scattering coefficient
		    __m512d psiv;          //    vegetation coefficient
		    __m512d beta;          //deg,surface slope
		    __m512d th3;           //deg,elevation beamwidth (half-power, i.e. 3db)
		    __m512d thmax;         //deg,beam axis elevation angle
		    __m512d kme;           //    tracking error slope in elevation channel
		    __m512d Ns;            //    refractivity at the radar site
		    // int32_t flucts;      //    fluctuations switch, 1==low fluctuations,2==average fluctuations,3==high fluctuations in the -
		                          //    - refractive index
		    // Removed from the member list shall be passed as an function argument.			 
	    } PropagationErrSIMD_R8_8;


	    // Athmospheric attenuation loss models (Blake Model)
	  /*  typedef struct __ATTR_ALIGN__(64) BlakeModelAos_R4_1 {

                    float theta;
		    float R;         //nm, target range in vacuum
		    float ha;        //m, heigth...?
		    float freq;      //Mhz, radar frequency
		    float Rkm;       //km, target range
		    float Rm;        //m,  target range
		    float hm;        //m,  target height above antenna
		    float hkm;       //km, target heigth
		    float dh;        //km, target height increment
		    float L1a;       //dB, 2-way atmospheric attenuation loss
		    int32_t K;       // number of points needed for approximation of absorption integral
		    
	    } BlakeModelAos_R4_1;*/
            // Oxygen resonance frequencies (Blake Model usage) 
            const float __ATTR_ALIGN__(64) f_N_plus[48] = {56.2648f,56.2648f,58.4466f,58.4466f,
	                                                   59.5910f,59.5910f,60.4348f,60.4348f,
							   61.1506f,61.1506f,61.8002f,61.8002f,
							   62.4212f,62.4212f,62.9980f,62.9980f,
							   63.5658f,63.5658f,64.1272f,64.1272f,
							   64.6779f,64.6779f,65.2240f,65.2240f,
							   65.7626f,65.7626f,66.2978f,66.2978f,
							   66.8313f,66.8313f,67.3627f,67.3627f,
							   67.8923f,67.8923f,68.4205f,68.4205f,
							   68.9478f,68.9478f,69.4741f,69.4741f,
							   70.0f,70.0f,70.5249f,70.55249f,
							   71.0497f,71.0497f,0.0f,0.0f,0.0f};

	    const float __ATTR_ALIGN__(64) f_N_minus[48] = {118.7505f,118.7505f,62.4862f,62.4862f,
	                                                    60.3061f,60.3061f,59.1642f,59.1642f,
							    58.3239f,58.3239f,57.6125f,57.6125f,
							    56.9682f,56.9682f,56.3634f,56.3634f,
							    55.7839f,55.7839f,55.2214f,55.2214f,
							    54.6728f,54.6728f,54.1294f,54.1294f,
							    53.5960f,53.5960f,53.0695f,53.0695f,
							    52.5458f,52.5458f,52.0259f,52.0259f,
							    51.5091f,51.5091f,50.9949f,50.9949f,
							    50.4830f,50.4830f,49.9730f,49.9730f,
							    49.4648f,49.4648f,48.9582f,48.9582f,
							    48.4530f,48.4530f,0.0f,0.0f,0.0f};

            const float __ATTR_ALIGN__(64) Z[48]         = {1.0f,0.0f,1.0f,0.0f,1.0f,0.0f,1.0f,
                                                            0.0f,1.0f,0.0f,1.0f,0.0f,1.0f,0.0f,
                                                            1.0f,0.0f,1.0f,0.0f,1.0f,0.0f,1.0f,
                                                            0.0f,1.0f,0.0f,1.0f,0.0f,1.0f,0.0f,
                                                            1.0f,0.0f,1.0f,0.0f,1.0f,0.0f,1.0f,
                                                            0.0f,1.0f,0.0f,1.0f,0.0f,1.0f,0.0f,
                                                            1.0f,0.0f,1.0f,0.0f,0.0f,0.0f};

            const float __ATTR_ALIGN__(64) f_N_plusr8[48] = {56.2648,56.2648,58.4466,58.4466,
	                                                   59.5910,59.5910,60.4348,60.4348,
							   61.1506,61.1506,61.8002,61.8002,
							   62.4212,62.4212,62.9980,62.9980,
							   63.5658,63.5658,64.1272,64.1272,
							   64.6779,64.6779,65.2240,65.2240,
							   65.7626,65.7626,66.2978,66.2978,
							   66.8313,66.8313,67.3627,67.3627,
							   67.8923,67.8923,68.4205,68.4205,
							   68.9478,68.9478,69.4741,69.4741,
							   70.0,70.0,70.5249,70.55249,
							   71.0497,71.0497,0.0,0.0,0.0};

	    const float __ATTR_ALIGN__(64) f_N_minusr8[48] = {118.7505,118.7505,62.4862,62.4862,
	                                                    60.3061,60.3061,59.1642,59.1642,
							    58.3239,58.3239,57.6125,57.6125,
							    56.9682,56.9682,56.3634,56.3634,
							    55.7839,55.7839,55.2214,55.2214,
							    54.6728,54.6728,54.1294,54.1294,
							    53.5960,53.5960,53.0695,53.0695,
							    52.5458,52.5458,52.0259,52.0259,
							    51.5091,51.5091,50.9949,50.9949,
							    50.4830,50.4830,49.9730,49.9730,
							    49.4648,49.4648,48.9582,48.9582,
							    48.4530,48.4530,0.0,0.0,0.0};

            const float __ATTR_ALIGN__(64) Zr8[48]         = {1.0,0.0,1.0,0.0,1.0,0.0,1.0,
                                                            0.0,1.0,0.0,1.0,0.0,1.0,0.0,
                                                            1.0,0.0,1.0,0.0,1.0,0.0,1.0,
                                                            0.0,1.0,0.0,1.0,0.0,1.0,0.0,
                                                            1.0,0.0,1.0,0.0,1.0,0.0,1.0,
                                                            0.0,1.0,0.0,1.0,0.0,1.0,0.0,
                                                            1.0,0.0,1.0,0.0,0.0,0.0};
      }

}






#endif /*__GMS_RADAR_TYPES_H__*/
