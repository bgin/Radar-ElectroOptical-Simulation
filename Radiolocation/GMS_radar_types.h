

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

      }

}






#endif /*__GMS_RADAR_TYPES_H__*/
