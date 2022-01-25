

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


		

      }

}






#endif /*__GMS_RADAR_TYPES_H__*/
