

#ifndef __GMS_RADAR_JAMMING_HPP__
#define __GMS_RADAR_JAMMING_HPP__ 250120220951


namespace file_info {

     const unsigned int GMS_RADAR_JAMMING_MAJOR = 1;
     const unsigned int GMS_RADAR_JAMMING_MINOR = 1;
     const unsigned int GMS_RADAR_JAMMING_MICRO = 0;
     const unsigned int GMS_RADAR_JAMMING_FULLVER =
       1000U*GMS_RADAR_JAMMING_MAJOR+100U*GMS_RADAR_JAMMING_MINOR+
       10U*GMS_RADAR_JAMMING_MICRO;
     const char * const GMS_RADAR_JAMMING_CREATION_DATE = "25-01-2022 09:51 +00200 (TUE 25 JAN 2022 GMT+2)";
     const char * const GMS_RADAR_JAMMING_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_RADAR_JAMMING_SYNOPSIS      = "Radar Jamming Equations."

}


#include <cstdint>
#include <cmath> //for double precision
#include "GMS_cephes.h" // single precision
#include "GMS_config.h"
#include "GMS_radar_types.h"

namespace gms {

       namespace  radiolocation {


                 
                   namespace {
                        // speed of light
                        const float  c4    = 299792458.0f;
			const double c8    = 299792458.0;
			// Standard temp (K)
			const float  T04   = 290.0f;
			const double T08   = 290.0;
			// Effective Earth radius
			const float  k_e4  = 1.33333f;
			const double k_e8  = 1.333333333333333333;
			// Boltzmann constant
			const float  k_B4  = 1.38064852e-23f;
			const double k_B8  = 1.38064852e-23;
			// Earth radius (m)
			const float  a_e4  = 6378388.0f;
			const double a_e8  = 6378388.0;
			// PI constants
			const float  PI4   = 3.1415926535897932384626f;
			const double PI8   = 3.1415926535897932384626;
			const float  _4PI4 = 12.5663706143591729538506f; //4*PI
			const double _4PI8 = 12.5663706143591729538506;
			
		 }

		     // useful (dB) conversion functions.
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float to_dB_r4_1(const float x) {

		            constexpr float tm30 = 0.000000000000000000000000000001f; //1.00000000317107685097105134714e-30
			    return (10.0f*cephes_log10f(x+tm30));
			    
		     }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double to_dB_r8_1(const double x) {

		            constexpr double tm30 = 0.000000000000000000000000000001; //1.00000000000000008333642060759e-30
			    return (10.0*std::log10(x+tm30));
			    
		     }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float from_dB_r4_1(const float x) {

                            return (cephes_powf(10.0f,0.1f*x));
		     }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double from_dB_r8_1(const double x) {

                             return (std::pow(10.0,0.1*x));
		     }


		     // Auxilliary formulae computations
		     // Data type: RadarParamAoS_R4_1, RadarParamsAoS_R8_1
		     
		     // Number of pulses integrated
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float n_pulses_integ_r4_1(const float xtf,      //the names of arguments is the same as corresponding data structure
		                               const float xtr) {    // i.e. the name of data structure variable is in this case: tf,tr

                          const float invtr = 1.0f/xtr;
			  return (xtf*invtr);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double n_pulses_integ_r8_1(const double xtf,      //the names of arguments is the same as corresponding data structure
		                                const double xtr) {    // i.e. the name of data structure variable is in this case: tf,tr

                          const double invtr = 1.0/xtr;
			  return (xtf*invtr);
		    }


		    // Radar duty cycle
		   
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float duty_cycle_r4_1(const float xrho,
		                           const float xtr) {

                           const float invtr = 1.0f/xtr;
			   return (xrho*invtr);
		     }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double duty_cycle_r8_1(const double xrho,
		                            const double xtr) {

                           const double invtr = 1.0/xtr;
			   return (xrho*invtr);
		     }


		     // Radar average power (W)
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float radar_avg_power_r4_1(const float xPt,
		                                const float Dc) { // duty cycle argument

			    return (xPt*Dc);
                     }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float radar_avg_power_r8_1(const double xPt,
		                                const double Dc) { // duty cycle argument

			    return (xPt*Dc);
                     }

                     // Azimuth beam-width
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float azimuth_bw_r4_1(const float xKth,
		                           const float xgamm,
					   const float xw) {

                            return (xKth*(xgamm/xw));
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double azimuth_bw_r8_1(const double xKth,
		                            const double xgamm,
					    const double xw) {

                            return (xKth*(xgamm/xw));
		    }


		    // Elevation beam-width
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float elevation_bw_r4_1(const float xKth,
		                             const float xgamm,
					     const float xh) {

                            return (xKth*(xgamm/xh));
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double elevation_bw_r8_1(const double xKth,
		                              const double xgamm,
					      const double xh) {

                            return (xKth*(xgamm/xh));
		    }

		    // Radar antenna gain
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float radar_ant_gain_r4_1(const float tha, //rad, azimuth bw
		                               const float the, //rad, elevation bw
					       const float xLn) {

                            const float den = tha*the*xLn;
			    return (_4PI4/den);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double radar_ant_gain_r8_1(const double tha, //rad, azimuth bw
		                                const double the, //rad, elevation bw
					        const double xLn) {

                            const double den = tha*the*xLn;
			    return (_4PI8/den);
		    }


		    //Radar noise density W/Hz
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float noise_density_r4_1(const float xTs) {

                             return (xTs*k_B4);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double noise_density_r8_1(const double xTs) {

                             return (xTs*k_B8);
		    }

                
		    

 

		     
		     // Effect of thermal noise on Radar range
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float thermal_noise_RR_r4_1(const RadarParamAoS_R4_1   &rp,
		                                 const JammerParamAoS_R4_1  &jp) {

		           /* const float dummy = rp.gamm;
			    const float rcs   = jp.sig;
		            const float dc    = duty_cycle_r4_1(rp.rho,
			                                      rp.tr);
                            const float Pav   = radar_avg_power_r4_1(rp.Pt,
			                                           dc);
			    const float ag    = radar_ant_gain_r4_1(azimuth_bw_r4_1(rp.Kth,
			                                                            rp.gamm,
										    rp.h),
								    elevation_bw_r4_1(rp.Kth,
								                      rp.gamm,
										      rp.w),
										      rp.Ln);
			    const float den   = 1984.4017075391884912304967f*rp.Dx*rp.Lt*rp.La;
			    const float t1    = rp.gamm*rp.gamm;
			    const float t2    = Pav*rp.tf;
			    const float t3    = ag*ag;
			    const float t4    = rcs*rp.Frdr*rp.Frdr*rp.Fp*rp.Fp;
			    const float t5    = rp.F*rp.F*rp.F*rp.F*rp.Flen*rp.Flen;
			    const float num   = t1*t2*t3*t4*t5;
			    return (cephes_powf(num/den,0.25f));*/
			    // More efficient implementation
			    const float rcs  = jp.sig;
			    const float xgam = rp.gamm;
			    const float xrho = rp.rho;
			    const float xtr  = rp.tr;
			    const float xPt  = rp.Pt;
			    const float xKth = rp.Kth;
			    const float xw   = rp.w;
			    const float xh   = rp.h;
			    const float xLn  = rp.Ln;
			    const float xDx  = rp.Dx;
			    const float xLt  = rp.Lt;
			    const float xLa  = rp.La;
			    const float xtf  = rp.tf;
			    const float xFrd = rp.Frdr;
			    const float xFp  = rp.Fp;
			    const float xF   = rp.F;
			    const float xFlen= rp.Flen;
			    const float ratio= 0.0f;
			    float range      = 0.0f;
			    const float dc    = duty_cycle_r4_1(xrho,
			                                        xtr);
                            const float Pav   = radar_avg_power_r4_1(xPt,
			                                             dc);
			    const float ag    = radar_ant_gain_r4_1(azimuth_bw_r4_1(xKth,
			                                                            xgam,
										    xh),
								    elevation_bw_r4_1(xKth,
								                      xgam,
										      xw),
										      xLn);
			    const float den   = 1984.4017075391884912304967f*xDx*xLt*xLa;
			    const float t1    = xgam*xgam;
			    const float t2    = Pav*xtf;
			    const float t3    = ag*ag;
			    const float t4    = rcs*xFrd*xFrd*xFp*xFp;
			    const float t5    = xF*xF*xF*xF*xFlen*xFlen;
			    const float num   = t1*t2*t3*t4*t5;
			    ratio             = num/den;
			    return (cephes_powf(ratio,0.25f));
		    }
                    

		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double thermal_noise_RR_r8_1(const RadarParamAoS_R8_1   &rp,
		                                  const JammerParamAoS_R8_1  &jp) {

		           /* const float dummy = rp.gamm;
			    const float rcs   = jp.sig;
		            const float dc    = duty_cycle_r4_1(rp.rho,
			                                      rp.tr);
                            const float Pav   = radar_avg_power_r4_1(rp.Pt,
			                                           dc);
			    const float ag    = radar_ant_gain_r4_1(azimuth_bw_r4_1(rp.Kth,
			                                                            rp.gamm,
										    rp.h),
								    elevation_bw_r4_1(rp.Kth,
								                      rp.gamm,
										      rp.w),
										      rp.Ln);
			    const float den   = 1984.4017075391884912304967f*rp.Dx*rp.Lt*rp.La;
			    const float t1    = rp.gamm*rp.gamm;
			    const float t2    = Pav*rp.tf;
			    const float t3    = ag*ag;
			    const float t4    = rcs*rp.Frdr*rp.Frdr*rp.Fp*rp.Fp;
			    const float t5    = rp.F*rp.F*rp.F*rp.F*rp.Flen*rp.Flen;
			    const float num   = t1*t2*t3*t4*t5;
			    return (cephes_powf(num/den,0.25f));*/
			    // More efficient implementation
			    const double rcs  = jp.sig;
			    const double xgam = rp.gamm;
			    const double xrho = rp.rho;
			    const double xtr  = rp.tr;
			    const double xPt  = rp.Pt;
			    const double xKth = rp.Kth;
			    const double xw   = rp.w;
			    const double xh   = rp.h;
			    const double xLn  = rp.Ln;
			    const double xDx  = rp.Dx;
			    const double xLt  = rp.Lt;
			    const double xLa  = rp.La;
			    const double xtf  = rp.tf;
			    const double xFrd = rp.Frdr;
			    const double xFp  = rp.Fp;
			    const double xF   = rp.F;
			    const double xFlen= rp.Flen;
			    double ratio      = 0.0;
			    double range      = 0.0;
			    const double dc    = duty_cycle_r8_1(xrho,
			                                        xtr);
                            const double Pav   = radar_avg_power_r8_1(xPt,
			                                             dc);
			    const double ag    = radar_ant_gain_r8_1(azimuth_bw_r8_1(xKth,
			                                                            xgam,
										    xh),
								    elevation_bw_r8_1(xKth,
								                      xgam,
										      xw),
										      xLn);
			    const double den   = 1984.4017075391884912304967*xDx*xLt*xLa;
			    const double t1    = xgam*xgam;
			    const double t2    = Pav*xtf;
			    const double t3    = ag*ag;
			    const double t4    = rcs*xFrd*xFrd*xFp*xFp;
			    const double t5    = xF*xF*xF*xF*xFlen*xFlen;
			    const double num   = t1*t2*t3*t4*t5;
			    ratio             = num/den;
			    return (std::pow(ratio,0.25));
		    }
                     

    }

}












#endif /*__GMS_RADAR_JAMMING_HPP__*/
