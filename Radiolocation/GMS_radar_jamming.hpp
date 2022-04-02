

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


	           /*
                          IEEE Standard 100, The Authoritative Dictionary of IEEE Standards Terms, 7th ed., New York:
                          IEEE Press, 2000.
                          [2] Barton, D. K., Radar System Analysis and Modeling, Norwood, MA: Artech House, 2005.
                          [3] Blake, L. V., Radar Range-Performance Analysis, Lexington, MA: D. C. Heath, 1980; Dedham,
                                            MA: Artech House, 1986.
                          [4] Barton, D. K. and H. R. Ward, Handbook of Radar Measurement, Englewood Cliffs, NJ: Pren-
                                            tice-Hall, 1968; Dedham, MA: Artech House, 1984.
                          [5] Goldstein, H., “Sea Echo,” Section 6.6, and “Meteorological Echoes,” Chapter 7, in Propaga-
                                            tion of Short Radio Waves, Vol. 13 in MIT Radiation Laboratory Series, New York: McGraw-
                                            Hill, 1951. Reprinted (CD ROM edition), Artech House, 1999.
                          [6] Doviak, R. J. and D. S. Zrnic, Doppler Radar and Weather Observations, 2nd ed., New York:
                                            Academic Press, 1993.
                          [7] Shrader, W. W. and V. Gregers-Hansen, “MTI Radar,” Chap. 2 in Radar Handbook, 3rd ed.,
                                            (M. I. Skolnik, ed.), New York: McGraw-Hill, 2008.
                          [8] Schleher, D. C., Introduction to Electronic Warfare, Dedham, MA: Artech House, 1986.
                          [9] Schleher, D. C., Electronic Warfare in the Information Age, Norwood, MA: Artech House,
                     */
                 
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


		     // Prefetch the data to L1D.
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
                     bool prefetch_data_r4_1(const RadarParamAoS_R4_1   &rp,
		                             const JammerParamAoS_R4_1  &jp) {

                           _mm_prefetch((const char*)&rp,0);
			   _mm_prefetch((const char*)&jp,0);
			   return (true); // can not fail at software visible way.
		     }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
                     bool prefetch_data_r8_1(const RadarParamAoS_R8_1   &rp,
		                             const JammerParamAoS_R8_1  &jp) {

                           _mm_prefetch((const char*)&rp,0);
			   _mm_prefetch((const char*)&jp,0);
			   return (true); // can not fail at software visible way.
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

                
		    
	           // Initialize 'RadarParamsAoS_R4_1 data type
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void initRadarParamsAoS_R4_1(const float xgamm,
		                                  const float xtf,
						  const float xrho,
						  const float xw,
						  const float xKth,
						  const float xLn,
						  const float xTs,
						  const float xFp,
						  const float xLa,
						  const float xF,
						  const float xPt,
						  const float xtr,
						  const float xLt,
						  const float xh,
						  const float xha,
						  const float xFrdr,
						  const float xDx,
						  const float xBt,
						  const float xFlen,
		                                  RadarParamsAoS_R4_1 &rp) {

                            rp.gamm = xgamm;
			    rp.xtf  = xtf;
			    rp.rho  = xrho;
			    rp.w    = xw;
			    rp.Kth  = xKth;
			    rp.Ln   = xLn;
			    rp.Ts   = xTs;
			    rp.Fp   = xFp;
			    rp.La   = xLa;
			    rp.F    = xF;
			    rp.Pt   = xPt;
			    rp.tr   = xtr;
			    rp.Lt   = xLt;
			    rp.h    = xh;
			    rp.ha   = xha;
			    rp.Frdr = xFrdr;
			    rp.Dx   = xDx;
			    rp.xBt  = xBt;
			    rp.xFlen= xFlen;
		     }


		      // Initialize 'RadarParamsAoS_R8_1 data type
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void initRadarParamsAoS_R8_1(const double xgamm,
		                                  const double xtf,
						  const double xrho,
						  const double xw,
						  const double xKth,
						  const double xLn,
						  const double xTs,
						  const double xFp,
						  const double xLa,
						  const double xF,
						  const double xPt,
						  const double xtr,
						  const double xLt,
						  const double xh,
						  const double xha,
						  const double xFrdr,
						  const double xDx,
						  const double xBt,
						  const double xFlen,
		                                  RadarParamsAoS_R8_1 &rp) {

                            rp.gamm = xgamm;
			    rp.xtf  = xtf;
			    rp.rho  = xrho;
			    rp.w    = xw;
			    rp.Kth  = xKth;
			    rp.Ln   = xLn;
			    rp.Ts   = xTs;
			    rp.Fp   = xFp;
			    rp.La   = xLa;
			    rp.F    = xF;
			    rp.Pt   = xPt;
			    rp.tr   = xtr;
			    rp.Lt   = xLt;
			    rp.h    = xh;
			    rp.ha   = xha;
			    rp.Frdr = xFrdr;
			    rp.Dx   = xDx;
			    rp.xBt  = xBt;
			    rp.xFlen= xFlen;
		     }


		     // Initialize JammerParamAoS_R4_1 data type
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void initJammerParamsAoS_R4_1(const float xsig,
		                                   const float xPj,
						   const float xGj,
						   const float xQj,
						   const float xFlenj,
						   const float xRj,
						   const float xBj,
						   const float xLtj,
						   const float xFpj,
						   const float xRmj,
						   const float xFj,
						   const float xLaj,
						   JammerParamsAoS_R4_1 &jp) {

                              jp.sig  = xsig;
			      jp.Pj   = xPj;
			      jp.Gj   = xGj;
			      jp.Qj   = xQj;
			      jp.Flenj= xFlenj;
			      jp.Rj   = xRj;
			      jp.Bj   = xBj;
			      jp.Ltj  = xLtj;
			      jp.Fpj  = xFpj;
			      jp.Rmj  = xRmj;
			      jp.Fj   = xFj;
			      jp.Laj  = xLaj;
		    }


		    // Initialize JammerParamAoS_R8_1 data type
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void initJammerParamsAoS_R8_1(const double xsig,
		                                   const double xPj,
						   const double xGj,
						   const double xQj,
						   const double xFlenj,
						   const double xRj,
						   const double xBj,
						   const double xLtj,
						   const double xFpj,
						   const double xRmj,
						   const double xFj,
						   const double xLaj,
						   JammerParamsAoS_R8_1 &jp) {

                              jp.sig  = xsig;
			      jp.Pj   = xPj;
			      jp.Gj   = xGj;
			      jp.Qj   = xQj;
			      jp.Flenj= xFlenj;
			      jp.Rj   = xRj;
			      jp.Bj   = xBj;
			      jp.Ltj  = xLtj;
			      jp.Fpj  = xFpj;
			      jp.Rmj  = xRmj;
			      jp.Fj   = xFj;
			      jp.Laj  = xLaj;
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


	             // Barrage Noise Jamming equations

		     //Troposhperic power loss at specific range.
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float tropo_range_loss_r4_1(const RadarParamAoS_R4_1   &rp,
		                                 const JammerParamsAoS_R4_1 &jp) {

                             const float xLa = rp.La;
			     const float xRmj= jp.Rmj;
			     const float Rm  =  thermal_noise_RR_r4_1(rp,jp);
			     return (xLa*(xRmj/Rm));
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double tropo_range_loss_r8_1(const RadarParamAoS_R8_1   &rp,
		                                  const JammerParamsAoS_R8_1 &jp) {
                             _mm_prefetch((const char*)&rp,0);
			     _mm_prefetch((const char*)&jp,0);
                             const double xLa = rp.La;
			     const double xRmj= jp.Rmj;
			     const double Rm  =  thermal_noise_RR_r8_1(rp,jp);
			     return (xLa*(xRmj/Rm));
		    }


		    // Effective radiated power of Jammer (W)
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float jammer_erp_r4_1(const JammerParamAoS_R4_1 &jp) {

                            const float xPj = jp.Pj;
			    const float xGj = jp.Gj;
			    const float xLtj= jp.Ltj;
			    return ((xPj*xGj)/xLtj);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double jammer_erp_r8_1(const JammerParamAoS_R8_1 &jp) {

                            const double xPj = jp.Pj;
			    const double xGj = jp.Gj;
			    const double xLtj= jp.Ltj;
			    return ((xPj*xGj)/xLtj);
		    }


		    // Effective radiated Jammer noise power (W)
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float jammer_ernp_r4_1(const JammerParamAoS_R4_1 &jp) {

                            const float xQj = jp.Qj;
			    const float xPj = jp.Pj;
			    const float xGj = jp.Gj;
			    const float xFpj= jp.Fpj;
			    const float xLtj= jp.Ltj;
			    const float xFpj2 = xFpj*xFpj;
			    return ((xQj*xPj*xGj*xFpj2)/xLtj);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double jammer_ernp_r8_1(const JammerParamAoS_R8_1 &jp) {

                            const double xQj = jp.Qj;
			    const double xPj = jp.Pj;
			    const double xGj = jp.Gj;
			    const double xFpj= jp.Fpj;
			    const double xLtj= jp.Ltj;
			    const double xFpj2 = xFpj*xFpj;
			    return ((xQj*xPj*xGj*xFpj2)/xLtj);
		    }


		    // Jamming spectral density (W/Hz)
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float jammer_sd_r4_1(const JammerParamsAoS_R4_1 &jp,
		                          const RadarParamsAoS_R4_1  &rp) {

			  // _mm_prefetch((const char*)&rp,0);
			   //_mm_prefetch((const char*)&jp,0);
			   const float j0    = 0.0f;
			   const float xKth  = rp.Kth; //1st cache miss for jp load
			   const float xgam  = rp.gamm;
			   const float xh    = rp.h;
			   const float xw    = rp.w;
			   const float xLn   = rp.Ln;
			   const float g2    = rp.gamm*rp.gamm;
			   const float xGr   = radar_ant_gain_r4_1(azimuth_bw_r4_1(xKth,xgam,xh),
			                                   elevation_bw_r4_1(xKth,xgam,xw),xLn);
			   const float xQj   = jp.Qj;// 1st cache miss for rp load	                     
                           const float xPj   = jp.Pj;
			   const float xGj   = jp.Gj;
			   const float xFpj  = jp.Fpj;
			   const float xFlnsj= jp.Flenj;
			   const float xFj   = jp.Fj;
			   const float xRj   = jp.Rj;
			   const float xBj   = jp.Bj;
			   const float xLtj  = jp.Ltj;
			   const float xLaj  = jp.Laj;
			   const float PI42  = 157.9136704174297379013522f;
			   const float r0    = xRj*xBj*xLtj*xLaj;
			   const float t0    = xQj*xPj*xGj*xGr;
			   const float t1    = g2*xFpj*xFpj*xFlnsj*xFlnsj;
			   const float t2    = t0*t1*xFj*xFj;
			   j0                = t2/(PI42*r0);
			   return (j0);
		   }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double jammer_sd_r8_1(const JammerParamsAoS_R8_1 &jp,
		                           const RadarParamsAoS_R8_1  &rp) {

			   //_mm_prefetch((const char*)&rp,0);
			   //_mm_prefetch((const char*)&jp,0);
			   const double j0    = 0.0;
			   const double xKth  = rp.Kth; //1st cache miss for jp load
			   const double xgam  = rp.gamm;
			   const double xh    = rp.h;
			   const double xw    = rp.w;
			   const double xLn   = rp.Ln;
			   const double g2    = rp.gamm*rp.gamm;
			   const double xGr   = radar_ant_gain_r8_1(azimuth_bw_r8_1(xKth,xgam,xh),
			                                   elevation_bw_r8_1(xKth,xgam,xw),xLn);
			   const double xQj   = jp.Qj;// 1st cache miss for rp load	                     
                           const double xPj   = jp.Pj;
			   const double xGj   = jp.Gj;
			   const double xFpj  = jp.Fpj;
			   const double xFlnsj= jp.Flenj;
			   const double xFj   = jp.Fj;
			   const double xRj   = jp.Rj;
			   const double xBj   = jp.Bj;
			   const double xLtj  = jp.Ltj;
			   const double xLaj  = jp.Laj;
			   const double PI42  = 157.9136704174297379013522;
			   const double r0    = xRj*xBj*xLtj*xLaj;
			   const double t0    = xQj*xPj*xGj*xGr;
			   const double t1    = g2*xFpj*xFpj*xFlnsj*xFlnsj;
			   const double t2    = t0*t1*xFj*xFj;
			   j0                 = t2/(PI42*r0);
			   return (j0);
		   }


		   // Available jamming temperature single jammer scenario.
		   // Remark: the implementation is identical to function coded above.
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float jammer_t1_r4_1(const JammerParamsAoS_R4_1 &jp,
		                          const RadarParamsAoS_R4_1  &rp) {

			   //_mm_prefetch((const char*)&rp,0);
			   //_mm_prefetch((const char*)&jp,0);
			   const float jt1    = 0.0f;
			   const float xKth  = rp.Kth; //1st cache miss for jp load
			   const float xgam  = rp.gamm;
			   const float xh    = rp.h;
			   const float xw    = rp.w;
			   const float xLn   = rp.Ln;
			   const float g2    = rp.gamm*rp.gamm;
			   const float xGr   = radar_ant_gain_r4_1(azimuth_bw_r4_1(xKth,xgam,xh),
			                                   elevation_bw_r4_1(xKth,xgam,xw),xLn);
			   const float xQj   = jp.Qj;// 1st cache miss for rp load	                     
                           const float xPj   = jp.Pj;
			   const float xGj   = jp.Gj;
			   const float xFpj  = jp.Fpj;
			   const float xFlnsj= jp.Flenj;
			   const float xFj   = jp.Fj;
			   const float xRj   = jp.Rj;
			   const float xBj   = jp.Bj;
			   const float xLtj  = jp.Ltj;
			   const float xLaj  = jp.Laj;
			   const float PI42  = 157.9136704174297379013522f;
			   const float r0    = xRj*k_B4*xBj*xLtj*xLaj;
			   const float t0    = xQj*xPj*xGj*xGr;
			   const float t1    = g2*xFpj*xFpj*xFlnsj*xFlnsj;
			   const float t2    = t0*t1*xFj*xFj;
			   jt1               = t2/(PI42*r0);
			   return (jt1);
		   }

		     
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double jammer_t1_r8_1(const JammerParamsAoS_R8_1 &jp,
		                           const RadarParamsAoS_R8_1  &rp) {

			   //_mm_prefetch((const char*)&rp,0);
			   //_mm_prefetch((const char*)&jp,0);
			   const double jt1    = 0.0;
			   const double xKth  = rp.Kth; //1st cache miss for jp load
			   const double xgam  = rp.gamm;
			   const double xh    = rp.h;
			   const double xw    = rp.w;
			   const double xLn   = rp.Ln;
			   const double g2    = rp.gamm*rp.gamm;
			   const double xGr   = radar_ant_gain_r8_1(azimuth_bw_r8_1(xKth,xgam,xh),
			                                   elevation_bw_r8_1(xKth,xgam,xw),xLn);
			   const double xQj   = jp.Qj;// 1st cache miss for rp load	                     
                           const double xPj   = jp.Pj;
			   const double xGj   = jp.Gj;
			   const double xFpj  = jp.Fpj;
			   const double xFlnsj= jp.Flenj;
			   const double xFj   = jp.Fj;
			   const double xRj   = jp.Rj;
			   const double xBj   = jp.Bj;
			   const double xLtj  = jp.Ltj;
			   const double xLaj  = jp.Laj;
			   const double PI42  = 157.9136704174297379013522;
			   const double r0    = xRj*k_B8*xBj*xLtj*xLaj;
			   const double t0    = xQj*xPj*xGj*xGr;
			   const double t1    = g2*xFpj*xFpj*xFlnsj*xFlnsj;
			   const double t2    = t0*t1*xFj*xFj;
			   jt1                = t2/(PI42*r0);
			   return (jt1);
		   }


		   // Required signal jamming temperature.
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float jammer_treq_r4_1(const JammerParamsAoS_R4_1 &jp,
		                            const RadarParamsAoS_R4_1  &rp) {

                            //_mm_prefetch((const char*)&rp,0);
			    //_mm_prefetch((const char*)&jp,0);
			    float tj         = 0.0f;
			    const float xTs  = rp.Ts;
			    const float Rm   = thermal_noise_RR_r4_1(rp,jp);
			    const float xRmj = jp.Rmj;
			    const float xLa  = rp.La;
			    const float xLa1 = tropo_range_loss_r4_1(rp,jp);
			    const float xFlns= rp.Flen;
			    const float xFln1= cephes_sqrtf(xFlns);
			    const float lrat = xLa/xLa1;
			    const float mrat = xFlns/xFln1;
			    const float rrat = Rm/xRmj;
			    const float rrat4= (rrat*rrat*rrat*rrat)-1.0f;
			    tj               = xTs*lrat*(mrat*mrat)*rrat4;
			    return (tj);
		   }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double jammer_treq_r8_1(const JammerParamsAoS_R8_1 &jp,
		                             const RadarParamsAoS_R8_1  &rp) {

                           // _mm_prefetch((const char*)&rp,0);
			   // _mm_prefetch((const char*)&jp,0);
			    double tj         = 0.0;
			    const double xTs  = rp.Ts;
			    const double Rm   = thermal_noise_RR_r8_1(rp,jp);
			    const double xRmj = jp.Rmj;
			    const double xLa  = rp.La;
			    const double xLa1 = tropo_range_loss_r8_1(rp,jp);
			    const double xFlns= rp.Flen;
			    const double xFln1= std::pow(xFlns);
			    const double lrat = xLa/xLa1;
			    const double mrat = xFlns/xFln1;
			    const double rrat = Rm/xRmj;
			    const double rrat4= (rrat*rrat*rrat*rrat)-1.0;
			    tj               = xTs*lrat*(mrat*mrat)*rrat4;
			    return (tj);
		   }


	             __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float n_jammers_req_r4_1(const JammerParamAoS_R4_1 &jp,
		                              const RadarParamAoS_R4_1  &rp) {

                            const float xTj  = jammer_treq_r4_1(jp,rp);
			    const float xTj1 = jammer_t1_r4_1(rp,jp);
			    return (std::ceilf(xTj/xTj1));
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double n_jammers_req_r8_1(const JammerParamAoS_R8_1 &jp,
		                              const RadarParamAoS_R8_1  &rp) {

                            const double xTj  = jammer_treq_r8_1(jp,rp);
			    const double xTj1 = jammer_t1_r8_1(rp,jp);
			    return (std::ceil(xTj/xTj1));
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float n_jammers_margin_r4_1(const JammerParamAoS_R4_1 &jp,
		                                 const RadarParamAoS_R4_1  &rp,
						 const float nj) { // number of jammers required

                            const float xTj  = jammer_treq_r4_1(jp,rp);
			    const float xTj1 = jammer_t1_r4_1(rp,jp);
			    return ((nj*xTj1)/xTj);
		   }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double n_jammers_margin_r8_1(const JammerParamAoS_R8_1 &jp,
		                                  const RadarParamAoS_R8_1  &rp,
						  const double nj) { // number of jammers required

                            const double xTj  = jammer_treq_r8_1(jp,rp);
			    const double xTj1 = jammer_t1_r8_1(rp,jp);
			    return ((nj*xTj1)/xTj);
		   }


		   // Number of jammers range (km)
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float n_jammers_range_r4_1(const JammerParamAoS_R4_1 &jp,
		                                const RadarParamAoS_R4_1  &rp) {

			    const float xKth = rp.Kth;
			    const float xgam = rp.gamm*rp.gamm;
			    const float xh   = rp.h;
			    const float xw   = rp.w;
			    const float xLn  = rp.Ln;
			    const float xrho = rp.rho;
			    const float xtf  = rp.tf;
			    const float xtr  = rp.tr;
			    const float xPt  = rp.Pt;
                            const float Gr   = radar_ant_gain_r4_1(azimuth_bw_r4_1(xKth,xgam,xh),
			                                   elevation_bw_r4_1(xKth,xgam,xw),xLn);
			    const float Gr2  = Gr*Gr;
			    const float xFrd = rp.Frdr;
			    const float xFp  = rp.Fp*rp.Fp;
			    const float xF   = rp.F;
			    const float xFlen= rp.Flen;
			    const float xsig = jp.sig;
			    const float xTs  = rp.Ts;
			    const float xDx  = rp.Dx;
			    const float xLt  = rp.Lt;
			    const float xLa  = rp.La;
			    const float dc   = duty_cycle_r4_1(xrho,xtr);
			    const float pav  = radar_avg_power_r4_1(xPt,dc);                           
                            const float xN0  = noise_density_r4_1(xTs);
			    const float xF4  = xF*xF*xF*xF;
			    const float num  = pav*xtf*Gr2*xgamm*xsig*xF4*xFp;
			    const float den  = 1984.4017075391884912304842f*k_B4*xTs*xDx*xLt*xLa;
			    const float ratio= num/den;
			    return (cephes_powf(ratio,0.25f));             
		   }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double n_jammers_range_r8_1(const JammerParamAoS_R8_1 &jp,
		                                 const RadarParamAoS_R8_1  &rp) {

			    const double xKth = rp.Kth;
			    const double xgam = rp.gamm*rp.gamm;
			    const double xh   = rp.h;
			    const double xw   = rp.w;
			    const double xLn  = rp.Ln;
			    const double xrho = rp.rho;
			    const double xtf  = rp.tf;
			    const double xtr  = rp.tr;
			    const double xPt  = rp.Pt;
                            const double Gr   = radar_ant_gain_r8_1(azimuth_bw_r8_1(xKth,xgam,xh),
			                                   elevation_bw_r8_1(xKth,xgam,xw),xLn);
			    const double Gr2  = Gr*Gr;
			    const double xFrd = rp.Frdr;
			    const double xFp  = rp.Fp*rp.Fp;
			    const double xF   = rp.F;
			    const double xFlen= rp.Flen;
			    const double xsig = jp.sig;
			    const double xTs  = rp.Ts;
			    const double xDx  = rp.Dx;
			    const double xLt  = rp.Lt;
			    const double xLa  = rp.La;
			    const double dc   = duty_cycle_r8_1(xrho,xtr);
			    const double pav  = radar_avg_power_r8_1(xPt,dc);                           
                            const double xN0  = noise_density_r8_1(xTs);
			    const double xF4  = xF*xF*xF*xF;
			    const double num  = pav*xtf*Gr2*xgamm*xsig*xF4*xFp;
			    const double den  = 1984.4017075391884912304842*k_B8*xTs*xDx*xLt*xLa;
			    const double ratio= num_den;
			    return (std::pow(ratio,0.25));             
		   }


    }

}












#endif /*__GMS_RADAR_JAMMING_HPP__*/
