
#ifndef __GMS_RCS_HPP__
#define __GMS_RCS_HPP__ 230120221416

namespace file_info {

     const unsigned int GMS_RCS_MAJOR = 1;
     const unsigned int GMS_RCS_MINOR = 1;
     const unsigned int GMS_RCS_MICRO = 0;
     const unsigned int GMS_RCS_FULLVER =
       1000U*GMS_RCS_MAJOR+100U*GMS_RCS_MINOR+
       10U*GMS_RCS_MICRO;
     const char * const GMS_RCS_CREATION_DATE = "23-01-2022 14:16 +00200 (SUN 23 JAN 2022 GMT+2)";
     const char * const GMS_RCS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_RCS_SYNOPSIS      = "Calculation of various aspects of the RCS."

}

#include <cstdint>
#include <cmath> // for double precision version
#include "GMS_config.h"
#include "GMS_cephes.h" // for single precision version

namespace gms {


             namespace radiolocation {


	           namespace {
                        // Free-space target RCS [m^2]
                        const float   sig0_r4 = 1.0f;
                        const double  sig0_r8 = 1.0;
			const float   PIr4   = 3.14159265358979323846264338328f;
			const double  PIr8   = 3.14159265358979323846264338328;
			const float   PI2r4  = 6.283185307179586476925286766559f;
			const double  PI2r8  = 6.283185307179586476925286766559;
			const float   PI4r4  = 12.566370614359172953850573533118f;
			const double  PI4r8  = 12.566370614359172953850573533118;
                        const float   zr4    = 0.017453292519943295769236907685f;
			const double  zr8    = 0.017453292519943295769236907685; //coeff deg-to-rad conversion
		      
		   }

		   /*
                       This function computes the empirical K parameter when unknown apriori.
                       Forward-scatterer path when beta=180.
                    */

		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     float empirical_K_r4_1(const float sig0,  // m^2, RCS monostatic target
		                            const float A,     // m^2, Area of target projected to the normal radar beam
					    const float gamma) { // m,   wavelength
			   constexpr float PI24 = 0.74159265358979323846264338328f;
                           float K              = 0.0f;
			   const float A2       = A*A;
			   const float gamm2    = gamma*gamma*sig0;
			   const float t0       = PI4r4*(A2/gamm2);
			   K = cepehes_logf(t0)/PI24;
			   return (K);
		    }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     double empirical_K_r4_1(const double sig0,  // m^2, RCS monostatic target
		                             const double A,     // m^2, Area of target projected to the normal radar beam
					     const double gamma) { // m,   wavelength
			   constexpr double PI24 = 0.74159265358979323846264338328;
                           double K              = 0.0;
			   const double A2       = A*A;
			   const double gamm2    = gamma*gamma*sig0;
			   const double t0       = PI4r8*(A2/gamm2);
			   K = std::log(t0)/PI24;
			   return (K);
		    }

		                                              

		    /*
                        
                      */

	             __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void effective_rcs_r4_1(const float gamma, // m, wavelength
		                             const float R,     // nm, target range
					     const float h_a,   // ft, antenna height
					     const float h_t,   // ft, target height
					     float * __restrict sig_eff,    // m, RCS effective
					     float * __restrict sig_eff_db) { // db, RCS effective

                         const float t0 = PI2r4*h_a*h_t;
			 const float t1 = R*gamma;
			 const float F  = 2.0f*cephes_sin(t0/t1);
			 const float F4 = F*F*F*F;
			 *sig_eff       = sig0_r4*F4;
			 *sig_eff_db    = 10.0f*cephes_log10f(*sig_eff);
		   }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void effective_rcs_r8_1(const double gamma, // m, wavelength
		                             const double R,     // nm, target range
					     const double h_a,   // ft, antenna height
					     const double h_t,   // ft, target height
					     double * __restrict sig_eff,    // m, RCS effective
					     double * __restrict sig_eff_db) { // db, RCS effective

                         const double t0 = PI2r8*h_a*h_t;
			 const double t1 = R*gamma;
			 const double F  = 2.0*std::sin(t0/t1);
			 const double F4 = F*F*F*F;
			 *sig_eff        = sig0_r8*F4;
			 *sig_eff_db     = 10.0*std::log10(*sig_eff);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void
		    bistatic_target_rcs_r4_1(const float sig0, // m^2, RCS monostatic target
		                             const float K,    // empirical constant defined by target configuration and complexity  // use function emprirical_K for computation
					     const float Beta, // deg, bistatic angle (for input=1)
					     const float R1,   // m, transmitter - target range (input=2)
					     const float R2,   // m, receiver    - target range (input=2)
					     const float B,    // m, baseline
					     const int32_t type, // input switch (1,2)
					     float * __restrict sigma, // RCS                  
                                             float * __restrict sigma_db) {   // dBsm of Target
                         
                         if(type==1) {
                            const float alpha = Beta*zr4;
			    const float t0    = K*std::abs(alpha)-2.4f*K-1.0f;
			    const float t1    = 1.0f+cephes_expf(t0);
			    *sigma            = sig0*t1;
			    *sigma_db         = 10.0*cephes_log10f(*sigma+0.00000000001f);
			 }
			 else if(type==2) {
                            const float t0    = 1.0f/(2.0*R1*R2);
			    const float R12   = R1*R1;
			    const float R22   = R2*R2;
			    const float B2    = B*B;
			    const float t1    = R12+R22+B2;
			    const float gam   = cephes_acosf(t0*t1);
			    const float alpha = gam;
			    const float t2    = K*std::abs(alpha)-2.4f*K-1.0f;
			    const float t3    = 1.0f+cephes_expf(t0);
			    *sigma            = sig0*t1;
			    *sigma_db         = 10.0*cephes_log10f(*sigma+0.00000000001f);
			 }
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void
		    bistatic_target_rcs_r8_1(const double sig0, // m^2, RCS monostatic target
		                             const double K,    // empirical constant defined by target configuration and complexity  // use function emprirical_K for computation
					     const double Beta, // deg, bistatic angle (for input=1)
					     const double R1,   // m, transmitter - target range (input=2)
					     const double R2,   // m, receiver    - target range (input=2)
					     const double B,    // m, baseline
					     const int32_t type, // input switch (1,2)
					     double * __restrict sigma, // RCS                  
                                             double * __restrict sigma_db) {   // dBsm of Target
                         
                         if(type==1) {
                            const double alpha = Beta*zr8;
			    const double t0    = K*std::abs(alpha)-2.4*K-1.0;
			    const double t1    = 1.0+std::exp(t0);
			    *sigma            = sig0*t1;
			    *sigma_db         = 10.0*std::log10(*sigma+0.00000000001);
			 }
			 else if(type==2) {
                            const double t0    = 1.0/(2.0*R1*R2);
			    const double R12   = R1*R1;
			    const double R22   = R2*R2;
			    const double B2    = B*B;
			    const double t1    = R12+R22+B2;
			    const double gam   = std::acos(t0*t1);
			    const double alpha = gam;
			    const double t2    = K*std::abs(alpha)-2.4*K-1.0;
			    const double t3    = 1.0+std::exp(t0);
			    *sigma             = sig0*t1;
			    *sigma_db          = 10.0*std::log10(*sigma+0.00000000001);
			 }
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void antenna_rcs_r4_1(const float Ae,   // m^2, antenna efective aperture
		                          const float gam,  // m, wavelength
					  const float G,    // refelcetion coefficient
					  float * __restrict sigma, // m^2, radar cross section
					  float * __restrict sigma_db) { //dBsm, rcs with repsect to 1m^2

			  const float Ae2  = Ae*Ae;
			  const float gam2 = gam*gam;
			  const float t0   = Ae2/gam2;
			  *sigma           = PI4r4*t0*G;
			  *sigma_db        = 10.0f*cephes_log10f(*sigma);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
                    void antenna_rcs_r8_1(const double Ae,   // m^2, antenna efective aperture
		                          const double gam,  // m, wavelength
					  const double G,    // reflection coefficient
					  double * __restrict sigma, // m^2, radar cross section
					  double * __restrict sigma_db) { //dBsm, rcs with repsect to 1m^2

			  const double Ae2  = Ae*Ae;
			  const double gam2 = gam*gam;
			  const double t0   = Ae2/gam2;
			  *sigma            = PI4r8*t0*G;
			  *sigma_db         = 10.0*std::log10(*sigma);
		   }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    float bird_insect_rcs_r4_1(const float W) { //gram, weight of bird or insect

                           return (-46.0f+5.8f*cephes_log10f(W));
		    }


		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    static
		    inline
		    double bird_insect_rcs_r8_1(const double W) { //gram, weight of bird or insect

                           return (-46.0+5.8*std::log10(W));
		    }


		    


     }//radiolocation

} //gms











#endif /*__GMS_RCS_HPP__*/
