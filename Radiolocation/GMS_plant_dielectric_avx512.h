
#ifndef __GMS_PLANT_DIELECTRIC_AVX512_H__
#define __GMS_PLANT_DIELECTRIC_AVX512_H__ 050120220920



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


#include "GMS_avx512vecf32.h"
#include "GMS_avx512vecf64.h"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_complex_zmm8r8.hpp"
#include "GMS_config.h"


namespace gms {


         namespace radiolocation {





                         
	              /*
                               This kernel operates on 16 leaves.
                          */
	            
                      __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __ATTR_ALIGN__(32)
		      ZMM16c4
		      zmm16c4_leaf_dielectric(const AVX512Vec16 leaf_mg,
		                              const AVX512Vec16 leaf_rho,
					      const AVX512Vec16 leaf_dens,
					      const AVX512Vec16 leaf_tau,
					      const AVX512Vec16 water_tmp,
					      const AVX512Vec16 veg_tmp,
					      const AVX512Vec16 theta,
					      const AVX512Vec16 freq,  // frequency value is a scalar broadcast to ZMM register
					      const bool dry_dens); 


		     
                      __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __ATTR_ALIGN__(32)
		     ZMM8c8
		     zmm8c8_leaf_dielectric( const AVX512Vec8 leaf_mg,
		                              const AVX512Vec8 leaf_rho,
					      const AVX512Vec8 leaf_dens,
					      const AVX512Vec8 leaf_tau,
					      const AVX512Vec8 water_tmp,
					      const AVX512Vec8 veg_tmp,
					      const AVX512Vec8 theta,
					      const AVX512Vec8 freq,  // frequency value is a scalar broadcast to ZMM register
					      const bool dry_dens); 
					      
					      
		      __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __ATTR_ALIGN__(32)		      
		      ZMM16c4
		      zmm16c4_veg_dielectric_2(const AVX512Vec16 mg,
		                               const AVX512Vec16 veg_rho,
					       const AVX512Vec16 tempC,
					       const AVX512Vec16 theta,
					       const AVX512Vec16 freq); 

		      
                      __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __ATTR_ALIGN__(32)
		      ZMM8c8
		      zmm8c8_veg_dielectric_2( const AVX512Vec8 mg,
		                               const AVX512Vec8 veg_rho,
					       const AVX512Vec8 tempC,
					       const AVX512Vec8 theta,
					       const AVX512Vec8 freq); 


		      
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __ATTR_ALIGN__(32)
		      ZMM16c4
		     zmm16c4_veg_dielectric_1(const AVX512Vec16 mg,
		                               const AVX512Vec16 tempC,
					       const AVX512Vec16 theta,
					       const AVX512Vec16 freq); 

		      
                       __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __ATTR_ALIGN__(32)
		     ZMM16c4
		     zmm8c8_veg_dielectric_1( const AVX512Vec8 mg,
		                               const AVX512Vec8 tempC,
					       const AVX512Vec8 theta,
					       const AVX512Vec8 freq); 

		   
		   

		   

		  


		  













		    
     } // radiolocation


} // gms














#endif /*__GMS_PLANT_DIELECTRIC_AVX512_H__*/
