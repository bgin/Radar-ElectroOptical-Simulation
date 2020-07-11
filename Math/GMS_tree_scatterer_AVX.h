

#ifndef __GMS_TREE_SCATTERER_AVX_H__
#define __GMS_TREE_SCATTERER_AVX_H__


namespace file_info {

      const unsigned int gGMS_TREE_SCATTERER_AVX_MAJOR = 1;
      const unsigned int gGMS_TREE_SCATTERER_AVX_MINOR = 0;
      const unsigned int gGMS_TREE_SCATTERER_AVX_MICRO = 0;
      const unsigned int gGMS_TREE_SCATTERER_AVX_FULLVER = 1000U*gGMS_TREE_SCATTERER_AVX_MAJOR+100U*gGMS_TREE_SCATTERER_AVX_MINOR+
                                                       10U*gGMS_TREE_SCATTERER_AVX_MICRO;
      const char * const pgGMS_TREE_SCATTERER_AVX_CREATION_DATE = "24-12-2019 13:52 +00200 (TUE 24 DEC 2019 GMT+2)";
      const char * const pgGMS_TREE_SCATTERER_AVX_BUILD_DATE    = __DATE__ " " __TIME__;
      const char * const pgGMS_TREE_SCATTERER_AVX_AUTHOR        = " Programmer: Bernard Gingold, contact: beniekg@gmail.com";
      const char * const pgGMS_TREE_SCATTERER_AVX_SYNOPSYS      = " Model of single tree scatterer(AVX implementation) suitable for computation of radar backscatter.";

}

#include <cstdint>
#include <complex>
#include "GMS_config.h"
#include "GMS_avxvecf32.h"



namespace  gms {

     namespace  math {

             namespace {

	           const AVXVec8 VINC0  =   AVXVec8{1.0f,2.0f,3.0f,4.0f,
		                                  5.0f,6.0f,7.0f,8.0f};
		   const AVXVec8 VINC1  =   AVXVec8{9.0f,10.0f,11.0f,12.0f,
		                                  13.0f,14.0f,15.0f,16.0f};
		   const AVXVec8 VINC2  =   AVXVec8{17.0f,18.0f,19.0f,20.0f,
		                                   21.0f,22.0f,23.0f,24.0f};
		   const AVXVec8 VINC3  =   AVXVec8{25.0f,26.0f,27.0f,28.0f,
		                                   29.0f,30.0f,31.0f,32.0f};
                   const AVXVec8 ZERO   =   AVXVec8{};

		   const AVXVec8 TWO    =   AVXVec8{2.0f};

		   const AVXVec8 PI     =   AVXVec8{3.141592653589793f};
	     }

	     // Low temporal access and spatial locality (cold)  data structure (POD) type 

	             struct TSColdAVX_t {

                         // Number of leaves (represented as an ellipsoidal surface)
			 int32_t nleaves;
			 // Number of branches (represented as an cylindrical volumes)
			 int32_t nbranches;
			 //  Number of simulation steps it is equal to Radar PRF (pulse repetetive frequency)
			 int32_t nsteps;
			 //  Tree scatterer ordinal number (for the forest simulation)
			 int32_t ordinal;
			 // Number of parametric equation evaluation 'points' for the trunk cylindrical approximation
                         int32_t trunk_param_npts;
			 //  Number of parametric equation evaluation 'points' for the leaves elliptical approximation
			 int32_t leaves_param_npts;
			 //   Number of parametric equation evaluation 'points' for the branches cylindrical approximation
			 int32_t branches_param_npts;
			 //  Total height of tree
			 float   tree_height;
			 //   Height of the trunk only
			 float   trunk_height;
			 //   Radius of trunk  (averaged)
			 float   trunk_radius;
			 //    Height of crown only
			 float   crown_height;
			 //    Total crown area (approximated) as sum of leaves area
			 float   crown_area;
			 //    Trunk area (cylinder area)
			 float   trunk_area;
			 //    Total tree area
			 float   tree_area;
			 //     Tree geo-location latitude
			 float   tree_lat;
			 //     Tree geo-location longtitude
			 float   tree_lon; // 1st cache line ends here!
			 //     Tree elevation (above the sea level) (meters)
			 float   tree_elevation; // 2nd cache lines starts here
#if (USE_STRUCT_PADDING) == 1
                         PAD_TO_ALIGNED(4,0,4)
#endif
                         // Is water or rather moistness present or not on the branch surface (per n-leaves) (allowed values only 0,1)
			 int32_t * __restrict __ATTR_ALIGN__(8) leaves_moist;
			 // Is water or rather moistness present or not on the branch surface (per n-branches) (allowed values only 0,1)
			 int32_t * __restrict __ATTR_ALIGN__(8) branches_moist;
			 //  ! Trunk parametric equation (approximated as a cylindrical object)
                         //   PAOS type size of array -- npoints/8
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) trunk_xparam;
			 //   ! PAOS type size of array -- npoints/8
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) trunk_yparam;
			 //    ! PAOS type size of array -- npoints/8
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) trunk_zparam;
			 //    ! Leaves thicknes (micron) per leaf
                         //    ! PAOS type size of array nleaves/8
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) leaves_thick;
			 //    ! Leaves density (g/cm^3) per leaf
                         //    ! PAOS type size of array nleaves/8
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) leaves_dens;
		
			 //  ! Leaves parameteric equation (approximated as an ellipses)
                         //  ! Parameter x,(a*cos(t))
        
                         //  ! PAOS type size of arrays  1st dim (evaluation of x) ,
                         //  !  2nd dim (number of leaves)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) leaves_xparam;
			 //   ! Parameter y, (b*sin(t))
                         //   ! PAOS type size of arrays is -- npoints/8 1st dim (evaluation of y) ,
                         //   ! nleaves/8 2nd dim (number of leaves)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) leaves_yparam;
			 //   Branchess thickness
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) branches_thick;
			 //   Branches density
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) branches_dens;
		
			 //  ! Branches parametric equation (approximated as a cylindrical objects)
                         //  ! Parameter x, (r*cos(t))
                         //  ! PAOS type size of arrays is -- npoints/8 1st dim (evaluation of x) ,
                         //  ! nbranches/8 2nd dim (number of leaves)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) branches_xparam;
			 //  ! Parameter y, (r*sin(t))
                         //  !2nd dimension is a branches number, 1st dimension evaluation of parameter y
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) branches_yparam;
			 //  ! Parameter z, (height)
                         //  !2nd dimension is a branch  number, 1st dimension evaluation of parameter z
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) branches_zparam;
		} __ATTR_ALIGN__(64);

		//  ! This is a high termporal and spatial locality data type
                //  ! These data type members characteristics are varying between each sample of Radar PRF.

		struct TSHotAVX_t {

                         // Whole tree direction theta angle (rad)
			 float tree_dtheta;
			 // Whole tree direction phi angle (rad)
			 float three_dphi;
			 //
			 float tree_rcs;
			 //  Crown cross section approximated as sum of leaves cross section.
			 float crown_rcs;
			 // Trunk cross section
			 float trunk_rcs;
#if (USE_STRUCT_PADDING) == 1
                         PAD_TO_ALIGNED(4,0,28)
#endif
                         //  ! Leaves cross section (varying due to leaves vibrations)
                         //  ! ( 1st dimension cross section variation per leaf,2nd dimension is PRF/s )
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) leaves_rcs;
			 //  ! Leaves reflectivity (varying due to leaves vibration)
                         //  ! ( 1st dimension reflectivity(dbm) pers leaf, 2nd dimension is PRF/s,)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) leaves_reflect;
			 //  ! Branches cross section (varying due to branches vibrations)
                         // ! ( 1st dimension cross section variation per branch, 2nd dimension is PRF/s))
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) branches_rcs;
			 //   ! Branches reflectivity (varying due to leaves vibration)
                         //     ! ( 1st dimension reflectivity(dbm) pers branch, 2nd dimension is PRF/s))
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) branches_reflect;
		
		} __ATTR_ALIGN__(64);

                // This is high spatial and temporal data structure.
		struct LeavesPhase_t {
                          // Theta angle of incoming EM-wave (rad) per each leaf
		         float * __restrict __ATTR_ALIGN__(8)  theta_inc; 
                          // Phi angle of incoming EM-wave (rad) per each leaf
			 float * __restrict __ATTR_ALIGN__(8)  phi_inc;
			  // THeta angle of scattered EM-wave (rad) per each leaf
			 float * __restrict __ATTR_ALIGN__(8)  theta_scat
			 // Phi angle of scattered EM-wave (rad) per each leaf
			 float * __restrict __ATTR_ALIGN__(8)  phi_scat;
			 // Theta angle of leaf orientation (rad) per each leaf
			 float * __restrict __ATTR_ALIGN__(8)  theta_dir;
			 // Phi angle of leaf oerientation (rad) per each leaf
			 float * __restrict __ATTR_ALIGN__(8)  phi_dir;
			 // Allocate this array as a [4*4*4*nleaves]
			 float * __restrict __ATTR_ALIGN__(8)  l4x4phm;
			 // Alocate this array as a [2*2*2*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) sm2x2avg;
			 // Allocate this array as a [2*2*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) l2x2mp;
			 // Allocate this array as a [2*2*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) l2x2mn;
			 // Allocate this array as a [4*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) eig1x4lp;
			 // Allocate this array as a [4*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) eig1x4ln;
			 // Allocate this array as a [4*4*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) eig4x4mp;
			 // Allocate this array as a [4*4*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) eig4x4mn;
			 // Allocate this array as a [4*4*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) eig4x4mpi;
			 // Allocate this array as a [4*4*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) eig4x4mni;
			 // Allocate this array as a [4*4*nleaves]
			 float * __restrict __ATTR_ALIGN__(8) expa4x4mp;
			 // Allocate this array as a [4*4*nleaves]
			 float * __restrict __ATTR_ALIGN__(8) expa4x4mn;
			 // Allocate this array as a [4*4*nleaves]
			 float * __restrict __ATTR_ALIGN__(8) stokes4x4m;
			 // Allocate this array as a [2*2*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) scat2x2m;
		}  _ATTR_ALIGN__(64);


		struct TreeScattererAVX {


		       TSColdAVX_t m_tsc __ATTR_ALIGN__(64);

		       TSHotAVX_t  m_tsh __ATTR_ALIGN__(64);

		       LeavesPhase_t m_lp __ATTR_ALIGN__(64);


		       TreeScattererAVX() __ATTR_COLD__ __ATTR_ALIGN__(32);

		       TreeScattererAVX(const int32_t,
		                     const int32_t,
				     const int32_t,
				     const int32_t,
				     const int32_t,
				     const int32_t,
				     const int32_t,
				     const float,
				     const float,
				     const float,
				     const float,
				     const float,
				     const float,
				     const float)  __ATTR_COLD__ __ATTR_ALIGN__(32);

			TreeScattererAVX(const TreeScattererAVX &) = delete;

			TreeScattererAVX(TreeScattererAVX &&) = delete;

			~TreeScattererAVX() noexcept(true);

			TreeScattererAVX & operator=(const TreeScattererAVX &) = delete;

			TreeScattererAVX & operator=(TreeScattererAVX &&) = delete;

			void SetMoistnessMask() __ATTR_COLD__ __ATTR_ALIGN__(32);

			void ComputeTrunkParamEq_ymm8r4(const int32_t) __ATTR_COLD__ __ATTR_ALIGN__(32);

			void SetThickDensAng_ymm8r4(const AVXVec8 * __restrict __ATTR_ALIGN__(64)) __ATTR_COLD__ __ATTR_ALIGN__(32);

			void ComputeLeavesParamEq_ymm8r4( const AVXVec8, const AVXVec8 ) __ATTR_COLD__ __ATTR_ALIGN__(32);

			void ComputeBranchesParamEq_ymm8r4(const int32_t ) __ATTR_COLD__ __ATTR_ALIGN__(32);


			__ATTR_COLD__
			__ATTR_ALIGN__(16)
			inline
			void Set_leaf_quadrature_bounds(int32_t & nth1,
			                                float   & tr_start1,
							float   & tr_stop1,
							float   & dt_rad1,
							int32_t & nth2,
			                                float   & tr_start2,
							float   & tr_stop2,
							float   & dt_rad2,
							int32_t & nth3,
			                                float   & tr_start3,
							float   & tr_stop3,
							float   & dt_rad3,
							int32_t & nph1,
							float   & pr_start1,
							float   & pr_stop1,
							float   & dp_rad1,
							int32_t & nph2,
							float   & pr_start2,
							float   & pr_stop2,
							float   & dp_rad2,
							int32_t & nph3,
							float   & pr_start3,
							float   & pr_stop3,
							float   & dp_rad3) {
                                  float   td_start1,td_stop1,dt_deg1, 
                                          td_start2,td_stop2,dt_deg2, 
                                          td_start3,td_stop3,dt_deg3, 
                                          pd_start1,pd_stop1,dp_deg1, 
                                          pd_start2,pd_stop2,dp_deg2, 
                                          pd_start3,pd_stop3,dp_deg3;
				  const float t0 =  0.017453292519943f;
				  td_start1 = 2.5f;
                                  td_stop1 = 177.5f;
                                  dt_deg1  = 5.0f;
                                  nth1  = 35;
                                  tr_start1 = t0*td_start1
                                  tr_stop1  = t0*td_stop1
                                  dt_rad1   = t0*dt_deg1
                                  pd_start1 = 2.5f;
                                  pd_stop1 = 177.5f;
                                  dp_deg1 = 5.0f;
                                  nph1 = 36;
                                  pr_start1 = t0*pd_start1;
                                  pr_stop1  = t0*pd_stop1;
                                  dp_rad1   = t0*dp_deg1;
                                  td_start2 = 0.0f;
                                  td_stop2 = 0.0f;
                                  dt_deg2 = 0.0f;
                                  nth2 = 0;
                                  tr_start2 = 0.0f;
                                  tr_stop2  = 0.0f;
                                  dt_rad2   = 0.0f;
                                  pd_start2 = 0.0f;
                                  pd_stop2 = 0.0f;
                                  dp_deg2 = 0.0f;
                                  nph2 = 0;
                                  pr_start2 = 0.0f;
                                  pr_stop2 = 0.0f;
                                  dp_rad2 = 0.0f;
                                  td_start3 = 0.0f
                                  td_stop3 = 0.0f;
                                  dt_deg3 = 0.0f;
                                  nth3 = 0;
                                  tr_start3 = 0.0f;
                                  tr_stop3  = 0.0f;
                                  dt_rad3 = 0.0f;
                                  pd_start3 = 0.0f;
                                  pd_stop3 = 0.0f;
                                  pd_deg3 = 0.0f;
                                  nph3 = 0;
                                  pr_start3 = 0.0f;
                                  pr_stop3 = 0.0f;
                                  dp_rad3 = 0.0f;

		    }

			__ATTR_HOT__
			__ATTR_ALIGN__(64)
			inline
		        std::complex<float>
			Leaf_dielectric(const float leaf_mg,
			                const float leaf_rho,
			                const float leaf_dens,
					const float leaf_diam,
					const float leaf_tau,
					const bool dry_dens,
				        const float water_tmp,
		                        const float veg_tmp,
					const float theta,
					const float rad_freq)  {
                           if(dry_dens) {
                                return (Veg_dielectric_2(leaf_mg,
				                         leaf_rho,
							 veg_tmp,
							 theta,
							 rad_freq));
			   }
			   else {
                                return (Veg_dielectric_1(leaf_mg,
				                         veg_tmp,
							 theta,
							 rad_freq));
			   }
		     }

                        __ATTR_HOT__
			__ATTR_ALIGN__(64)
			inline
		        std::complex<float>
			Veg_dielectric_2(const float mg,
			                 const float veg_rho,
					 const float tempC,
					 const float theta,
					 const float rad_freq) {
                             std::complex<float> e,f,g,w;
			     std::complex<float> result;
			     float mv,a,b,c,d;
			     float top,fn,en,ein;
			     float t0;
			     mv  = mg*veg_rho/(1.0f*(1.0f-veg_rho));
			     t0  = mv*mv;
			     a   = 1.7f+3.20f*mv+6.5f*t0;
			     top =  1.1109e-10f+tempC*(-3.824e-12f+tempC* 
                                    (6.938e-14f-tempC*5.096e-16f));
			     b   = mv*(0.82f*mv+0.166f);
			     fn  = 1.0f/(top*1.09f);
			     e   = {1.0f,(rad_freq/fn)};
			     d   = 22.74f;
			     c   = 31.4f*t0/(59.5f*t0+1.0f);
			     f   = {0.0f,(d_rad_freq)};
			     en  =  88.045f+tempC*(-0.4147f+tempC*(6.295e-4f +
                                    tempC*1.075e-5f));
			     result = {};
			     ein = 4.9f;
			     w   = 0.707106781186548f*{1.0f,1.0f}*std::sqrt(rad_freq/0.18_sp);
			     g   = 1.0f*w;
			     result = a+b*(4.9f+(en-ein)/e-f)+c*(2.9f+55.0f/g);
			     return (result);
		     }
		  
                     __ATTR_HOT__
		     __ATTR_ALIGN__(64)
		     inline
		     std::complex<float>
		     Veg_dielectric_1(const float mg,
		                      const float tempC,
				      const float theta,
				      const float rad_freq) {
                         std::complex<float> e,f,g,w;
			 std::complex<float> result;
			 float top,fn,en,ein,t0;
			 float a,b,c,d;
			 t0  = mg*mg;
			 a   = 1.7f-0.74f*mg+6.16f*t0;
			 top = 1.1109e-10f+tempC*(-3.824e-12f+tempC * 
                               (6.938e-14f-tempC*5.096e-16f));
			 b   = mg*(0.55f*mg-0.076f);
			 fn  = 1.0f/(top*1.09f);
			 e   = {1.0f,(rad_freq/fn)};
			 c   = 4.64f*t0/(7.36f*t0+1.0f);
			 d   = 22.74f;
			 f   = {0.0f,(d/rad_freq)};
			 en  = 88.045f+tempC*(-0.4147f+tempC*(6.295e-4f +
                               tempC*1.075e-5f));
			 w   = 0.707106781186548f*{1.0f,1.0f}*std::sqrt(rad_freq/0.18_sp);
			 ein = 4.9f;
			 result = {};
			 g   = 1.0f*w;
			 result = a+b*(4.9f+(en-ein)/e-f)+c*(2.9f+55.0f/g);
			 return (result);
		    }

		    

		    __ATTR_HOT__				
                    __ATTR_ALIGN__(32)
		    inline
		    float
		    Leaf_ang_orientation(const float mu,
		                         const float nu,
					 const float th) {
                         float t0,t1,t2;
			 float result;
			 const float two_over_pi =  0.636619772367581f;
			 const float half_pi     =  1.570796326794897f;
			 t0 = two_over_pi*exp(gamma(mu+nu)-gamma(mu)-gamma(nu));
			 t1 = 1.0f-th/std::pow(half_pi,(nu-1.0f));
			 result = 0.0f;
			 t2 = th/std::pow(half_pi,(mu-1.0f));
			 result = t0*t1*t2;
		    }
							  


		} __ATTR_ALIGN__(64);
		


     } // math



} // gms



















#endif /*__GMS_TREE_SCATTERER_H__*/
