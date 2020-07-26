

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
		        
			 // Allocate this array as a [4*4*4*nleaves]
			 float * __restrict __ATTR_ALIGN__(8)  l4x4phm;
		         // Allocate this array as a [2*2*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) l2x2mp;
			 // Allocate this array as a [2*2*nleaves]
			 std::complex<float> * __restrict __ATTR_ALIGN__(8) l2x2mn;
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
							// Length of arrays is the number of leaves.
                        void ComputeLeafPhaseMatrices(const float * __restrict __ATTR_ALIGN__(64),
			                              const float * __restrict __ATTR_ALIGN__(64),
						      const float * __restrict __ATTR_ALIGN__(64),
						      const float * __restrict __ATTR_ALIGN__(64),
						      const float * __restrict __ATTR_ALIGN__(64),
						      const std::complex<float> * __restrict __ATTR_ALIGN__(64),
						      const int32_t * __restrict __ATTR_ALIGN__(64),
		                                      const float,
						      const float,
						      const float,
						      const float,
						      const float,
						      const float) __ATTR_HOT__ __ATTR_ALIGN__(32);
		

                      
		  
                   

		    

		 				  


		} __ATTR_ALIGN__(64);
		


     } // math



} // gms



















#endif /*__GMS_TREE_SCATTERER_H__*/
