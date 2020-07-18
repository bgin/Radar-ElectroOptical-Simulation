

#ifndef __GMS_TREE_SCATTERER_COMMON_H__
#define __GMS_TREE_SCATTERER_COMMON_H__



namespace file_info {

     const unsigned int gGMS_TREE_SCATTERER_COMMON_MAJOR = 1;
     const unsigned int gGMS_TREE_SCATTERER_COMMON_MINOR = 0;
     const unsigned int gGMS_TREE_SCATTERER_COMMON_MICRO = 0;
     const unsigned int gGMS_TREE_SCATTERER_COMMON_FULLVER =
       1000U*gGMS_TREE_SCATTERER_COMMON_MAJOR+100U*gGMS_TREE_SCATTERER_COMMON_MINOR+
       10U*gGMS_TREE_SCATTERER_COMMON_MICRO;
     const char * const pgGMS_TREE_SCATTERER_COMMON_CREATION_DATE = "12-07-2020 10:12 +00200 (SUN 12 JUL 2020 GMT+2)";
     const char * const pgGMS_TREE_SCATTERER_COMMON_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const pgGMS_TREE_SCATTERER_COMMON_SYNOPSIS      = "Leaf phase matrices implementation."
}

#include <cstdint>
#include <complex>
#include <math.h>
#include <immintrin.h>
#include "GMS_config.h"

#if !defined LEAF_PHASE_MATRICES_AUTOVECTORIZE
    #define LEAF_PHASE_MATRICES_AUTOVECTORIZE 1
#endif
namespace gms {

         namespace math {


	                __ATTR_HOT__
		        __ATTR_ALIGN__(32)
			void
                        Leaf_phase_matrices(float * __restrict __ATTR_ALIGN__(64),
			                    std::complex<float> * __restrict __ATTR_ALIGN__(32),
					    std::complex<float> * __restrict __ATTR_ALIGN__(32),
					    float * __restrict __ATTR_ALIGN__(64),
					    std::complex<float> * __restrict __ATTR_ALIGN__(32),
					    const float,
					    const float,
					    const float,
					    const float,
					    const float,
					    const float,
					    const float,
					    const float,
					    const float,
					    const float,
					    const float,
					    const std::complex<float>,
					    const float,
					    const float,
					    const float,
					    const int32_t);
					    
			
	               	__ATTR_COLD__
			__ATTR_ALIGN__(32)
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


		       
                       float
		       Compute_leaf_odf(const int32_t,
		                        const float) __ATTR_HOT__ __ATTR_ALIGN__(32);

		       __ATTR_HOT__				
                       __ATTR_ALIGN__(64)
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

		    void
		    Leaf_PO_approximation(const float,
		                          const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const float,
					  const std::complex<float>,
					  std::complex<float> * __restrict __ATTR_ALIGN__(32)) __ATTR_HOT__ __ATTR_ALIGN__(32);

		    void
		    Leaf_Rayleigh_scattering(const float,
		                             const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const std::complex<float>,
					     std::complex<float> * __restrict __ATTR_ALIGN__(32)) __ATTR_HOT__ __ATTR_ALIGN__(32);

		    __ATTR_HOT__
		    __ATTR_ALIGN__(64)
		    __ATTR_VECTORCALL__
		    inline
		    void
		    vec1x3_smooth(float * __restrict __ATTR_ALIGN__(16) in,
		                  float * __restrict __ATTR_ALIGN__(16) out) {
			 __m128 tmp;
                         float mag;
			 mag = 0.0f;
#if defined __INTEL_COMPILER
                         __assume_aligned(in,16);
			 __assume_aligned(out,16);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                         in = (float*)__builtin_assume_aligned(in,16);
			 out = (float*)__builtin_assume_aligned(out,16);
#endif
			 mag = sqrtf(in[0]*in[0]+in[1]*in[1]+in[2]*in[2]);
			 if(mag!=0.0f) {
			    tmp = _mm_setzero_ps();
			    tmp = _mm_load_ps(&in[0]);
                            _mm_store_ps(&in[0],_mm_div_ps(tmp,_mm_set1_ps(mag));
			    check_mag_tol(in,out);
			    mag = sqrtf(in[0]*in[0]+in[1]*in[1]+in[2]*in[2]);
			    _mm_store_ps(&out[0],_mm_div_ps(_mm_load_ps(&in[0]),
			                                     _mm_set1_ps(mag)));
			 }
			 else {
                            _mm_store_ps(&out[0],_mm_load_ps(&in[0]));
			 }
		  }

		 __ATTR_HOT__
                 __ATTR_ALIGN__(64)
		 __ATTR_VECTORCALL__
		 inline
		 void
		 check_mag_tol(const float * __restrict __ATTR_ALIGN__(16) in,
		               float * __restrict __ATTR_ALIGN__(16) out) {
		      __m128 tmp,tol;
                      tol = _mm_set1_ps(0.00001f);
		      tmp = _mm_setzero_ps();
		      tmp = _mm_cmplt_ps(_mm_load_ps(&in[0]),tol);
		      if(!_mm_test_all_ones(_mm_castps_si128(tmp))) {
                          _mm_store_ps(&out[0],_mm_setzero_ps());
		      }
		      else {
                          _mm_store_ps(&out[0],_mm_load_ps(&in[0]));
		      }
		 }

                 __ATTR_HOT__
		 __ATTR_ALIGN__(64)
		 __ATTR_VECTORCALL__
		 inline
		 void
		 cross_prod(float * __restrict __ATTR_ALIGN__(16) a,
		            float * __restrict __ATTR_ALIGN__(16) b,
			    float * __restrict __ATTR_ALIGN__(16) c) {
#if defined __INTEL_COMPILER
                         __assume_aligned(a,16);
			 __assume_aligned(b,16);
			 __assume_aligned(c,16);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                        a = (float*)__builtin_assume_aligned(a,16);
			b = (float*)__builtin_assume_aligned(b,16);
			c = (float*)__builtin_assume_aligned(c,16);
#endif
                        c[0] = a[1]*b[2]-a[2]*b[1];
			c[1] = a[2]*b[0]-a[0]*b[2];
			c[2] = a[0]*b[1]-a[1]*b[0];
		}

		__ATTR_HOT__
		__ATTR_ALIGN__(64)
		__ATTR_VECTORCALL__
		inline
		float dot_prod(const float * __restrict __ATTR_ALIGN__(16) a,
		               const float * __restrict __ATTR_ALIGN__(16) b) {
#if defined __INTEL_COMPILER
                         __assume_aligned(a,16);
			 __assume_aligned(b,16);
			
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                        a = (const float*)__builtin_assume_aligned(a,16);
			b = (const float*)__builtin_assume_aligned(b,16);
		
#endif
                        float res;
			res = 0.0f;
			res = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
			return (res);
		}

		__ATTR_HOT__
		__ATTR_ALIGN__(64)
		__ATTR_VECTORCALL__
		inline
        	void
		stokes_matrix(std::complex<float>* __restrict __ATTR_ALIGN__(32) scat_mat,
		              float * __restrict __ATTR_ALIGN__(64) stokes_mat) {

		      std::complex<float> CW1121C,CW1122C,CW1112C,CW2122C,
                                          CW1221C,CW1222C,CW1,CW2;
		      float  w1,w2,w3,w4;
#if defined __INTEL_COMPILER
                      __assume_aligned(scat_mat,32);
		      __assume_aligned(stokes_mat,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                      scat_mat   = (std::complex<float>*)__builtin_assume_aligned(scat_mat,32);
		      stokes_mat = (float*)__builtin_assume_aligned(stokes_mat,64);
#endif
		      w1 = std::abs(scat_mat[0]);
		      stokes_mat[0] = w1*w1;
		      CW1121C = scat_mat[0]*std::conj(scat_mat[2]);
		      w2 = std::abs(scat_mat[1]);
		      stokes_mat[1] = w2*w2;
		      CW1222C = scat_mat[1]*std::conj(scat_mat[3]);
		      w3 = std::abs(scat_mat[2]);
		      stokes_mat[2] = w3*w3;
		      CW1112C = scat_mat[0]*std::conj(scat_mat[1]);
		      w4 = std::abs(scat_mat[3]);
		      stokes_mat[3] = w4*w4;
		      CW2122C = scat_mat[2]*std::conj(scat_mat[3]);
		      CW1122C = scat_mat[0]*std::conj(scat_mat[3]);
		      CW1221C = scat_mat[1]*std::conj(scat_mat[2]);
		      CW1 =  CW1122C + CW1221C;
		      CW2 =  CW1122C - CW1221C;
		      stokes_mat[4] = 2.0f*CW1121C.real();
		      stokes_mat[5] = 2.0f*CW1121C.imag();
		      stokes_mat[6] = 2.0f*CW1222C.real();
                      stokes_mat[7] = 2.0f*CW1222C.imag();
		      stokes_mat[8] = CW1112C.real();
		      stokes_mat[9] = CW2122C.imag();
		      stokes_mat[10] = CW1.real();
		      stokes_mat[11] = CW1.imag();
		      stokes_mat[12] = -CW1112C.imag();
		      stokes_mat[13] = -CW2122C.imag();
		      stokes_mat[14] = -CW2.imag();
		      stokes_mat[15] = CW2.real();
		}

     } // math

} // gms




#endif /*__GMS_TREE_SCATTERER_COMMON_H__*/
