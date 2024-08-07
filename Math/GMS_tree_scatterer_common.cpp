
#include <limits>
#if defined __GNUC__ && !defined __INTEL_COMPILER
#include <omp.h>
#endif
#include "GMS_tree_scatterer_common.h"
#include "GMS_indices.h"


void
gms::math
::Leaf_phase_matrices(float * __restrict __ATTR_ALIGN__(64) l4x4phm,
		      std::complex<float> * __restrict __ATTR_ALIGN__(32) l2x2mp,
		      std::complex<float> * __restrict __ATTR_ALIGN__(32) l2x2mn,
		      float * __restrict __ATTR_ALIGN__(64) stokes4x4m,
		      std::complex<float> * __restrict __ATTR_ALIGN__(32) scat2x2m,
		      const float ldiam,
		      const float lthick,
		      const float lmg,
		      const float lrho,
		      const float ldens,
		      const float theta,
		      const float ctheta,
		      const float stheta,
		      const std::complex<float> epsr,
		      const float rad_freq,
		      const float rad_wv,
		      const float rad_k0,
		      const int32_t leaf_orient) {

     __attribute__((aligned(64))) float l4x4phm_t1[4][4][4] = {};
     __attribute__((aligned(64))) float l4x4phm_t2[4][4][4] = {};
     __attribute__((aligned(64))) float l4x4phm_t3[4][4][4] = {};
     __attribute__((aligned(64))) std::complex<float> sm2x2avg_t1[2][2][2] = {0.0f,0.0f};
     __attribute__((aligned(64))) std::complex<float> sm2x2avg_t2[2][2][2] = {0.0f,0.0f};
     __attribute__((aligned(64))) std::complex<float> sm2x2avg_t3[2][2][2] = {0.0f,0.0f};
     std::complex<float> j,cwork;
     float  tr_start1,tr_stop1,dt_rad1,
            tr_start2,tr_stop2,dt_rad2,
            tr_start3,tr_stop3,dt_rad3,
            pr_start1,pr_stop1,dp_rad1,
            pr_start2,pr_stop2,dp_rad2,
            pr_start3,pr_stop3,dp_rad3;
     float  dp1t1,dp2t2,dp3t3,orient_distr,work;
     float  t0,norm,thinc,phinc,thsc,phsc,thdr,phdr;
     float  t1,t2,t3,t4,t5,t6,t7,t8;
     float  t9,t10,t11,t12,t13,t14,t15,t16;
     int32_t nth1,nth2,nth3,nph1,nph2,nph3;
     int32_t ii,j,l,k,jj;
     bool  po;
     t0 = ldiam/100.0f;
     norm = 6.283185307179586f;
     nth1 = 0;
     tr_start1 = 0.0f;
     tr_stop1  = 0.0f;
     dt_rad1   = 0.0f;
     nph1 = 0;
     pr_start1 = 0.0f;
     pr_stop1  = 0.0f;
     dp_rad1   = 0.0f;
     if((rad_wv/t0)<1.5f)
        po = true;
     else
        po = false
     Set_leaf_quadrature_bounds(nth1,
			        tr_start1,
				tr_stop1,
				dt_rad1,
			        nth2,
			        tr_start2,
				tr_stop2,
			        dt_rad2,
				nth3,
			        tr_start3,
				tr_stop3,
				dt_rad3,
				nph1,
				pr_start1,
				pr_stop1,
				dp_rad1,
				nph2,
				pr_start2,
				pr_stop2,
				dp_rad2,
				nph3,
				pr_start3,
				pr_stop3,
				dp_rad3);
       dp1t1 = dp_rad1*dt_rad1;
       dp2t2 = dp_rad2*dt_rad2;
       dp3t3 = dp_rad3*dt_rad3;
       orient_distr = 0.0f;
       if((nth1!=0) && (nph1!=0)) {
           t1=0.0f;
           t2=0.0f;
           t3=0.0f;
           t4=0.0f;
           t5=0.0f;
           t6=0.0f;
           t7=0.0f;
           t8=0.0f;
           t9=0.0f;
           t10=0.0f;
           t11=0.0f;
           t12=0.0f;
           t13=0.0f;
           t14=0.0f;
           t15=0.0f;
           t16=0.0f;
	   thdr = 0.0f;
           for(jj=0; jj != nth1; ++jj) {
               thdr = tr_start1+dt_rad1*static_cast<float>(jj);
	       orient_distr = Compute_leaf_odf(leaf_orient,thdr);
	       if(orient_distr>0.0f){
	          phdr  = 0.0f;
		  thinc = 0.0f;
		  thsc  = 0.0f;
		  phinc = 0.0f;
		  phsc  = 0.0f;
                  for(ii=0; ii != nph1; ++ii) {
                      phdr  = tr_start1+dt_rad1*static_cast<float>(ii);
		      thinc = theta;
		      thsc  = 3.141592653589793f-theta;
		      phinc = 3.141592653589793f;
		      if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		      }
		       else {
                             Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		      stokes_matrix(&scat2x2m[0],
		                    &stokes4x4m[0]);
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1		      
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][0]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][0] = t1;
			   }
		       }
#else
                        t1 = l4x4phm_t1[0][0][0]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t1[0][0][0] = t1;
                        t2 =  l4x4phm_t1[0][1][0]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t1[0][1][0] = t2;
                        t3 =  l4x4phm_t1[0][2][0]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t1[0][2][0] = t3;
                        t4 =  l4x4phm_t1[0][3][0]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t1[0][3][0] = t4;
                        t5 =  l4x4phm_t1[1][0][0]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t1[1][0][0] = t5;
                        t6 =  l4x4phm_t1[1][1][0]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t1[1][1][0] = t6;
                        t7 =  l4x4phm_t1[1][2][0]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t1[1][2][0] = t7;
                        t8 =  l4x4phm_t1[1][3][0]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t1[1][3][0] = t8;
                        t9 =  l4x4phm_t1[2][0][0]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t1[2][0][0] = t9;
                        t10 = l4x4phm_t1[2][1][0]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t1[2][1][0] = t10;
                        t11 = l4x4phm_t1[2][2][0]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t1[2][2][0] = t11;
                        t12 = l4x4phm_t1[2][3][0]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t1[2][3][0] = t12;
                        t13 = l4x4phm_t1[3][0][0]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t1[3][0][0] = t13;
                        t14 = l4x4phm_t1[3][1][0]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t1[3][1][0] = t14;
                        t15 = l4x4phm_t1[3][2][0]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t1[3][2][0] = t15;
                        t16 = l4x4phm_t1[3][3][0]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t1[3][3][0] = t16;
                      
#endif
                       scat2x2m[1] = -scat2x2m[1];
		       scat2x2m[2] = -scat2x2m[2];
		       stokes_matrix(&scat2x2m[0],
		                     &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][0]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][0] = t1;
			   }
		       }
#else
                           t1 = l4x4phm_t1[0][0][0]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t1[0][0][0] = t1;
                        t2 =  l4x4phm_t1[0][1][0]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t1[0][1][0] = t2;
                        t3 =  l4x4phm_t1[0][2][0]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t1[0][2][0] = t3;
                        t4 =  l4x4phm_t1[0][3][0]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t1[0][3][0] = t4;
                        t5 =  l4x4phm_t1[1][0][0]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t1[1][0][0] = t5;
                        t6 =  l4x4phm_t1[1][1][0]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t1[1][1][0] = t6;
                        t7 =  l4x4phm_t1[1][2][0]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t1[1][2][0] = t7;
                        t8 =  l4x4phm_t1[1][3][0]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t1[1][3][0] = t8;
                        t9 =  l4x4phm_t1[2][0][0]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t1[2][0][0] = t9;
                        t10 = l4x4phm_t1[2][1][0]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t1[2][1][0] = t10;
                        t11 = l4x4phm_t1[2][2][0]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t1[2][2][0] = t11;
                        t12 = l4x4phm_t1[2][3][0]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t1[2][3][0] = t12;
                        t13 = l4x4phm_t1[3][0][0]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t1[3][0][0] = t13;
                        t14 = l4x4phm_t1[3][1][0]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t1[3][1][0] = t14;
                        t15 = l4x4phm_t1[3][2][0]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t1[3][2][0] = t15;
                        t16 = l4x4phm_t1[3][3][0]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t1[3][3][0] = t16;

#endif

                  // Case 2
		  thinc =  3.141592653589793f-theta;
                  thsc  =  3.141592653589793f-theta;
                  phinc =  0.0f
                  phsc  =  3.141592653589793f;
		  if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		    stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][1]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][1] = t1;
			   }
		       }
#else
                         t1 = l4x4phm_t1[0][0][1]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t1[0][0][1] = t1;
                         t2 =  l4x4phm_t1[0][1][1]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t1[0][1][1] = t2;
                         t3 =  l4x4phm_t1[0][2][1]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t1[0][2][1] = t3;
                         t4 =  l4x4phm_t1[0][3][1]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t1[0][3][1] = t4;
                         t5 =  l4x4phm_t1[1][0][1]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t1[1][0][1] = t5;
                         t6 =  l4x4phm_t1[1][1][1]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t1[1][1][1] = t6;
                         t7 =  l4x4phm_t1[1][2][1]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t1[1][2][1] = t7;
                         t8 =  l4x4phm_t1[1][3][1]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t1[1][3][1] = t8;
                         t9 =  l4x4phm_t1[2][0][1]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t1[2][0][1] = t9;
                         t10 = l4x4phm_t1[2][1][1]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t1[2][1][1] = t10;
                         t11 = l4x4phm_t1[2][2][1]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t1[2][2][1] = t11;
                         t12 = l4x4phm_t1[2][3][1]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t1[2][3][1] = t12;
                         t13 = l4x4phm_t1[3][0][1]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t1[3][0][1] = t13;
                         t14 = l4x4phm_t1[3][1][1]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t1[3][1][1] = t14;
                         t15 = l4x4phm_t1[3][2][1]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t1[3][2][1] = t15;
                         t16 = l4x4phm_t1[3][3][1]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t1[3][3][1] = t16;

#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][1]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][1] = t1;
			   }
		       }

#else
                        t1 = l4x4phm_t1[0][0][1]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t1[0][0][1] = t1;
                         t2 =  l4x4phm_t1[0][1][1]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t1[0][1][1] = t2;
                         t3 =  l4x4phm_t1[0][2][1]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t1[0][2][1] = t3;
                         t4 =  l4x4phm_t1[0][3][1]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t1[0][3][1] = t4;
                         t5 =  l4x4phm_t1[1][0][1]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t1[1][0][1] = t5;
                         t6 =  l4x4phm_t1[1][1][1]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t1[1][1][1] = t6;
                         t7 =  l4x4phm_t1[1][2][1]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t1[1][2][1] = t7;
                         t8 =  l4x4phm_t1[1][3][1]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t1[1][3][1] = t8;
                         t9 =  l4x4phm_t1[2][0][1]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t1[2][0][1] = t9;
                         t10 = l4x4phm_t1[2][1][1]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t1[2][1][1] = t10;
                         t11 = l4x4phm_t1[2][2][1]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t1[2][2][1] = t11;
                         t12 = l4x4phm_t1[2][3][1]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t1[2][3][1] = t12;
                         t13 = l4x4phm_t1[3][0][1]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t1[3][0][1] = t13;
                         t14 = l4x4phm_t1[3][1][1]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t1[3][1][1] = t14;
                         t15 = l4x4phm_t1[3][2][1]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t1[3][2][1] = t15;
                         t16 = l4x4phm_t1[3][3][1]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t1[3][3][1] = t16;

#endif
                    // Case 3
		    thinc = theta;
                    thsc  = theta;
                    phinc = 3.141592653589793f;
                    phsc  = 0.0f;
		    if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		    stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][2]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][2] = t1;
			   }
		       }
#else
                         t1 = l4x4phm_t1[0][0][2]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t1[0][0][2] = t1;
                         t2 =  l4x4phm_t1[0][1][2]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t1[0][1][2] = t2;
                         t3 =  l4x4phm_t1[0][2][2]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t1[0][2][2] = t3;
                         t4 =  l4x4phm_t1[0][3][2]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t1[0][3][2] = t4;
                         t5 =  l4x4phm_t1[1][0][2]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t1[1][0][2] = t5;
                         t6 =  l4x4phm_t1[1][1][2]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t1[1][1][2] = t6;
                         t7 =  l4x4phm_t1[1][2][2]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t1[1][2][2] = t7;
                         t8 =  l4x4phm_t1[1][3][2]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t1[1][3][2] = t8;
                         t9 =  l4x4phm_t1[2][0][2]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t1[2][0][2] = t9;
                         t10 = l4x4phm_t1[2][1][2]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t1[2][1][2] = t10;
                         t11 = l4x4phm_t1[2][2][2]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t1[2][2][2] = t11;
                         t12 = l4x4phm_t1[2][3][2]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t1[2][3][2] = t12;
                         t13 = l4x4phm_t1[3][0][2]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t1[3][0][2] = t13;
                         t14 = l4x4phm_t1[3][1][2]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t1[3][1][2] = t14;
                         t15 = l4x4phm_t1[3][2][2]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t1[3][2][2] = t15;
                         t16 = l4x4phm_t1[3][3][2]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t1[3][3][2] = t16;

#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][2]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][2] = t1;
			   }
		       }
#else
                       t1 = l4x4phm_t1[0][0][2]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t1[0][0][2] = t1;
                         t2 =  l4x4phm_t1[0][1][2]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t1[0][1][2] = t2;
                         t3 =  l4x4phm_t1[0][2][2]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t1[0][2][2] = t3;
                         t4 =  l4x4phm_t1[0][3][2]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t1[0][3][2] = t4;
                         t5 =  l4x4phm_t1[1][0][2]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t1[1][0][2] = t5;
                         t6 =  l4x4phm_t1[1][1][2]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t1[1][1][2] = t6;
                         t7 =  l4x4phm_t1[1][2][2]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t1[1][2][2] = t7;
                         t8 =  l4x4phm_t1[1][3][2]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t1[1][3][2] = t8;
                         t9 =  l4x4phm_t1[2][0][2]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t1[2][0][2] = t9;
                         t10 = l4x4phm_t1[2][1][2]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t1[2][1][2] = t10;
                         t11 = l4x4phm_t1[2][2][2]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t1[2][2][2] = t11;
                         t12 = l4x4phm_t1[2][3][2]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t1[2][3][2] = t12;
                         t13 = l4x4phm_t1[3][0][2]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t1[3][0][2] = t13;
                         t14 = l4x4phm_t1[3][1][2]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t1[3][1][2] = t14;
                         t15 = l4x4phm_t1[3][2][2]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t1[3][2][2] = t15;
                         t16 = l4x4phm_t1[3][3][2]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t1[3][3][2] = t16;
#endif

                    // Case 4:
		    thinc = 3.141592653589793f-theta;
                    thsc  = theta;
                    phinc = 0.0f;
                    phsc  = 3.141592653589793f;
		    if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		    stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][3]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][3] = t1;
			   }
		       }
#else
                        t1 = l4x4phm_t1[0][0][3]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t1[0][0][3] = t1;
                        t2 =  l4x4phm_t1[0][1][3]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t1[0][1][3] = t2;
                        t3 =  l4x4phm_t1[0][2][3]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t1[0][2][3] = t3;
                        t4 =  l4x4phm_t1[0][3][3]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t1[0][3][3] = t4;
                        t5 =  l4x4phm_t1[1][0][3]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t1[1][0][3] = t5;
                        t6 =  l4x4phm_t1[1][1][3]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t1[1][1][3] = t6;
                        t7 =  l4x4phm_t1[1][2][3]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t1[1][2][3] = t7;
                        t8 =  l4x4phm_t1[1][3][3]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t1[1][3][3] = t8;
                        t9 =  l4x4phm_t1[2][0][3]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t1[2][0][3] = t9;
                        t10 = l4x4phm_t1[2][1][3]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t1[2][1][3] = t10;
                        t11 = l4x4phm_t1[2][2][3]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t1[2][2][3] = t11;
                        t12 = l4x4phm_t1[2][3][3]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t1[2][3][3] = t12;
                        t13 = l4x4phm_t1[3][0][3]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t1[3][0][3] = t13;
                        t14 = l4x4phm_t1[3][1][3]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t1[3][1][3] = t14;
                        t15 = l4x4phm_t1[3][2][3]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t1[3][2][3] = t15;
                        t16 = l4x4phm_t1[3][3][3]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t1[3][3][3] = t16;
#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t1[k][l][3]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][3] = t1;
			   }
		       }

#else
                        t1 = l4x4phm_t1[0][0][3]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t1[0][0][3] = t1;
                        t2 =  l4x4phm_t1[0][1][3]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t1[0][1][3] = t2;
                        t3 =  l4x4phm_t1[0][2][3]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t1[0][2][3] = t3;
                        t4 =  l4x4phm_t1[0][3][3]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t1[0][3][3] = t4;
                        t5 =  l4x4phm_t1[1][0][3]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t1[1][0][3] = t5;
                        t6 =  l4x4phm_t1[1][1][3]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t1[1][1][3] = t6;
                        t7 =  l4x4phm_t1[1][2][3]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t1[1][2][3] = t7;
                        t8 =  l4x4phm_t1[1][3][3]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t1[1][3][3] = t8;
                        t9 =  l4x4phm_t1[2][0][3]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t1[2][0][3] = t9;
                        t10 = l4x4phm_t1[2][1][3]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t1[2][1][3] = t10;
                        t11 = l4x4phm_t1[2][2][3]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t1[2][2][3] = t11;
                        t12 = l4x4phm_t1[2][3][3]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t1[2][3][3] = t12;
                        t13 = l4x4phm_t1[3][0][3]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t1[3][0][3] = t13;
                        t14 = l4x4phm_t1[3][1][3]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t1[3][1][3] = t14;
                        t15 = l4x4phm_t1[3][2][3]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t1[3][2][3] = t15;
                        t16 = l4x4phm_t1[3][3][3]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t1[3][3][3] = t16;

#endif
                    // Extinction matrix: case 1
		    thinc = theta;
                    thsc  = thinc;
                    phinc = 3.141592653589793f;
                    phsc  = phinc;
		    if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                          for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                               __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                               scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                  for(l=0; l != 2; ++l) {
                                      t1 = sm2x2avg_t1[k][l][0]+orient_distr*scat2x2m[k][l];
		                      sm2x2avg_t1[k][l][0] = t1;
		                  }
	                  }
#else
               t1 = sm2x2avg_t1[0][0][0]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t1[0][0][0] = t1;
	       t2 = sm2x2avg_t1[0][1][0]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t1[0][1][0] = t1;
	       t3 = sm2x2avg_t1[1][0][0]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t1[1][0][0] = t3;
	       t4 = sm2x2avg_t1[1][1][0]+orient_distr*scat2x2m[1][1];
               sm2x2avg_t1[1][1][0] = t4;
#endif

                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                          for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                               __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                               scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                  for(l=0; l != 2; ++l) {
                                      t1 = sm2x2avg_t1[k][l][0]+orient_distr*scat2x2m[k][l];
		                      sm2x2avg_t1[k][l][0] = t1;
		                  }
	                  }
#else
               t1 = sm2x2avg_t1[0][0][0]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t1[0][0][0] = t1;
	       t2 = sm2x2avg_t1[0][1][0]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t1[0][1][0] = t1;
	       t3 = sm2x2avg_t1[1][0][0]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t1[1][0][0] = t3;
	       t4 = sm2x2avg_t1[1][1][0]+orient_distr*scat2x2m[1][1];
               sm2x2avg_t1[1][1][0] = t4;
#endif	     

                    // Extinction matrix: case 2
		    thinc =   3.141592653589793f-theta;
                    thsc  =   3.141592653589793f-theta;
                    phinc =   0.0f;
                    phsc  =   phinc;
		    if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                               for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                     __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                   scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                      for(l=0; l != 2; ++l) {
                                          t1 = sm2x2avg_t1[k][l][1]+orient_distr*scat2x2m[k][l];
		                          sm2x2avg_t1[k][l][1] = t1;
		                      }
	                        }
#else
               t1 = sm2x2avg_t1[0][0][1]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t1[0][0][1] = t1;
	       t2 = sm2x2avg_t1[0][1][1]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t1[0][1][1] = t1;
	       t3 = sm2x2avg_t1[1][0][1]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t1[1][0][1] = t3;
               t4 = sm2x2avg_t1[1][1][1]+orient_distr*scat2x2m[1][1];
	       sm2x2avg_t1[1][1][1] = t4;
#endif

                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                               for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                     __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                   scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                      for(l=0; l != 2; ++l) {
                                          t1 = sm2x2avg_t1[k][l][1]+orient_distr*scat2x2m[k][l];
		                          sm2x2avg_t1[k][l][1] = t1;
		                      }
	                        }
#else
               t1 = sm2x2avg_t1[0][0][1]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t1[0][0][1] = t1;
	       t2 = sm2x2avg_t1[0][1][1]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t1[0][1][1] = t1;
	       t3 = sm2x2avg_t1[1][0][1]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t1[1][0][1] = t3;
               t4 = sm2x2avg_t1[1][1][1]+orient_distr*scat2x2m[1][1];
	       sm2x2avg_t1[1][1][1] = t4;
#endif		     
	     
		  } // end for(ii=0 loop
	       } // end if(orient_distr block
	   } // end for(jj=0 loop
       } // end nth1!=0 block

       if((nth2!=0) && (nph2!=0)) {
           t1=0.0f;
           t2=0.0f;
           t3=0.0f;
           t4=0.0f;
           t5=0.0f;
           t6=0.0f;
           t7=0.0f;
           t8=0.0f;
           t9=0.0f;
           t10=0.0f;
           t11=0.0f;
           t12=0.0f;
           t13=0.0f;
           t14=0.0f;
           t15=0.0f;
           t16=0.0f;
	   for(jj=0; jj != nth2; ++jj) {
               thdr = tr_start2+dt_rad2*static_cast<float>(jj);
	       orient_distr = Compute_leaf_odf(leaf_orient,thdr);
	       if(orient_distr>0.0f) {
                  for(ii=0; ii != nph2; ++ii) {
                      phdr = pr_start2+dp_rad2*static_cast<float>(ii);
		      thinc = theta;
                      thsc  = 3.141592653589793f-theta;
                      phinc = 3.141592653589793f;
                      phsc  = 0.0f;
		      if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		       stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1				  
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][0]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][0] = t1;
			   }
		       }
#else
                         t1 = l4x4phm_t2[0][0][0]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t2[0][0][0] = t1;
                         t2 =  l4x4phm_t2[0][1][0]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t2[0][1][0] = t2;
                         t3 =  l4x4phm_t2[0][2][0]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t2[0][2][0] = t3;
                         t4 =  l4x4phm_t2[0][3][0]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t2[0][3][0] = t4;
                         t5 =  l4x4phm_t2[1][0][0]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t2[1][0][0] = t5;
                         t6 =  l4x4phm_t2[1][1][0]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t2[1][1][0] = t6;
                         t7 =  l4x4phm_t2[1][2][0]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t2[1][2][0] = t7;
                         t8 =  l4x4phm_t2[1][3][0]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t2[1][3][0] = t8;
                         t9 =  l4x4phm_t2[2][0][0]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t2[2][0][0] = t9;
                         t10 = l4x4phm_t2[2][1][0]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t2[2][1][0] = t10;
                         t11 = l4x4phm_t2[2][2][0]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t2[2][2][0] = t11;
                         t12 = l4x4phm_t2[2][3][0]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t2[2][3][0] = t12;
                         t13 = l4x4phm_t2[3][0][0]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t2[3][0][0] = t13;
                         t14 = l4x4phm_t2[3][1][0]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t2[3][1][0] = t14;
                         t15 = l4x4phm_t2[3][2][0]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t2[3][2][0] = t15;
                         t16 = l4x4phm_t2[3][3][0]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t2[3][3][0] = t16;
#endif
		       
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1				   
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][0]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][0] = t1;
			   }
		       }

#else
                        t1 = l4x4phm_t2[0][0][0]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t2[0][0][0] = t1;
                         t2 =  l4x4phm_t2[0][1][0]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t2[0][1][0] = t2;
                         t3 =  l4x4phm_t2[0][2][0]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t2[0][2][0] = t3;
                         t4 =  l4x4phm_t2[0][3][0]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t2[0][3][0] = t4;
                         t5 =  l4x4phm_t2[1][0][0]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t2[1][0][0] = t5;
                         t6 =  l4x4phm_t2[1][1][0]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t2[1][1][0] = t6;
                         t7 =  l4x4phm_t2[1][2][0]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t2[1][2][0] = t7;
                         t8 =  l4x4phm_t2[1][3][0]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t2[1][3][0] = t8;
                         t9 =  l4x4phm_t2[2][0][0]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t2[2][0][0] = t9;
                         t10 = l4x4phm_t2[2][1][0]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t2[2][1][0] = t10;
                         t11 = l4x4phm_t2[2][2][0]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t2[2][2][0] = t11;
                         t12 = l4x4phm_t2[2][3][0]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t2[2][3][0] = t12;
                         t13 = l4x4phm_t2[3][0][0]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t2[3][0][0] = t13;
                         t14 = l4x4phm_t2[3][1][0]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t2[3][1][0] = t14;
                         t15 = l4x4phm_t2[3][2][0]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t2[3][2][0] = t15;
                         t16 = l4x4phm_t2[3][3][0]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t2[3][3][0] = t16;
#endif
                     // Phase matrix: case 2
		  thinc = 3.141592653589793f-theta;
                  thsc  = 3.141592653589793f-theta;
                  phinc = 0.0f;
                  phsc  = 3.141592653589793f;
		  if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		       stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1	
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][1]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][1] = t1;
			   }
		       }

#else
                        t1 = l4x4phm_t2[0][0][1]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t2[0][0][1] = t1;
                        t2 =  l4x4phm_t2[0][1][1]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t2[0][1][1] = t2;
                        t3 =  l4x4phm_t2[0][2][1]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t2[0][2][1] = t3;
                        t4 =  l4x4phm_t2[0][3][1]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t2[0][3][1] = t4;
                        t5 =  l4x4phm_t2[1][0][1]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t2[1][0][1] = t5;
                        t6 =  l4x4phm_t2[1][1][1]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t2[1][1][1] = t6;
                        t7 =  l4x4phm_t2[1][2][1]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t2[1][2][1] = t7;
                        t8 =  l4x4phm_t2[1][3][1]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t2[1][3][1] = t8;
                        t9 =  l4x4phm_t2[2][0][1]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t2[2][0][1] = t9;
                        t10 = l4x4phm_t2[2][1][1]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t2[2][1][1] = t10;
                        t11 = l4x4phm_t2[2][2][1]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t2[2][2][1] = t11;
                        t12 = l4x4phm_t2[2][3][1]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t2[2][3][1] = t12;
                        t13 = l4x4phm_t2[3][0][1]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t2[3][0][1] = t13;
                        t14 = l4x4phm_t2[3][1][1]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t2[3][1][1] = t14;
                        t15 = l4x4phm_t2[3][2][1]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t2[3][2][1] = t15;
                        t16 = l4x4phm_t2[3][3][1]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t2[3][3][1] = t16;
#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1					   
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][1]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][1] = t1;
			   }
		       }
#else
                         t1 = l4x4phm_t2[0][0][1]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t2[0][0][1] = t1;
                        t2 =  l4x4phm_t2[0][1][1]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t2[0][1][1] = t2;
                        t3 =  l4x4phm_t2[0][2][1]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t2[0][2][1] = t3;
                        t4 =  l4x4phm_t2[0][3][1]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t2[0][3][1] = t4;
                        t5 =  l4x4phm_t2[1][0][1]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t2[1][0][1] = t5;
                        t6 =  l4x4phm_t2[1][1][1]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t2[1][1][1] = t6;
                        t7 =  l4x4phm_t2[1][2][1]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t2[1][2][1] = t7;
                        t8 =  l4x4phm_t2[1][3][1]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t2[1][3][1] = t8;
                        t9 =  l4x4phm_t2[2][0][1]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t2[2][0][1] = t9;
                        t10 = l4x4phm_t2[2][1][1]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t2[2][1][1] = t10;
                        t11 = l4x4phm_t2[2][2][1]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t2[2][2][1] = t11;
                        t12 = l4x4phm_t2[2][3][1]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t2[2][3][1] = t12;
                        t13 = l4x4phm_t2[3][0][1]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t2[3][0][1] = t13;
                        t14 = l4x4phm_t2[3][1][1]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t2[3][1][1] = t14;
                        t15 = l4x4phm_t2[3][2][1]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t2[3][2][1] = t15;
                        t16 = l4x4phm_t2[3][3][1]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t2[3][3][1] = t16;

#endif
		       
                    // Phase matrix: case 3
		    thinc = theta;
                    thsc  = theta;
                    phinc =  3.141592653589793f;
                    phsc  = 0.0f;
		    if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		       stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1			      
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][2]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][2] = t1;
			   }
		       }
#else
                        t1 = l4x4phm_t2[0][0][2]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t2[0][0][2] = t1;
                        t2 =  l4x4phm_t2[0][1][2]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t2[0][1][2] = t2;
                        t3 =  l4x4phm_t2[0][2][2]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t2[0][2][2] = t3;
                        t4 =  l4x4phm_t2[0][3][2]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t2[0][3][2] = t4;
                        t5 =  l4x4phm_t2[1][0][2]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t2[1][0][2] = t5;
                        t6 =  l4x4phm_t2[1][1][2]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t2[1][1][2] = t6;
                        t7 =  l4x4phm_t2[1][2][2]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t2[1][2][2] = t7;
                        t8 =  l4x4phm_t2[1][3][2]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t2[1][3][2] = t8;
                        t9 =  l4x4phm_t2[2][0][2]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t2[2][0][2] = t9;
                        t10 = l4x4phm_t2[2][1][2]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t2[2][1][2] = t10;
                        t11 = l4x4phm_t2[2][2][2]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t2[2][2][2] = t11;
                        t12 = l4x4phm_t2[2][3][2]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t2[2][3][2] = t12;
                        t13 = l4x4phm_t2[3][0][2]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t2[3][0][2] = t13;
                        t14 = l4x4phm_t2[3][1][2]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t2[3][1][2] = t14;
                        t15 = l4x4phm_t2[3][2][2]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t2[3][2][2] = t15;
                        t16 = l4x4phm_t2[3][3][2]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t2[3][3][2] = t16;
#endif
		       
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][2]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][2] = t1;
			   }
		       }

#else
                         t1 = l4x4phm_t2[0][0][2]+orient_distr*stokes4x4m[0][0];
                        l4x4phm_t2[0][0][2] = t1;
                        t2 =  l4x4phm_t2[0][1][2]+orient_distr*stokes4x4m[0][1];
                        l4x4phm_t2[0][1][2] = t2;
                        t3 =  l4x4phm_t2[0][2][2]+orient_distr*stokes4x4m[0][2];
                        l4x4phm_t2[0][2][2] = t3;
                        t4 =  l4x4phm_t2[0][3][2]+orient_distr*stokes4x4m[0][3];
                        l4x4phm_t2[0][3][2] = t4;
                        t5 =  l4x4phm_t2[1][0][2]+orient_distr*stokes4x4m[1][0];
                        l4x4phm_t2[1][0][2] = t5;
                        t6 =  l4x4phm_t2[1][1][2]+orient_distr*stokes4x4m[1][1];
                        l4x4phm_t2[1][1][2] = t6;
                        t7 =  l4x4phm_t2[1][2][2]+orient_distr*stokes4x4m[1][2];
                        l4x4phm_t2[1][2][2] = t7;
                        t8 =  l4x4phm_t2[1][3][2]+orient_distr*stokes4x4m[1][3];
                        l4x4phm_t2[1][3][2] = t8;
                        t9 =  l4x4phm_t2[2][0][2]+orient_distr*stokes4x4m[2][0];
                        l4x4phm_t2[2][0][2] = t9;
                        t10 = l4x4phm_t2[2][1][2]+orient_distr*stokes4x4m[2][1];
                        l4x4phm_t2[2][1][2] = t10;
                        t11 = l4x4phm_t2[2][2][2]+orient_distr*stokes4x4m[2][2];
                        l4x4phm_t2[2][2][2] = t11;
                        t12 = l4x4phm_t2[2][3][2]+orient_distr*stokes4x4m[2][3];
                        l4x4phm_t2[2][3][2] = t12;
                        t13 = l4x4phm_t2[3][0][2]+orient_distr*stokes4x4m[3][0];
                        l4x4phm_t2[3][0][2] = t13;
                        t14 = l4x4phm_t2[3][1][2]+orient_distr*stokes4x4m[3][1];
                        l4x4phm_t2[3][1][2] = t14;
                        t15 = l4x4phm_t2[3][2][2]+orient_distr*stokes4x4m[3][2];
                        l4x4phm_t2[3][2][2] = t15;
                        t16 = l4x4phm_t2[3][3][2]+orient_distr*stokes4x4m[3][3];
                        l4x4phm_t2[3][3][2] = t16;

#endif
                    // Phase matrix: case 4
                    thinc =  3.141592653589793f-theta;
                    thsc  = theta;
                    phinc = 0.0f;
                    phsc  =  3.141592653589793f;
		    if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		      }
		       stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][3]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][3] = t1;
			   }
		       }
#else
                         t1 = l4x4phm_t2[0][0][3]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t2[0][0][3] = t1;
                         t2 =  l4x4phm_t2[0][1][3]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t2[0][1][3] = t2;
                         t3 =  l4x4phm_t2[0][2][3]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t2[0][2][3] = t3;
                         t4 =  l4x4phm_t2[0][3][3]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t2[0][3][3] = t4;
                         t5 =  l4x4phm_t2[1][0][3]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t2[1][0][3] = t5;
                         t6 =  l4x4phm_t2[1][1][3]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t2[1][1][3] = t6;
                         t7 =  l4x4phm_t2[1][2][3]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t2[1][2][3] = t7;
                         t8 =  l4x4phm_t2[1][3][3]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t2[1][3][3] = t8;
                         t9 =  l4x4phm_t2[2][0][3]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t2[2][0][3] = t9;
                         t10 = l4x4phm_t2[2][1][3]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t2[2][1][3] = t10;
                         t11 = l4x4phm_t2[2][2][3]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t2[2][2][3] = t11;
                         t12 = l4x4phm_t2[2][3][3]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t2[2][3][3] = t12;
                         t13 = l4x4phm_t2[3][0][3]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t2[3][0][3] = t13;
                         t14 = l4x4phm_t2[3][1][3]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t2[3][1][3] = t14;
                         t15 = l4x4phm_t2[3][2][3]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t2[3][2][3] = t15;
                         t16 = l4x4phm_t2[3][3][3]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t2[3][3][3] = t16;
#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t2[k][l][3]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t2[k][l][3] = t1;
			   }
		       }

#else
                        t1 = l4x4phm_t2[0][0][3]+orient_distr*stokes4x4m[0][0];
                         l4x4phm_t2[0][0][3] = t1;
                         t2 =  l4x4phm_t2[0][1][3]+orient_distr*stokes4x4m[0][1];
                         l4x4phm_t2[0][1][3] = t2;
                         t3 =  l4x4phm_t2[0][2][3]+orient_distr*stokes4x4m[0][2];
                         l4x4phm_t2[0][2][3] = t3;
                         t4 =  l4x4phm_t2[0][3][3]+orient_distr*stokes4x4m[0][3];
                         l4x4phm_t2[0][3][3] = t4;
                         t5 =  l4x4phm_t2[1][0][3]+orient_distr*stokes4x4m[1][0];
                         l4x4phm_t2[1][0][3] = t5;
                         t6 =  l4x4phm_t2[1][1][3]+orient_distr*stokes4x4m[1][1];
                         l4x4phm_t2[1][1][3] = t6;
                         t7 =  l4x4phm_t2[1][2][3]+orient_distr*stokes4x4m[1][2];
                         l4x4phm_t2[1][2][3] = t7;
                         t8 =  l4x4phm_t2[1][3][3]+orient_distr*stokes4x4m[1][3];
                         l4x4phm_t2[1][3][3] = t8;
                         t9 =  l4x4phm_t2[2][0][3]+orient_distr*stokes4x4m[2][0];
                         l4x4phm_t2[2][0][3] = t9;
                         t10 = l4x4phm_t2[2][1][3]+orient_distr*stokes4x4m[2][1];
                         l4x4phm_t2[2][1][3] = t10;
                         t11 = l4x4phm_t2[2][2][3]+orient_distr*stokes4x4m[2][2];
                         l4x4phm_t2[2][2][3] = t11;
                         t12 = l4x4phm_t2[2][3][3]+orient_distr*stokes4x4m[2][3];
                         l4x4phm_t2[2][3][3] = t12;
                         t13 = l4x4phm_t2[3][0][3]+orient_distr*stokes4x4m[3][0];
                         l4x4phm_t2[3][0][3] = t13;
                         t14 = l4x4phm_t2[3][1][3]+orient_distr*stokes4x4m[3][1];
                         l4x4phm_t2[3][1][3] = t14;
                         t15 = l4x4phm_t2[3][2][3]+orient_distr*stokes4x4m[3][2];
                         l4x4phm_t2[3][2][3] = t15;
                         t16 = l4x4phm_t2[3][3][3]+orient_distr*stokes4x4m[3][3];
                         l4x4phm_t2[3][3][3] = t16;

#endif
                     // Extinction matrix: case 1
		    thinc = theta;
                    thsc  = thinc;
                    phinc = 3.141592653589793f;
                    phsc  = phinc;
		    if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		     }

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                  for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                        __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                       for(l=0; l != 2; ++l) {
                                           t1 = sm2x2avg_t2[k][l][0]+orient_distr*scat2x2m[k][l];
		                           sm2x2avg_t2[k][l][0] = t1;
		                        }
	                        }
#else
               t1 = sm2x2avg_t2[0][0][0]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t2[0][0][0] = t1;
	       t2 = sm2x2avg_t2[0][1][0]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t2[0][1][0] = t1;
	       t3 = sm2x2avg_t2[1][0][0]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t2[1][0][0] = t3;
	       t4 = sm2x2avg_t2[1][1][0]+orient_distr*scat2x2m[1][1];
               sm2x2avg_t2[1][1][0] = t4;
#endif		     

                         scat2x2m[1] = -scat2x2m[1];
		         scat2x2m[2] = -scat2x2m[2];
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                  for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                        __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                       for(l=0; l != 2; ++l) {
                                           t1 = sm2x2avg_t2[k][l][0]+orient_distr*scat2x2m[k][l];
		                           sm2x2avg_t2[k][l][0] = t1;
		                        }
	                        }
#else
               t1 = sm2x2avg_t2[0][0][0]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t2[0][0][0] = t1;
	       t2 = sm2x2avg_t2[0][1][0]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t2[0][1][0] = t1;
	       t3 = sm2x2avg_t2[1][0][0]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t2[1][0][0] = t3;
	       t4 = sm2x2avg_t2[1][1][0]+orient_distr*scat2x2m[1][1];
               sm2x2avg_t2[1][1][0] = t4;
#endif				 

                      // Extinction matrix: case 2
		        thinc = 3.141592653589793f-theta;
                        thsc  = 3.141592653589793f-theta;
                        phinc = 0.0f;
                        phsc  = phinc;
			if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		     }

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                    for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                        __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                       for(l=0; l != 2; ++l) {
                                           t1 = sm2x2avg_t2[k][l][1]+orient_distr*scat2x2m[k][l];
		                           sm2x2avg_t2[k][l][1] = t1;
		                       }
	                        }
#else
               t1 = sm2x2avg_t2[0][0][1]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t2[0][0][1] = t1;
	       t2 = sm2x2avg_t2[0][1][1]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t2[0][1][1] = t1;
	       t3 = sm2x2avg_t2[1][0][1]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t2[1][0][1] = t3;
               t4 = sm2x2avg_t2[1][1][1]+orient_distr*scat2x2m[1][1];
	       sm2x2avg_t2[1][1][1] = t4;
#endif		     
 
                         scat2x2m[1] = -scat2x2m[1];
		         scat2x2m[2] = -scat2x2m[2];
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                    for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                        __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                       for(l=0; l != 2; ++l) {
                                           t1 = sm2x2avg_t2[k][l][1]+orient_distr*scat2x2m[k][l];
		                           sm2x2avg_t2[k][l][1] = t1;
		                       }
	                        }
#else
               t1 = sm2x2avg_t2[0][0][1]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t2[0][0][1] = t1;
	       t2 = sm2x2avg_t2[0][1][1]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t2[0][1][1] = t1;
	       t3 = sm2x2avg_t2[1][0][1]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t2[1][0][1] = t3;
               t4 = sm2x2avg_t2[1][1][1]+orient_distr*scat2x2m[1][1];
	       sm2x2avg_t2[1][1][1] = t4;
#endif	
			 

		  } // end for (ii=0 loop
	       } // if(orient_distr block
	   } // end of for(jj=0 loop
       } // end nth2 != 0 block

       if((nth3!=0) && (nph3!=0)) {
           t1=0.0f;
           t2=0.0f;
           t3=0.0f;
           t4=0.0f;
           t5=0.0f;
           t6=0.0f;
           t7=0.0f;
           t8=0.0f;
           t9=0.0f;
           t10=0.0f;
           t11=0.0f;
           t12=0.0f;
           t13=0.0f;
           t14=0.0f;
           t15=0.0f;
           t16=0.0f;
	   for(jj=0; jj != nth3; ++jj) {
               thdr = tr_start3+dt_rad3*static_cast<float>(jj);
	       orient_distr = Compute_leaf_odf(leaf_orient,thdr);
	       if(orient_distr>0.0f) {
                  for(ii=0; ii != nph3; ++ii) {
                      phdr  = pr_start3+dp_rad3*static_cast<float>(ii);
		      thinc = theta;
                      thsc  = 3.141592653589793f-theta;
                      phinc = 3.141592653589793f;
                      phsc  = 0.0f;
		      if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		     }
		     stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1				  
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][0]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][0] = t1;
			   }
		       }

#else
                         t1 = l4x4phm_t3[0][0][0]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][0] = t1;
     t2 =  l4x4phm_t3[0][1][0]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][0] = t2;
     t3 =  l4x4phm_t3[0][2][0]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][0] = t3;
     t4 =  l4x4phm_t3[0][3][0]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][0] = t4;
     t5 =  l4x4phm_t3[1][0][0]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][0] = t5;
     t6 =  l4x4phm_t3[1][1][0]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][0] = t6;
     t7 =  l4x4phm_t3[1][2][0]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][0] = t7;
     t8 =  l4x4phm_t3[1][3][0]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][0] = t8;
     t9 =  l4x4phm_t3[2][0][0]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][0] = t9;
     t10 = l4x4phm_t3[2][1][0]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][0] = t10;
     t11 = l4x4phm_t3[2][2][0]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][0] = t11;
     t12 = l4x4phm_t3[2][3][0]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][0] = t12;
     t13 = l4x4phm_t3[3][0][0]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][0] = t13;
     t14 = l4x4phm_t3[3][1][0]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][0] = t14;
     t15 = l4x4phm_t3[3][2][0]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][0] = t15;
     t16 = l4x4phm_t3[3][3][0]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][0] = t16;

#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1					   
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][0]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][0] = t1;
			   }
		       }

#else
                                           t1 = l4x4phm_t3[0][0][0]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][0] = t1;
     t2 =  l4x4phm_t3[0][1][0]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][0] = t2;
     t3 =  l4x4phm_t3[0][2][0]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][0] = t3;
     t4 =  l4x4phm_t3[0][3][0]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][0] = t4;
     t5 =  l4x4phm_t3[1][0][0]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][0] = t5;
     t6 =  l4x4phm_t3[1][1][0]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][0] = t6;
     t7 =  l4x4phm_t3[1][2][0]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][0] = t7;
     t8 =  l4x4phm_t3[1][3][0]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][0] = t8;
     t9 =  l4x4phm_t3[2][0][0]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][0] = t9;
     t10 = l4x4phm_t3[2][1][0]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][0] = t10;
     t11 = l4x4phm_t3[2][2][0]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][0] = t11;
     t12 = l4x4phm_t3[2][3][0]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][0] = t12;
     t13 = l4x4phm_t3[3][0][0]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][0] = t13;
     t14 = l4x4phm_t3[3][1][0]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][0] = t14;
     t15 = l4x4phm_t3[3][2][0]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][0] = t15;
     t16 = l4x4phm_t3[3][3][0]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][0] = t16;  

#endif
                     // Phase matrix: case 2
		       thinc = 3.141592653589793f-theta;
                       thsc  = 3.141592653589793f-theta;
                       phinc = 0.0f;
                       thsc  = 3.141592653589793f;
		       if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
			                      
						 
		       }
		       else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                             rad_freq,rad_k0,rad_wv,lmg,lrho,
						     ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		       }
		     stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1				  
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][1]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][1] = t1;
			   }
		       }

#else
                        t1 = l4x4phm_t3[0][0][1]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][1] = t1;
     t2 =  l4x4phm_t3[0][1][1]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][1] = t2;
     t3 =  l4x4phm_t3[0][2][1]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][1] = t3;
     t4 =  l4x4phm_t3[0][3][1]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][1] = t4;
     t5 =  l4x4phm_t3[1][0][1]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][1] = t5;
     t6 =  l4x4phm_t3[1][1][1]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][1] = t6;
     t7 =  l4x4phm_t3[1][2][1]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][1] = t7;
     t8 =  l4x4phm_t3[1][3][1]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][1] = t8;
     t9 =  l4x4phm_t3[2][0][1]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][1] = t9;
     t10 = l4x4phm_t3[2][1][1]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][1] = t10;
     t11 = l4x4phm_t3[2][2][1]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][1] = t11;
     t12 = l4x4phm_t3[2][3][1]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][1] = t12;
     t13 = l4x4phm_t3[3][0][1]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][1] = t13;
     t14 = l4x4phm_t3[3][1][1]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][1] = t14;
     t15 = l4x4phm_t3[3][2][1]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][1] = t15;
     t16 = l4x4phm_t3[3][3][1]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][1] = t16;

#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1				   
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][1]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][1] = t1;
			   }
		       }

#else

                                            t1 = l4x4phm_t3[0][0][1]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][1] = t1;
     t2 =  l4x4phm_t3[0][1][1]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][1] = t2;
     t3 =  l4x4phm_t3[0][2][1]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][1] = t3;
     t4 =  l4x4phm_t3[0][3][1]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][1] = t4;
     t5 =  l4x4phm_t3[1][0][1]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][1] = t5;
     t6 =  l4x4phm_t3[1][1][1]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][1] = t6;
     t7 =  l4x4phm_t3[1][2][1]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][1] = t7;
     t8 =  l4x4phm_t3[1][3][1]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][1] = t8;
     t9 =  l4x4phm_t3[2][0][1]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][1] = t9;
     t10 = l4x4phm_t3[2][1][1]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][1] = t10;
     t11 = l4x4phm_t3[2][2][1]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][1] = t11;
     t12 = l4x4phm_t3[2][3][1]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][1] = t12;
     t13 = l4x4phm_t3[3][0][1]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][1] = t13;
     t14 = l4x4phm_t3[3][1][1]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][1] = t14;
     t15 = l4x4phm_t3[3][2][1]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][1] = t15;
     t16 = l4x4phm_t3[3][3][1]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][1] = t16;
#endif
		       // Phase matrix: case 3
		        thinc = theta;
                        thsc  = theta;
                        phinc = 3.141592653589793f;
                        phsc  = 0.0f;
			if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
			                      
						 
		       }
		       else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                             rad_freq,rad_k0,rad_wv,lmg,lrho,
						     ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		       }
		     stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1					  
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][2]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][2] = t1;
			   }
		       }

#else
                         t1 = l4x4phm_t3[0][0][2]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][2] = t1;
     t2 =  l4x4phm_t3[0][1][2]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][2] = t2;
     t3 =  l4x4phm_t3[0][2][2]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][2] = t3;
     t4 =  l4x4phm_t3[0][3][2]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][2] = t4;
     t5 =  l4x4phm_t3[1][0][2]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][2] = t5;
     t6 =  l4x4phm_t3[1][1][2]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][2] = t6;
     t7 =  l4x4phm_t3[1][2][2]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][2] = t7;
     t8 =  l4x4phm_t3[1][3][2]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][2] = t8;
     t9 =  l4x4phm_t3[2][0][2]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][2] = t9;
     t10 = l4x4phm_t3[2][1][2]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][2] = t10;
     t11 = l4x4phm_t3[2][2][2]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][2] = t11;
     t12 = l4x4phm_t3[2][3][2]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][2] = t12;
     t13 = l4x4phm_t3[3][0][2]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][2] = t13;
     t14 = l4x4phm_t3[3][1][2]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][2] = t14;
     t15 = l4x4phm_t3[3][2][2]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][2] = t15;
     t16 = l4x4phm_t3[3][3][2]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][2] = t16;

#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][2]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][2] = t1;
			   }
		       }
#else
                         t1 = l4x4phm_t3[0][0][2]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][2] = t1;
     t2 =  l4x4phm_t3[0][1][2]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][2] = t2;
     t3 =  l4x4phm_t3[0][2][2]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][2] = t3;
     t4 =  l4x4phm_t3[0][3][2]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][2] = t4;
     t5 =  l4x4phm_t3[1][0][2]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][2] = t5;
     t6 =  l4x4phm_t3[1][1][2]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][2] = t6;
     t7 =  l4x4phm_t3[1][2][2]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][2] = t7;
     t8 =  l4x4phm_t3[1][3][2]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][2] = t8;
     t9 =  l4x4phm_t3[2][0][2]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][2] = t9;
     t10 = l4x4phm_t3[2][1][2]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][2] = t10;
     t11 = l4x4phm_t3[2][2][2]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][2] = t11;
     t12 = l4x4phm_t3[2][3][2]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][2] = t12;
     t13 = l4x4phm_t3[3][0][2]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][2] = t13;
     t14 = l4x4phm_t3[3][1][2]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][2] = t14;
     t15 = l4x4phm_t3[3][2][2]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][2] = t15;
     t16 = l4x4phm_t3[3][3][2]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][2] = t16;

#endif
		       
                     // Phase matrix: case 4
		       thinc = 3.141592653589793f-theta;
                       thsc  = theta;
                       phinc = 0.0f;
                       phsc  = 3.141592653589793f;
		       if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
			                      
						 
		       }
		       else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                             rad_freq,rad_k0,rad_wv,lmg,lrho,
						     ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		       }
		       stokes_matrix(&scat2x2m[0],
		                  &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1				  
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][3]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][3] = t1;
			   }
		       }

#else

                       t1 = l4x4phm_t3[0][0][3]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][3] = t1;
     t2 =  l4x4phm_t3[0][1][3]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][3] = t2;
     t3 =  l4x4phm_t3[0][2][3]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][3] = t3;
     t4 =  l4x4phm_t3[0][3][3]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][3] = t4;
     t5 =  l4x4phm_t3[1][0][3]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][3] = t5;
     t6 =  l4x4phm_t3[1][1][3]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][3] = t6;
     t7 =  l4x4phm_t3[1][2][3]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][3] = t7;
     t8 =  l4x4phm_t3[1][3][3]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][3] = t8;
     t9 =  l4x4phm_t3[2][0][3]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][3] = t9;
     t10 = l4x4phm_t3[2][1][3]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][3] = t10;
     t11 = l4x4phm_t3[2][2][3]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][3] = t11;
     t12 = l4x4phm_t3[2][3][3]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][3] = t12;
     t13 = l4x4phm_t3[3][0][3]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][3] = t13;
     t14 = l4x4phm_t3[3][1][3]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][3] = t14;
     t15 = l4x4phm_t3[3][2][3]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][3] = t15;
     t16 = l4x4phm_t3[3][3][3]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][3] = t16;
#endif
                     scat2x2m[1] = -scat2x2m[1];
		     scat2x2m[2] = -scat2x2m[2];
		     stokes_matrix(&scat2x2m[0],
		                   &stokes4x4m[0]);

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1				   
                       for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
        __assume_aligned(stokes4x4m,64);
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        stokes4x4m = (float*)__builtin_assume_aligned(stokes4x4m,64);
#pragma omp simd
#endif
                           for(l=0; l != 4; ++l) {
                               t1 = l4x4phm_t3[k][l][3]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t3[k][l][3] = t1;
			   }
		       }

#else

                                         t1 = l4x4phm_t3[0][0][3]+orient_distr*stokes4x4m[0][0];
     l4x4phm_t3[0][0][3] = t1;
     t2 =  l4x4phm_t3[0][1][3]+orient_distr*stokes4x4m[0][1];
     l4x4phm_t3[0][1][3] = t2;
     t3 =  l4x4phm_t3[0][2][3]+orient_distr*stokes4x4m[0][2];
     l4x4phm_t3[0][2][3] = t3;
     t4 =  l4x4phm_t3[0][3][3]+orient_distr*stokes4x4m[0][3];
     l4x4phm_t3[0][3][3] = t4;
     t5 =  l4x4phm_t3[1][0][3]+orient_distr*stokes4x4m[1][0];
     l4x4phm_t3[1][0][3] = t5;
     t6 =  l4x4phm_t3[1][1][3]+orient_distr*stokes4x4m[1][1];
     l4x4phm_t3[1][1][3] = t6;
     t7 =  l4x4phm_t3[1][2][3]+orient_distr*stokes4x4m[1][2];
     l4x4phm_t3[1][2][3] = t7;
     t8 =  l4x4phm_t3[1][3][3]+orient_distr*stokes4x4m[1][3];
     l4x4phm_t3[1][3][3] = t8;
     t9 =  l4x4phm_t3[2][0][3]+orient_distr*stokes4x4m[2][0];
     l4x4phm_t3[2][0][3] = t9;
     t10 = l4x4phm_t3[2][1][3]+orient_distr*stokes4x4m[2][1];
     l4x4phm_t3[2][1][3] = t10;
     t11 = l4x4phm_t3[2][2][3]+orient_distr*stokes4x4m[2][2];
     l4x4phm_t3[2][2][3] = t11;
     t12 = l4x4phm_t3[2][3][3]+orient_distr*stokes4x4m[2][3];
     l4x4phm_t3[2][3][3] = t12;
     t13 = l4x4phm_t3[3][0][3]+orient_distr*stokes4x4m[3][0];
     l4x4phm_t3[3][0][3] = t13;
     t14 = l4x4phm_t3[3][1][3]+orient_distr*stokes4x4m[3][1];
     l4x4phm_t3[3][1][3] = t14;
     t15 = l4x4phm_t3[3][2][3]+orient_distr*stokes4x4m[3][2];
     l4x4phm_t3[3][2][3] = t15;
     t16 = l4x4phm_t3[3][3][3]+orient_distr*stokes4x4m[3][3];
     l4x4phm_t3[3][3][3] = t16;
#endif
                     // Extinction matrix: case 1
		       thinc = theta;
                       thsc  = thinc;
                       phinc = 3.141592653589793f;
                       phsc  = phinc;
		       if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
			                      
						 
		       }
		       else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                             rad_freq,rad_k0,rad_wv,lmg,lrho,
						     ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		       }

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                          for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                               __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                               scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                                for(l=0; l != 2; ++l) {
                                                     t1 = sm2x2avg_t3[k][l][0]+orient_distr*scat2x2m[k][l];
		                                     sm2x2avg_t3[k][l][0] = t1;
		                               }
	                                }
#else
               t1 = sm2x2avg_t3[0][0][0]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t3[0][0][0] = t1;
	       t2 = sm2x2avg_t3[0][1][0]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t3[0][1][0] = t1;
	       t3 = sm2x2avg_t3[1][0][0]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t3[1][0][0] = t3;
	       t4 = sm2x2avg_t3[1][1][0]+orient_distr*scat2x2m[1][1];
               sm2x2avg_t3[1][1][0] = t4;
#endif		       

                           scat2x2m[1] = -scat2x2m[1];
		           scat2x2m[2] = -scat2x2m[2];

#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                          for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                               __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                               scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                                for(l=0; l != 2; ++l) {
                                                     t1 = sm2x2avg_t3[k][l][0]+orient_distr*scat2x2m[k][l];
		                                     sm2x2avg_t3[k][l][0] = t1;
		                               }
	                                }
#else
               t1 = sm2x2avg_t3[0][0][0]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t3[0][0][0] = t1;
	       t2 = sm2x2avg_t3[0][1][0]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t3[0][1][0] = t1;
	       t3 = sm2x2avg_t3[1][0][0]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t3[1][0][0] = t3;
	       t4 = sm2x2avg_t3[1][1][0]+orient_distr*scat2x2m[1][1];
               sm2x2avg_t3[1][1][0] = t4;
#endif				   

                      // Extinction matrix: case 2
		        thinc =  3.141592653589793f-theta;
                        thsc  =  3.141592653589793f-theta;
                        phinc = 0.0f;
                        phsc  = phinc;
			if(po) {
                           Leaf_PO_approximation(thinc,phinc,thsc,phsc,thdr,phdr,
			                         rad_freq,rad_k0,rad_wv,lmg,lrho,
						 ldens,ldiam,lthick,epsr,&scat2x2m[0]);
			                      
						 
		       }
		       else {
                            Leaf_Rayleigh_scattering(thinc,phinc,thsc,phsc,thdr,phdr,
			                             rad_freq,rad_k0,rad_wv,lmg,lrho,
						     ldens,ldiam,lthick,epsr,&scat2x2m[0]);
		       }


#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                  for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                      __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                          for(l=0; l != 2; ++l) {
                                              t1 = sm2x2avg_t3[k][l][1]+orient_distr*scat2x2m[k][l];
		                              sm2x2avg_t3[k][l][1] = t1;
		                          }
	                           }
#else
               t1 = sm2x2avg_t3[0][0][1]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t3[0][0][1] = t1;
	       t2 = sm2x2avg_t3[0][1][1]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t3[0][1][1] = t1;
	       t3 = sm2x2avg_t3[1][0][1]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t3[1][0][1] = t3;
               t4 = sm2x2avg_t3[1][1][1]+orient_distr*scat2x2m[1][1];
	       sm2x2avg_t3[1][1][1] = t4;
#endif		       

                           scat2x2m[1] = -scat2x2m[1];
		           scat2x2m[2] = -scat2x2m[2];
#if (LEAF_PHASE_MATRICES_AUTOVECTORIZE) == 1

                                  for(k=0; k != 2; ++k) {
#if defined __INTEL_COMPILER
                                      __assume_aligned(scat2x2m,64);
#pragma vector always
#pragma code_align(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      scat2x2m = (std::complex<float>*)__builtin_assume_aligned(scat2x2m,64);
#pragma omp simd
#endif
                                          for(l=0; l != 2; ++l) {
                                              t1 = sm2x2avg_t3[k][l][1]+orient_distr*scat2x2m[k][l];
		                              sm2x2avg_t3[k][l][1] = t1;
		                          }
	                           }
#else
               t1 = sm2x2avg_t3[0][0][1]+orient_distr*scat2x2m[0][0];
	       sm2x2avg_t3[0][0][1] = t1;
	       t2 = sm2x2avg_t3[0][1][1]+orient_distr*scat2x2m[0][1];
	       sm2x2avg_t3[0][1][1] = t1;
	       t3 = sm2x2avg_t3[1][0][1]+orient_distr*scat2x2m[1][0];
	       sm2x2avg_t3[1][0][1] = t3;
               t4 = sm2x2avg_t3[1][1][1]+orient_distr*scat2x2m[1][1];
	       sm2x2avg_t3[1][1][1] = t4;
#endif	
			   

		  } // end of for(ii=0 loop
	       } // if(orient_distr block
	   } // end of for(jj=0 loop
       } // end nth3 != 0 block

        //! Phase and M matrices
        work = 1.0f/norm;
	for(i=0; i != 4; ++i) {
#if defined __INTEL_COMPILER
#pragma unroll_and_jam (4)
#endif
             for(k=0; k != 4; ++k) {
#if defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd
#endif
                     for(l=0; l != 4; ++l) {
		          const float t0 = dp1t1*l4x4phm_t1[i][k][l];
			  const float t1 = dp2t2*l4x4phm_t2[i][k][l];
			  const float t2 = dp3t3*l4x4phm_t3[i][k][l];
                          l4x4phm[Ix3D(i,4,k,4,l)] = work*(t0+t1+t2);
		     }
	     }
	}
	// Loops fused, may inhibit LSD detection.
	cwork = j*6.283185307179586f/(norm*rad_k0);
	for(k=0; k != 2; ++k) {
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd
#endif
                for(l=0; l != 2; ++l) {
                    const std::complex<float> t0 = dp1t1*sm2x2avg_t1[k][l][0];
		    const std::complex<float> t1 = dp2t2*sm2x2avg_t2[k][l][0];
		    const std::complex<float> t2 = dp3t3*sm2x2avg_t3[k][l][0];
		    l2x2mp[Ix2D(k,2,l)] = cwork*(t0+t1+t2);
		    const std::complex<float> t3 = dp1t1*sm2x2avg_t1[k][l][1];
		    const std::complex<float> t4 = dp2t2*sm2x2avg_t2[k][l][1];
		    const std::complex<float> t5 = dp3t3*sm2x2avg_t3[k][l][1];
		    l2x2mn[Ix2D(k,2,l)] = cwork*(t3+t4+t5);
		}
	}
}


float
gms::math
::Compute_leaf_odf(const int32_t leaf_orient,
		   const float th) {

     
     float mu,nu,pi_tol
     const float tol = 8.727e-7f;
     pi_tol = 1.570796326794897f-tol;
     switch(leaf_orient) {

          case 1: {
                     if(th>1.570796326794897f)
		        return (0.0f);
		     else if(th==1.570796326794897f)
		        return (0.5f);
		     else
		        return (sinf(th));
	   	 }
		  break;
	   case 2: {
                     if(th>pi_tol)
		        return (0.0f);
		     else if(th<=tol)
		        return (0.0f);
		     else {
		        mu = 2.770f;
			nu = 1.172f;
		        return (Leaf_ang_orientation(mu,nu,th,0.0f));
		      }
	         }
		   break;
	    case 3: {
                      if(th>=pi_tol)
		         return (0.0f);
		      else if(th<=tol)
		         return (0.0f);
		      else {
                         nu = 1.1720f;
			 mu = 2.7700f;
			 return (Leaf_ang_orientation(mu,nu,th,0.0f));
		      }
	         }
                   break;
	     case 4: {
                       if(th>=pi_tol)
		          return (0.0f);
		       else if(th<=tol)
		          return (0.0f);
		       else {
                          nu = 3.326f;
			  mu = 3.326f;
			  return (Leaf_ang_orientation(mu,nu,th,0.0f));
		       }
	         }
		   break;
	      case 5: {
                        if(th>=pi_tol)
		          return (0.0f);
		        else if(th<=tol)
		          return (0.0f);
		        else {
                          nu = 0.433f;
			  mu = 0.433f;
			  return (Leaf_ang_orientation(mu,nu,th,0.0f));
		       }
	         }
		   break;
	      case 6: {
                        if(th>=pi_tol)
			   return (0.0f);
			else if(th<=tol)
			   return (0.0f);
			else {
                          nu = 1.0f;
			  mu = 1.0f;
			  return (Leaf_ang_orientation(mu,nu,th,0.0f));
			}
	         }
		   break;
	       case 7: {
                         if(th>=pi_tol)
			    return (0.0f);
			 else if(th<=tol)
			    return (0.0f);
			 else {
                           nu = 1.101f;
			   mu = 1.930f;
			   return (Leaf_ang_orientation(mu,nu,th,0.0f));
			 }
	         }
		   break;
	       default : {
                             if(leaf_orient<1 || leaf_orient>7){
                                 return (std::numeric_limits<float>::quiet_NaN());
                               }
	         }
        }
}

void
gms::math::
Leaf_PO_approximation(const float thinc,
		      const float phinc,
		      const float thsc,
		      const float phsc,
		      const float thdr,
		      const float phdr,
		      const float rad_freq,
		      const float rad_k0,
		      const float rad_wv,
		      const float leaf_mg,
		      const float leaf_rho,
		      const float leaf_dens,
		      const float leaf_diam,
		      const float leaf_tau,
		      const std::complex<float> epsr,
		      std::complex<float> * __restrict __ATTR_ALIGN__(32) scat_mat) {

           std::complex<float> j,eps,res,constant;
	   std::complex<float> gamh,game,gamhe_c1,gamhe_c2;
	   float tau,q,p,a,b,tol;
	   float cosphi1;
	   float thj,phij;
	   float sin_thi,cos_thi,sin_thj,cos_thj;
	   float sin_ths,cos_ths,cos_phj,sin_phij,sin_phji;
	   float cos_phij,cos_phji,sin_phsj,sin_phjs,cos_phjs;
	   float cos_beta,sin_beta,sin_phi,cos_phi,cos_beta_phi;
	   float sin_phpr,cos_phpr,sin_betapr;
	   float u,v,sinu,sinv,sinu_u,sinv_v;
	   float cosb_p,sintj_ti,sints_tj;
	   float costj_sinpjs,costs_sinpsj;
	   float sinpij_cospsj,cpij_cti_ctj,cpsj_cts_ctj;
	   float  s1,s2,s3,s4,s5,w1,w2,t0,t1,t2,t3,t4;
           tol = 0.0001f;
           t0  = sinf(thdr);
           j = {0.0f,1.0f};
           t1 = sinf(thinc);
           tau = leaf_tau/100.0f;
           t2 = cosf(phdr-phinc);
           eps = epsr;
           t3 = cosf(thinc);
           a = 0.5f*1.772453850905516f*leaf_diam/100.0f;
           t4 = cosf(thdr);
           b = a;
           cosphi1 = -(t0*t1*t2*t3*t4);
           if(cosphi1<0.0f) {
              thj = 3.141592653589793f-thinc;
              phij = 3.141592653589793f+phinc;
              cosphi1 = -cosphi1;
	   }
           else {
              thj = thinc;
              phij = phinc;
           }
	   sin_thi = t1;
           cos_thi = t3;
           sin_thj = sinf(thj);
           cos_thj = cosf(thj);
           sin_ths = sinf(thsc);
           cos_ths = cosf(thsc);
           cos_phij = cosf(phij);
           sin_phij = sinf(phinc-phij);
           sin_phji = -sin_phij;
           w1 = sin_thi*sin_phji;
           cos_phij = cosf(phinc-phij);
           q = 1.0f/sqrt(1.0f-w1*w1);
           cos_phji = cos_phij;
           cos_beta = q*cosphi1;
           sin_phsj = sinf(phsc-phij);
           sin_beta = q*(-cos_thj*sin_thi*cos_phji+
                         cos_thi*sin_thj); 
           sin_phjs = -sin_phsj;
           sin_phi  = sin_thi*sin_phij;
           cos_phsj = cosf(phsc-phij);
           cos_phi  = sqrtf(1.0f-sin_phi*sin_phi);
           cos_beta_phi = cos_beta*cos_phi;
           sin_phpr = sin_ths*sin_phsj;
           cos_phpr = sqrtf(1.0f-sin_phpr*sin_phpr);
           sin_betapr = (cos_ths*sin_thi- 
                         cos_thj*sin_ths*cos_phsj)/cos_phpr;
           p = 1.0f/sqrtf(1.0f-cos_beta_phi*cos_beta_phi);
           res = j/(rad_k0*tau*(eps-1.0f));
           gamh = 1.0f/(1.0f+2.0f*res/cosphi1);
           game = 1.0f/(1.0f+2.0f*res*cosphi1);
           u = 0.5f*rad_k0*a*(sin_phi-sin_phpr);
           sinu = sinf(u);
           v = 0.5f*rad_k0*b*(sin_beta*cos_phi-sin_betapr*cos_phpr);
           sinv = sinf(v);
	   if(fabsf(sinu)<=tol){
              sinu_u = 1.0f
	   }
           else {
              sinu_u = sinu/u;
           }
           cosb_p = cos_beta*cos_phi;
           gamhe_c1 = (gamh-game)*cosb_p;
           gamhe_c2 = gamh-cosb_p*cosb_p*game;
           !
           if(abs(sinv)<=tol){
              sinv_v = 1.0f;
	   }
           else{
              sinv_v = sinv/v;
           }
	   constant = -j*rad_k0*a*b*sinv_v*sinu_u*p*p/6.283185307179586f;
           sintj_ti = sin_thj*sin_thi;
           sints_tj = sin_ths*sin_thj;
           costj_sinpjs = cos_thj*sin_phjs;
           s5 = sin_phij*costj_sinpjs;
           costs_sinpsj = cos_ths*sin_phsj;
           s3 = sin_phij*costs_sinpsj;
           sinpij_cospsj = sin_phij*cos_phsj;
           cpij_cti_ctj  = cos_phij*cos_thi*cos_thj;
           s1 = sintj_ti+cpij_cti_ctj;
           cpsj_cts_ctj  = cos_phsj*cos_ths*cos_thj;
           s2 = sints_tj+cpsj_cts_ctj;
           w1 = s1*s2+cos_thi*s3;
           s4 = sin_phij*s2;
           w2 = cos_phij*s2+cos_thj*s3;
#if defined __INTEL_COMPILER
           __assume_aligned(scat_mat,32);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
           scat_mat = (std::complex<float>*)__builtin_assume_aligned(scat_mat,32);
#endif
           scat_mat[0] = constant*(w1*gamhe_c1+w2*gamhe_c2);
           w1 = -cos_thj*s4+cos_phij*costs_sinpsj;
           w2 = -cos_thi*s4+s1*costs_sinpsj
           scat_mat[1] = constant*(w1*gamhe_c1+w2*gamhe_c2);
           w1 = s1-costj_sinpjs+cos_thi*sinpij_cospsj;
           w2 = cos_phij*costj_sinpjs+cos_thj*sinpij_cospsj;
           scat_mat[2] = constant*(w1*gamhe_c1+w2*gamhe_c2);
           w1 = -cos_phij*s5+cos_phij*cos_phsj;
           w2 = -cos_thi*s5+s1*cos_phsj;
           scat_mat[3] = constant*(w1+gamhe_c1+w2*gamhe_c2);
}

void
gms::math::
Leaf_Rayleigh_scattering( const float thinc,
		          const float phinc,
		          const float thsc,
		          const float phsc,
		          const float thdr,
		          const float phdr,
		          const float rad_freq,
		          const float rad_k0,
		          const float rad_wv,
		          const float leaf_mg,
		          const float leaf_rho,
		          const float leaf_dens,
		          const float leaf_diam,
		          const float leaf_tau,
		          const std::complex<float> epsr,
		          std::complex<float> * __restrict __ATTR_ALIGN__(32) scat_mat) {

     __attribute__((aligned(16))) float xhat[4];
     __attribute__((aligned(16))) float yhat[4];
     __attribute__((aligned(16))) float zhat[4];
     __attribute__((aligned(16))) float xhatl[4];
     __attribute__((aligned(16))) float yhatl[4];
     __attribute__((aligned(16))) float zhatl[4];
     __attribute__((aligned(16))) float khati[4];
     __attribute__((aligned(16))) float khats[4];
     __attribute__((aligned(16))) float hhati[4];
     __attribute__((aligned(16))) float vhati[4];
     __attribute__((aligned(16))) float hhats[4];
     __attribute__((aligned(16))) float vhats[4];
     std::complex<float> eps,cdum;
     std::complex<float> Vd,cduma,cdumb,cdumc;
     float dum,dum2,sumcheck,t_lf, d_lf;
     float vhsdxhl,vhsdyhl,vhsdzhl,hhsdxhl,hhsdyhl,hhsdzhl;
     float vhidxhl,vhidyhl,vhidzhl,hhidxhl,hhidyhl,hhidzhl;
     float Ae,Be,Ce,Ac,Ab,Aa,Vo;
     float t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
     xhat[4] = {};
     _mm_store_ps(&xhat[0],_mm_setr_ps(0.0f,0.0f,0.0f,1.0f));
     t0 = sinf(thinc);
     t1 = cosf(phinc);
     yhat[4] = {};
     khati[0] = t0*t1;
     _mm_store_ps(&yhat[0],_mm_setr_ps(0.0f,1.0f,0.0f,0.0f));
     t6 = cosf(thdr);
     t7 = cosf(phdr);
     t8 = sinf(thdr);
     khati[1] = t0*t1;
     zhat[4] = {};
     _mm_store_ps(&zhat[0],_mm_setr_ps(0.0f,0.0f,1.0f,0.0f));
     t9 = sinf(phdr);
     khati[2] = t1;
     khati[3] = 1.0f;
     t2 = sinf(thsc);
     t3 = cosf(phsc);
     vec1x3_smooth(&khati[0],&khati[0]);
     t4 = sinf(phsc);
     khats[4] = {};
     t5 = sinf(phinc);
     _mm_store_ps(&khats[0],_mm_setr_ps(1.0f,cosf(thsc),t2*t4,t2*t3));
     vec1x3_smooth(&khats[0],&khats[0]);
     hhats[4] = {};
     _mm_store_ps(&hhati[0],_mm_setr_ps(1.0f,0.0f,t1,-t5));
     vec1x3_smooth(&hhati[0],&hhati[0]);
     _mm_store_ps(&hhats[0],_mm_setr_ps(-t4,t3,0.0f,1.0f));
     vhati[4] = {};
     cross_prod(&hhati[0],&khati[0],&vhati[0]);
     vec1x3_smooth(&vhati[0],&vhati[0]);
     vec1x3_smooth(&hhats[0],&hhats[0]);
     vhats[4] = {};
     cross_prod(&hhats[0],&khats[0],&vhats[0]);
     vec1x3_smooth(&vhats[0],&vhats[0]);
     _mm_store_ps(&xhatl[0],_mm_setr_ps(t6*t7,t6*t8,-t8,1.0f));
     vec1x3_smooth(&xhatl[0],&xhatl[0]);
     _mm_store_ps(&yhatl[0],_mm_setr_ps(-t8,t7,0.0f,1.0f));
     vec1x3_smooth(&yhatl[0],&yhatl[0]);
     _mm_store_ps(&zhatl[0],_mm_setr_ps(t8*t7,t8*t9,t6,1.0f));
     vec1x3_smooth(&zhatl[0],&zhatl[0]);
     vhsdxhl = 0.0f;
     t_lf = leaf_tau;
     vhsdxhl = dot_prod(vhats,xhatl);
     vhsdyhl = 0.0f;
     d_lf = leaf_diam;
     vhsdyhl = dot_prod(vhats,yhatl);
     vhsdzhl = 0.0f;
     eps = epsrc;
     vhsdzhl = dot_prod(vhats,zhatl);
     dum = powf(1.5f,0.33333333333333f);
     Ce  = (t_lf/200.0f)*dum;
     Ae  = (d_lf/200.0f)*dum;
     Be  = Ae;
     hhsdxhl = 0.0f;
     hhsdxhl = dot_prod(hhats,xhatl);
     dum = Ae*Ae-Ce*Ce;
     hhsdyhl = 0.0f;
     hhsdyhl = dot_prod(hhats,yhatl);
     dum2 = (sqrtf(dum))/Ce;
     Ac =  2.0f/(powf(dum,1.5f))*(dum2-atan(dum2));
     hhsdzhl = 0.0f;
     hhsdzhl = dot_prod(hhats,zhatl);
     Ab =  (2.0f/(Ae*Be*Ce)-Ac)*0.5f;
     Aa = Ab;
     vhidxhl = 0.0f;
     vhidxhl = dot_prod(vhati,xhatl);
     Vo = 12.566370614359173f*Ae*Be*Ce*0.3333333333333f;
     vhidyhl = 0.0f;
     vhidyhl = dot_prod(vhati,yhatl);
     Vd = (Ae*Be*Ce*0.5f)*(eps-1.0f);
     vhidzhl = 0.0f;
     vhidzhl = dot_prod(vhati,zhatl);
     hhidxhl = 0.0f;
     hhidxhl = dot_prod(hhati,xhatl);
     cdum = (rad_k0*rad_k0/12.566370614359173f)*Vo*(eps-1.0f);
     hhidyhl = 0.0f;
     hhidyhl = dot_prod(hhati,yhatl);
     cduma = 1.0f + Vd*Aa;
     cdumb = cduma;
     hhidzhl = 0.0f;
     hhidzhl = dot_prod(hhati,zhatl);
     cdumc = 1.0f + Vd*Ac;
#if defined __INTEL_COMPILER
     __assume_aligned(scat_mat,32);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     scat_mat = (std::complex<float>*)__builtin_assume_aligned(scat_mat,32);
#endif
     scat_mat[0] = cdum*(vhsdxhl*vhidxhl/cduma + 
		   vhsdyhl*vhidyhl/cdumb + vhsdzhl*vhidzhl/cdumc);

     scat_mat[1] = cdum*(hhsdxhl*vhidxhl/cduma + 
                   hhsdyhl*vhidyhl/cdumb + hhsdzhl*vhidzhl/cdumc);

     scat_mat[2] = cdum*(vhsdxhl*hhidxhl/cduma + 
                   vhsdyhl*hhidyhl/cdumb + vhsdzhl*hhidzhl/cdumc);

     scat_mat[3] = cdum*(hhsdxhl*hhidxhl/cduma + 
                   hhsdyhl*hhidyhl/cdumb + hhsdzhl*hhidzhl/cdumc);
}

void
gms::math
::Cylinder_scattering(const float thetai,
		      const float phii,
		      const float thetas,
		      const float phis,
		      const float thetac,
		      const float phic,
		      const float diam,
		      const float length,
		      const float rad_k0,
		      const std::complex<float> epsr,
		      std::complex<float> * __restrict __ATTR_ALIGN__(32) scat_mat) {
     __attribute__((aligned(24))) std::complex<float> pdoth[3] = {};
     __attribute__((aligned(24))) std::complex<float> pdotv[3] = {};
     __attribute__((aligned(24))) std::complex<float> Sh[3] = {};
     __attribute__((aligned(24))) std::complex<float> Sv[3] = {};
     __attribute__((aligned(16))) float xhat[4]   = {};
     __attribute__((aligned(16))) float yhat[4]   = {};
     __attribute__((aligned(16))) float zhat[4]   = {};
     __attribute__((aligned(16))) float xhatl[4]  = {};
     __attribute__((aligned(16))) float yhatl[4]  = {};
     __attribute__((aligned(16))) float zhatl[4]  = {};
     __attribute__((aligned(16))) float khati[4]  = {};
     __attribute__((aligned(16))) float khats[4]  = {};
     __attribute__((aligned(16))) float hhati[4]  = {};
     __attribute__((aligned(16))) float vhati[4]  = {};
     __attribute__((aligned(16))) float hhats[4]  = {};
     __attribute__((aligned(16))) float vhats[4]  = {};
     __attribute__((aligned(16))) float khatip[4] = {};
     __attribute__((aligned(16))) float khatsp[4] = {};
     __attribute__((aligned(16))) float vhatip[4] = {};
     __attribute__((aligned(16))) float hhatip[4] = {};
     __attribute__((aligned(16))) float vhatsp[4] = {};
     __attribute__((aligned(16))) float hhatsp[4] = {}; 
     __attribute__((aligned(16))) float xhatp[4]  = {};
     __attribute__((aligned(16))) float yhatp[4]  = {};
     __attribute__((aligned(16))) float zhatp[4]  = {};
     std::complex<float> Pxx,Pyy,Pzz,Svv,Svh,Shh;
     float area,U,sfact,mfact,khsdzhl,cosbeta,radius;
     float t0,t1,t2,t3,t4,t5,t6,t7;
     float t8,t9,t10,t11,t12,t13,t14;
     radius = diam*0.5f;
     _mm_store_ps(&xhat[0],_mm_setr_ps(0.0f,0.0f,1.0f,1.0f)); // 4th element is unused
     t0 = sinf(thetai);
     t1 = cosf(phii);
     khati[0] = t0*t1;
     _mm_store_ps(&yhat[0],_mm_setr_ps(0.0f,1.0f,0.0f,1.0f));
     t4 = sinf(phii);
     khati[1] = t0*t4;
     _mm_store_ps(&zhat[0],_mm_store_ps(0.0f,0.0f,1.0f,1.0f));
     khati[2] = cosf(thetai);
     khati[3] = 1.0f;
     vec1x3_smooth(&khati[0],&khati[0]);
     t2 = sinf(thetas);
     t3 = cosf(phis);
     khats[0] = t2*t3;
     khats[1] = t2*sinf(phis);
     khats[2] = cosf(thetas);
     khats[3] = 1.0f;
     vec1x3_smooth(&khats[0],&khats[0]);
     hhati[0] = -t4;
     hhati[1] = t1;
     hhati[2] = 0.0f;
     hhati[3] = 1.0f;
     vec1x3_smooth(&hhati[0],&hhati[0]);
     hhats[0] = -sinf(phiis);
     hhats[1] = t3;
     cross_prod(&hhati[0],&khati[0],&vhati[0]);
     vec1x3_smooth(&vhati[0],&vhati[0]);
     hhats[2] = 0.0f;
     hhats[3] = 1.0f;
     vec1x3_smooth(&hhats[0],&hhats[0]);
     t5 = cosf(thetac);
     t6 = cosf(phic);
     cross_prod(&hhats[0],&khats[0],&vhats[0]);
     t7 = sinf(thetac);
     t8 = sinf(phic);
     vec1x3_smooth(&vhats[0],&vhats[0]);
     _mm_store_ps(&xhatl[0],_mm_setr_ps(1.0f,-t7,t5*t8,t5*t6));
     vec1x3_smooth(&xhatl[0],&xhatl[0]);
     _mm_store_ps(&yhatl[0],_mm_setr_ps(1.0f,0.0f,t6,-t8));
     vec1x3_smooth(&yhatl[0],&yhatl[0]);
     t9 = sinf(thetac);
     _mm_store_ps(&zhatl[0],_mm_setr_ps(1.0f,t5,t9*t8,t9*t6));
     vec1x3_smooth(&zhatl[0],&zhatl[0]);
     xhatp[0] = xhatl[0];
     vhatip[0] = dot_prod(vhati,xhatl);
     vhatsp[0] = dot_prod(vhats,xhatl);
     hhatip[0] = dot_prod(hhati,xhatl);
     hhatsp[0] = dot_prod(hhats,xhatl);
     yhatp[0] = xhatl[1];
     zhatp[0] = xhatl[2];
     vhatip[1] = dot_prod(vhati,yhatl);
     vhatsp[1] = dot_prod(vhats,yhatl);
     hhatip[1] = dot_prod(hhati,yhatl);
     hhatsp[1] = dot_prod(hhats,yhatl);
     xhatp[1] = yhatl[1];
     yhatp[1] = yhatl[1];
     zhatp[1] = xhatl[2];
     vhatip[2] = dot_prod(vhati,zhatl);
     vhatsp[2] = dot_prod(vhats,zhatl);
     hhatip[2] = dot_prod(hhati,zhatl);
     hhatsp[2] = dot_prod(hhats,zhatl);
     vhatip[3] = 1.0f;
     vhatsp[3] = 1.0f;
     hhatip[3] = 1.0f;
     hhatsp[3] = 1.0f;
     xhatp[2] = zhatl[1];
     yhatp[2] = zhatl[1];
     zhatp[2] = zhatl[2];
     xhatp[3] = 1.0f;
     yhatp[3] = 1.0f;
     zhatp[3] = 1.0f;
     cross_prod(hhatip,vhatip,khatip);
     khsdzhl = 0.0f;
     khsdzhl = dot_prod(khats,zhatl);
     cross_prod(hhatsp,vhatsp,khatsp);
     cosbeta = 0.0f;
     cosbeta = dot_prod(khati,zhatl);
     area = 3.141592654f*radius*radius/(10000.0f);
     Pzz = area*(epsr-1.0f);
     Pxx = 2.0f*Pzz/(epsr+1.0f);
     Pyy = Pxx;
     U = 0.5f*rad_k0*(length/100.0f)*(khsdzhl-cosbeta);
     if(absf(U)>0.01f) {
        sfact = sinf(U)/u;
     } else {
        sfact = 1.0f;
     }
     pdotv[0] = Pxx*vhatip[0];
     pdotv[1] = Pyy*vhatip[1];
     pdotv[2] = Pzz*vhatip[2];
     t10 = khatsp[1]*khatsp[1];
     t11 = khatsp[2]*khatsp[2];
     t12 = khatsp[0]*khatsp[1];
     t13 = khatsp[0]*khatsp[2];
     t14 = khatsp[0]*khatsp[0];
     Sv[0] = pdotv[0]*(-(t10+t11))+
             pdotv[1]*t12+
	     pdotv[2]*t13;
     pdoth[0] = Pxx*hhatip[0];
     Sv[1] = pdotv[0]*t12+
             pdotv[1]*(-(t11+t14))+
	     pdotv[2]*(khatsp[2]*khatsp[1]);
     pdoth[1] = Pyy*hhatip[1];
     Sv[2] = pdotv[0]*t13+
             pdotv[1]*(khatsp[1]*khatsp[2])+
	     pdotv[2]*(-(t14+t10));
     pdoth[2] = Pzz*hhatip[2];
     Sh[0] = pdoth[0]*(-(t10+t11))+
             pdoth[1]*t12+
	     pdoth[2]*t13;
     Sh[1] = pdoth[0]*t12+
             pdoth[1]*(-(t11+t14))+
	     pdoth[2]*(khatsp[2]*khatsp[1]);
     Sh[2] = pdoth[0]*t13+
             pdoth[1]*(khatsp[1]*khatsp[2])+
	     pdoth[2]*(-t14+t11);
     Svv = Sv[0]*vhatsp[0]+Sv[1]*vhatsp[1]+Sv[2]*vhatsp[2];
     Shv = Sv[0]*hhatsp[0]+Sv[1]*hhatsp[1]+Sv[2]*hhatsp[2];
     Svh = Sh[0]*vhatsp[0]+Sh[1]*vhatsp[1]+Sh[2]*vhatsp[2];
     Shh = Sh[0]*hhatsp[0]+Sh[1]*hhatsp[1]+Sh[2]*hhatsp[2];
     mfact = -sfact*0.25f*rad_k0*rad_k*length/(1256.637061435917295f);
#if defined __INTEL_COMPILER
     __assume_aligned(scat_mat,32);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     scat_mat = (std::complex<float>*)__builtin_assume_aligned(scat_mat,32);
#endif
     scat_mat[0] = Svv*mfact;
     scat_mat[1] = Shv*mfact;
     scat_mat[2] = Svh*mfact;
     scat_mat[3] = Shh*mfact;
}


