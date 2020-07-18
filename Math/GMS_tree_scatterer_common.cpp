
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
		      const float lvar,
		      const float lvarv,
		      const float dens,
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
                           Leaf_PO_approximation(thinc,
			                         phinc,
						 thsc,
						 phsc,
						 thdr,
						 phdr,
						 rad_freq,
						 rad_k0,
						 rad_wv,
						 lmg,
						 lrho,
						 ldens,
						 ldiam,
						 lthick,
						 epsr,
						 &scat2x2m[0]);
						 
		      }
		       else {
                             Leaf_Rayleigh_scattering(thinc,
			                              phinc,
						      thsc,
						      phsc,
						      thdr,
						      phdr,
						      rad_freq,
						      rad_k0,
						      rad_wv,
						      lmg,
						      lrho,
						      ldens,
						      ldiam,
						      lthick,
						      epsr,
						      &scat2x2m[0]);
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
                               t1 = l4x4phm_t1[k][l][0]+orient_distr*stokes4x4m[Ix2D(k,4,l)];
			       l4x4phm_t1[k][l][0] = t1;
			   }
		       }
#else
#include "l4x4phm_t1_1_1_1.c"
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
#include "l4x4phm_t1_1_1_1.c"
#endif
                  // Case 2
		  thinc =  3.141592653589793f-theta;
                  thsc  =  3.141592653589793f-theta;
                  phinc =  0.0f
                  phsc  =  3.141592653589793f;
		  if(po) {
                           Leaf_PO_approximation(thinc,
			                         phinc,
						 thsc,
						 phsc,
						 thdr,
						 phdr,
						 rad_freq,
						 rad_k0,
						 rad_wv,
						 lmg,
						 lrho,
						 ldens,
						 ldiam,
						 lthick,
						 epsr,
						 &scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,
			                              phinc,
						      thsc,
						      phsc,
						      thdr,
						      phdr,
						      rad_freq,
						      rad_k0,
						      rad_wv,
						      lmg,
						      lrho,
						      ldens,
						      ldiam,
						      lthick,
						      epsr,
						      &scat2x2m[0]);
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
#include "l4x4phm_t1_1_1_1.c"
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
#include "l4x4phm_t1_1_1_1.c"
#endif
                    // Case 3
		    thinc = theta;
                    thsc  = theta;
                    phinc = 3.141592653589793f;
                    phsc  = 0.0f;
		    if(po) {
                           Leaf_PO_approximation(thinc,
			                         phinc,
						 thsc,
						 phsc,
						 thdr,
						 phdr,
						 rad_freq,
						 rad_k0,
						 rad_wv,
						 lmg,
						 lrho,
						 ldens,
						 ldiam,
						 lthick,
						 epsr,
						 &scat2x2m[0]);
						 
		    }
		    else {
                            Leaf_Rayleigh_scattering(thinc,
			                              phinc,
						      thsc,
						      phsc,
						      thdr,
						      phdr,
						      rad_freq,
						      rad_k0,
						      rad_wv,
						      lmg,
						      lrho,
						      ldens,
						      ldiam,
						      lthick,
						      epsr,
						      &scat2x2m[0]);
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
#include "l4x4phm_t1_1_1_2.c"
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
#include "l4x4phm_t1_1_1_2.c"
#endif
                    // Case 4:
		    
		  } // end for(ii=0 loop
	       } // end if(orient_distr block
	   } // end for(jj=0 loop
       } // end nth1!=0 block
     
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
           scat_mat[3] = const*(w1+gamhe_c1+w2*gamhe_c2);
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

     __attribute__((aligned(16))) xhat[4];
     __attribute__((aligned(16))) yhat[4];
     __attribute__((aligned(16))) zhat[4];
     __attribute__((aligned(16))) xhatl[4];
     __attribute__((aligned(16))) yhatl[4];
     __attribute__((aligned(16))) zhatl[4];
     __attribute__((aligned(16))) khati[4];
     __attribute__((aligned(16))) khats[4];
     __attribute__((aligned(16))) hhati[4];
     __attribute__((aligned(16))) vhati[4];
     __attribute__((aligned(16))) hhats[4];
     __attribute__((aligned(16))) vhats[4];
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


