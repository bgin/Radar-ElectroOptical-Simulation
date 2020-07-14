
#include <limits>
#include "GMS_tree_scatterer_common.h"


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
     
}
