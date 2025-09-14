
#include <cuda_runtime.h>
#include <cuda.h>
#include "GMS_cpu_config.cuh"
#include "GMS_rcs_kernels.cuh"



#define sig0    1.0f
#define PI      3.14159265358979323846264338328f
#define PI2     6.283185307179586476925286766559f
#define PI4     12.566370614359172953850573533118f
#define z       0.017453292519943295769236907685f 
#define PI16    50.2654824574366918154023f 

__global__
void empirical_K_kernel1(const float * __restrict rcs, //// m^2, RCS monostatic target
                        const float * __restrict A,   //   m^2, Area of target projected to the normal radar beam 
                        const float              gam, //   m,   wavelength
                        float *       __restrict K,   //   Empirical K factor
                        const uint32_t n_threads   ) {

     uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
     if(tid < n_threads) {

        const float PI24 = 0.74159265358979323846264338328f;
        const float tA   = A[tid];
        const float A2   = tA*tA;
        const float trcs = rcs[tid];
        const float gam2 = gam*gam*trcs;
        const float t0   = PI4*(A2/gam2);
        K[tid]           = logf(t0)/PI24;
     }
}


__global__
void empirical_K_kernel2const float * __restrict rcs, //// m^2, RCS monostatic target
                        const float * __restrict A,   //   m^2, Area of target projected to the normal radar beam 
                        const float              gam, //   m,   wavelength
                        float *       __restrict K,   //   Empirical K factor
                        const uint32_t n   ) {

    uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
    uint32_t stride = blockDim.x*gridDim.x;
    for(uint32_t i = tid; i < n; i += stride) {
         const float PI24 = 0.74159265358979323846264338328f;
         const float tA   = A[i];
         const float A2   = tA*tA;
         const float trcs = rcs[i];
         const float gam2 = gam*gam*trcs;
         const float t0   = PI4*(A2/gam2);
         K[i]           = logf(t0)/PI24; 
    }
}


void empirical_K_cuda(const float * __restrict rcs,
                      const float * __restrict A,
                      const float gam,
                      float * __restrict K,
                      const uint32_t n_threads,
                      const uint32_t n,
                      const uint32_t type) {
        
        if(type==1) {
           uint32_t threadsBlock = 32;
           uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
           emprirical_K_kernel1<<<blocksGrid,threadsBlock>>>(rcs,A,gam,K,n_threads);
        }else if(type==2) {
           uint32_t threadsBlock = 256;
           uint32_t blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
           emprirical_K_kernel1<<<blocksGrid,threadsBlock>>>(rcs,A,gam,K,n);
        }
}




__forceinline__ __device__ 
float wavenumber_r4(const float gamma) { //m, wavelength
       return (PI2/gamma);
}


__global__
void effective_rcs_kernel1(const float gamma, //m, wavelength
                           const float * __restrict__ R, //nm, target range
                           const float h_a,              //ft, antenna height
                           const float h_t,              //ft, target  height
                           float * __restrict__ rcs_eff, //m,  effective rcs
                           float * __restrict__ rcs_eff_db, //m, effective rcs,
                           const uint32_t n_threads) {

       uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
       if(tid < n_threads) {
          const float t0 = PI2*ha_*h_t;
          const float t1 = R[tid]*gamma;
          const float F  = sinf(t0/t1);
          const float F4 = F*F*F*F;
          rcs_eff[tid]   = sig0*F4;
          rcs_eff_db[tid]= 10.0f*logf(rcs_eff[tid]);
       }
}


__global__
void effective_rcs_kernel2(const float gamma, //m, wavelength
                           const float * __restrict__ R, //nm, target range
                           const float h_a,              //ft, antenna height
                           const float h_t,              //ft, target  height
                           float * __restrict__ rcs_eff, //m,  effective rcs
                           float * __restrict__ rcs_eff_db, //m, effective rcs,
                           const uint32_t n) {
    
     uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
     uint32_t stride = blockDim.x*gridDim.x;
     for(uint32_t i = tid; i < n; i += stride) {  
          const float t0 = PI2*ha_*h_t;
          const float t1 = R[i]*gamma;
          const float F  = sinf(t0/t1);
          const float F4 = F*F*F*F;
          rcs_eff[i]   = sig0*F4;
          rcs_eff_db[tid]= 10.0f*logf(rcs_eff[i]);
     }
}


void effective_rcs_cuda(const float gamma
                        const float * __restrict__ R, 
                        const float h_a, 
                        const float h_t, 
                        float * __restrict__ rcs_eff,
                        float * __restrict__ rcs_eff_db,
                        const uint32_t n_threads,
                        const uint32_t n,
                        const uint32_t type) {

       if(type==1) {
            uint32_t threadsBlock = 32;
            uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
            effective_rcs_kernel1<<<blocksGrid,threadsBlock>>>(gamma,R,h_a,h_t,rcs_eff,rcs_eff_db,n);
       }else if(type==2) {
             uint32_t threadsBlock = 256;
             uint32_t blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
             effective_rcs_kernel2<<<blocksGrid,threadsBlock>>>(gamma,R,h_a,h_t,rcs_eff,rcs_eff_db,n);
      }
}


__global__
void bistatic_target_rcs_kernel1(const float sig0, // m^2, RCS monostatic target
		                 const float K,    // empirical constant defined by target configuration and complexity  // use function emprirical_K for computation
	                         const float * __restrict__ Beta, // deg, bistatic angle (for input=1)
		                 const float * __restrict__ R1,   // m, transmitter - target range (input=2)
			         const float * __restrict__ R2,   // m, receiver    - target range (input=2)
			         const float B,    // m, baseline
			         const int32_t type, // input switch (1,2)
			         float * __restrict sigma, // RCS                  
                                 float * __restrict sigma_db,
                                 const uint32_t n_threads) {

       uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
       if(type==1) {
          if(tid < n_threads)
             const float alpha = Beta[tid]*z;
             const float t0    = fabsf(alpha)-2.4f*K-1.0f;
	     const float t1    = 1.0f+expf(t0);
	     sigma[tid]        = sig0*t1;
	     sigma_db[tid]     = 10.0f*log10f(sigma[tid]+0.00000000001f);
            }
	}
	else if(type==2) {
           if(tid < n_threads) {
             const float tR1   = R1[tid];
             const float tR2   = R2[tid];
             const float t0    = 1.0f/(2.0f*tR1*tR2);
             const float tB    = B[tid];
	     const float R12   = tR1*tR1;
	     const float R22   = tR2*tR2;
	     const float B2    = tB*tB;
	     const float t1    = R12+R22+B2;
	     const float gam   = acosf(t0*t1);
	     const float alpha = gam;
	     const float t2    = K*fabsf(alpha)-2.4f*K-1.0f;
	     const float t3    = 1.0f+expf(t0);
	     sigma[tid]        = sig0*t1;
	     sigma_db[tid]     = 10.0f*log10f(sigma[tid]+0.00000000001f);


          }
       }
}   



__global__
void bistatic_target_rcs_kernel2(const float sig0, // m^2, RCS monostatic target
		                 const float K,    // empirical constant defined by target configuration and complexity  // use function emprirical_K for computation
	                         const float * __restrict__ Beta, // deg, bistatic angle (for input=1)
		                 const float * __restrict__ R1,   // m, transmitter - target range (input=2)
			         const float * __restrict__ R2,   // m, receiver    - target range (input=2)
			         const float B,    // m, baseline
			         const int32_t type, // input switch (1,2)
			         float * __restrict sigma, // RCS                  
                                 float * __restrict sigma_db,
                                 const uint32_t n) {

       uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
       uint32_t stride = blockDim.x*gridDim.x;
       if(type==1) {
          for(int32_t i = tid; i < n; i += stride) {
              const float alpha = Beta[i]*z;
              const float t0    = fabsf(alpha)-2.4f*K-1.0f;
	      const float t1    = 1.0f+expf(t0);
	      sigma[i]        = sig0*t1;
	      sigma_db[i]     = 10.0f*log10f(sigma[i]+0.00000000001f);
            }
	}
	else if(type==2) {
          for(uint32_t i = tid; i < n; i += stride) {
              const float tR1   = R1[i];
              const float tR2   = R2[i];
              const float t0    = 1.0f/(2.0f*tR1*tR2);
              const float tB    = B[i];
	      const float R12   = tR1*tR1;
	      const float R22   = tR2*tR2;
	      const float B2    = tB*tB;
	      const float t1    = R12+R22+B2;
	      const float gam   = acosf(t0*t1);
	      const float alpha = gam;
	      const float t2    = K*fabsf(alpha)-2.4f*K-1.0f;
	      const float t3    = 1.0f+expf(t0);
	      sigma[i]          = sig0*t1;
	      sigma_db[i]       = 10.0f*log10f(sigma[i]+0.00000000001f);


          }
       }
}  


void bistatic_target_rcs_cuda(const float sig0, 
		              const float K, 
	                      const float * __restrict__ Beta,
		              const float * __restrict__ R1,
			      const float * __restrict__ R2,
			      const float B,   
			      const int32_t type, 
			      float * __restrict sigma,               
                              float * __restrict sigma_db,
                              const uint32_t n_threads,
                              const uint32_t kernel_type,
                              const uint32_t n) {

       if(kernel_type==1) {
            uint32_t threadsBlock = 32;
            uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
            bistatic_target_rcs_kernel1<<<blocksGrid,threadsBlock>>>(sig0,K,Beta,R1,R2,B,type,sigma,sigma_db,n_threads);
       }else if(kernel_type==2) {
             uint32_t threadsBlock = 256;
             uint32_t blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
             bistatic_target_rcs_kernel2<<<blocksGrid,threadsBlock>>>(sig0,K,Beta,R1,R2,B,type,sigma,sigma_db,n);
      }
}


__global__
void antenna_rcs_kernel1(const float * __restrict__ Ae, //m^2, antenna effective apperture
                         const float gamma,             //m,   wavelength
                         const float G,                 // reflection coeff
                         float * __restrict__ sigma,    // m^2, rcs
                         float * __restrict__ sigma_db,   // dBsm, rcs with respect to area of 1 m^2
                         const uint32_t n_threads) {

        uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
        if(tid < n_threads) {
           const float tAe = Ae[tid];
           const float Ae2 = tAe*tAe;
           const float gamma2 = gamma*gamma;
           const float t0  = Ae2/gamma2;
           sigma[tid]    = PI4*t0*G;
           sigma_db[tid] = 10.0f*log10f(sigma[tid]);
        }
}


__global__
void antenna_rcs_kernel2(const float * __restrict__ Ae, //m^2, antenna effective apperture
                         const float gamma,             //m,   wavelength
                         const float G,                 // reflection coeff
                         float * __restrict__ sigma,    // m^2, rcs
                         float * __restrict__ sigma_db,   // dBsm, rcs with respect to area of 1 m^2
                         const uint32_t n) {

       uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
       uint32_t stride = blockDim.x*gridDim.x;
       for(uint32_t i = tid; i < n; i += stride) {
           const float tAe = Ae[i];
           const float Ae2 = tAe*tAe;
           const float gamma2 = gamma*gamma;
           const float t0  = Ae2/gamma2;
           sigma[i]    = PI4*t0*G;
           sigma_db[i] = 10.0f*log10f(sigma[i]);
      }  
}


void antenna_rcs_cuda(   const float * __restrict__ Ae,
                         const float gamma,            
                         const float G,                
                         float * __restrict sigma,   
                         float * __restrict__ sigma_db,   
                         const uint32_t n_threads,
                         const uint32_t type,
                         const uint32_t n) {

     if(type==1) {
         uint32_t threadsBlock = 32;
         uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
         antenna_rcs_kernel1<<<blocksGrid,threadsBlock>>>(Ae,gamma,G,sigma,sigma_db,n_threads);
     }else if(type==2) {
          uint32_t threadsBlock = 256;
          uint32_t blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
          antenna_rcs_kernel2<<<blocksGrid,threadsBlock>>>(Ae,gamma,G,sigma,sigma_db,n);
    }
}


__global__
void bird_insect_rcs_kernel(const float * __restrict__ W,
                            float * __restrict__ rcs
                            const uint32_t n_threads) {// gram, bird/insect weight
     
      uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
      if(tid < n_threads) {
         const float tW = W[tid];
         rcs[tid] = -46.0f+5.8f*log10f(tW);
      } 
}


void bird_insect_rcs_cuda(const float * __restrict__ W,
                          float * __restrict__ rcs,
                          const uint32_t n_threads) {

     uint32_t threadsBlock = 32;
     uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
     bird_insect_rcs_kernel<<<blocksGrid,threadsBlock>>>(W,rcs,n_threads);
}


__global__
void cone_ogive_rcs_kernel1(const float gamma, //m, wavelength
                            const float * __restrict__ cha, // deg, cone half-angle, incidence nose-on
                            float * __restrict__ sig1,
                            const uint32_t n_threads) {
 // This model applies to infintely large cone, or ogive larger than incident wavelength
      uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
      if(tid < n_threads) {
         const float tcha = cha[tid];
         const float gamma2 = gamma*gamma;
         const float ratio  = gamm2/PI16;
         const float targ   = z*tcha
         const float targ4  = targ*targ*targ*targ;
         const float tan4   = tanf(targ4);
         sig1[tid]          = ratio*tan4;
      }
}


__global__
void cone_ogive_rcs_kernel2(const float gamma, //m, wavelength
                            const float * __restrict__ cha, // deg, cone half-angle, incidence nose-on
                            float * __restrict__ sig1,
                            const uint32_t n) {

      uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
      uint32_t stride = blockDim.x*gridDim.x;
      for(uint32_t i = tid; i < n; i += stride) {
          const float tcha = cha[i];
          const float gamma2 = gamma*gamma;
          const float ratio  = gamm2/PI16;
          const float targ   = z*tcha
          const float targ4  = targ*targ*targ*targ;
          const float tan4   = tanf(targ4);
          sig1[i]          = ratio*tan4;
      }
}


void cone_ogive_rcs_cuda(const float gamma, 
                         const float * __restrict__ cha,
                         float * __restrict__ sig1,
                         const uint32_t n_threads,
                         const uint32_t type,
                         const uint32_t n) {

     if(type==1) {
        uint32_t threadsBlock = 32;
        uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
        cone_ogive_rcs_kernel1<<<blocksGrid,threadsBlock>>>(gamma,cha,sig1,n_threads);
     }
     else if(type==2) {
        uint32_t threadsBlock = 256;
        uint32_t blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
        cone_ogive_rcs_kernel2<<<blocksGrid,threadsBlock>>>(gamma,cha,sig1,n);
     }
}


__global__
void cylinder_rcs_kernel1(const float * __restrict__ Rcyl, //m, radius of cylinder
                          const float * __restrict__ Lcyl, //m, length of cylinder
                          const float gamma,               //m, wavelength
                          const float theta,               //deg, angle between cylinder axis and radar line-of-sight
                          float * __restrict__ sig,
                          const uint32_t n_threads) {

       uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
       if(tid < n_threads) {
          const float tRcyl = Rcyl[tid];
          const float tLcyl = Lcyl[tid];
          const float t0    = PI4*tRcyl;
          const float zth   = z*theta;
          const float sinz  = sinf(zth);
          const float term  = (t0*sinz)/gamma;
          if(term<1.0f) sig[tid] = 0.0f;
          const float cosz  = cosf(zth);
          const float term1 = (tRcyl*gamma*sinz)/PI2;
          const float sarg1 = (PI2*tLcyl)/gamma;
          const float sarg2 = sarg1*cosz;
	  const float term2 = sinf(sarg2/cosz);
	  const float term22= term2*term2;
          sig[tid]          = term1*term22;
       } 
}


__global__
void cylinder_rcs_kernel2(const float * __restrict__ Rcyl, //m, radius of cylinder
                          const float * __restrict__ Lcyl, //m, length of cylinder
                          const float gamma,               //m, wavelength
                          const float theta,               //deg, angle between cylinder axis and radar line-of-sight
                          float * __restrict__ sig,
                          const uint32_t n) {

      uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
      uint32_t stride = blockDim.x*gridDim.x;
      for(uint32_t i = tid; i < n; i += stride) {
          const float tRcyl = Rcyl[i];
          const float tLcyl = Lcyl[i];
          const float t0    = PI4*tRcyl;
          const float zth   = z*theta;
          const float sinz  = sinf(zth);
          const float term  = (t0*sinz)/gamma;
          if(term<1.0f) sig[i] = 0.0f;
          const float cosz  = cosf(zth);
          const float term1 = (tRcyl*gamma*sinz)/PI2;
          const float sarg1 = (PI2*tLcyl)/gamma;
          const float sarg2 = sarg1*cosz;
	  const float term2 = sinf(sarg2/cosz);
	  const float term22= term2*term2;
          sig[i]          = term1*term22;
       } 
}


void cylinder_rcs_cuda(const float * __restrict__ Rcyl,
                       const float * __restrict__ Lcyl,
                       const float gamma,              
                       const float theta,            
                       float * __restrict__ sig,
                       const uint32_t n_threads,
                       const uint32_t type,
                       const uint32_t n) {

     if(type==1) {
        uint32_t threadsBlock = 32;
        uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
        cylinder_rcs_kernel1<<<blocksGrid,threadsBlock>>>(Rcyl,Lcyl,gamma,theta,sig,n_threads);
     }
     else if(type==2) {
        uint32_t threadsBlock = 256;
        uint32_t blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
        cylinder_rcs_kernel1<<<blocksGrid,threadsBlock>>>(Rcyl,Lcyl,gamma,theta,sig,n);
     }
}


__global__
void disk_rcs_kernel1(const float * __restrict__ Rd, //m^2, radius if disk
                      const float * __restrict__ theta, //deg, angle relative to disk normal
                      const float gamma,             //m, wavelength
                      const int32_t type, //// disk type: 1=electrically large disk, 2=electrically small circular disk
                      float * __restrict__ sig,
                      const uint32_t n_threads) {

        uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
        
           if(type==1) {
             if(tid < n_threads) {
                 const float tRd = Rd[tid];
                 const float tth = theta[tid];
                 const float k = wavenumber_r4(gamma);
                 const float gamma2 = gamma*gamma;
                 const float carg = z*theta;
                 const float coszt= cosf(carg*carg);
                 const float t0 = PI*tRd*tRd;
                 const float tleft = t1/gamma;
                 const float sinzt = sinf(carg);
                 const float jarg  = 2.0f*k*tRd*sinzt;
                 const float bj1   = j1f(jarg);
                 const float term  = 2.0f*bj1/jarg;
                 const float tright= term*term;
                 sig[tid] = tleft*tright*coszt;
               }
            }
            else if(type==2) {
               if(tid < n_threads) {
                  const float tRd   = Rd[tid];
                  const float Rd2   = tRd*tRd;
                  const float tleft = 3534.2917352885173932704738f*Rd2;
                  const float ratio = tRd/gamma;
                  const float ratio4= ratio*ratio*ratio*ratio;
                  sig[tid] = tleft*ratio4;
            }
        }    
}


__global__
void disk_rcs_kernel2(const float * __restrict__ Rd, //m^2, radius if disk
                      const float * __restrict__ theta, //deg, angle relative to disk normal
                      const float gamma,             //m, wavelength
                      const int32_t type, //// disk type: 1=electrically large disk, 2=electrically small circular disk
                      float * __restrict__ sig,
                      const uint32_t n) {

           uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
           uint32_t stride = blockDim.x*gridDim.x;
           if(type==1) {
              for(uint32_t i = tid; i < n; i += stride) {
                  const float tRd = Rd[i];
                  const float tth = theta[i];
                  const float k = wavenumber_r4(gamma);
                  const float gamma2 = gamma*gamma;
                  const float carg = z*theta;
                  const float coszt= cosf(carg*carg);
                  const float t0 = PI*tRd*tRd;
                  const float tleft = t1/gamma;
                  const float sinzt = sinf(carg);
                  const float jarg  = 2.0f*k*tRd*sinzt;
                  const float bj1   = j1f(jarg);
                  const float term  = 2.0f*bj1/jarg;
                  const float tright= term*term;
                  sig[i] = tleft*tright*coszt;
               }
            }
            else if(type==2) {
               for(uint32_t i = tid; i < n; i += stride) {
               const float tRd   = Rd[i];
               const float Rd2   = tRd*tRd;
               const float tleft = 3534.2917352885173932704738f*Rd2;
               const float ratio = tRd/gamma;
               const float ratio4= ratio*ratio*ratio*ratio;
               sig[i] = tleft*ratio4;
            }
        }    
}


void disk_rcs_cuda(const float * __restrict__ Rd, 
                   const float * __restrict__ theta, 
                   const float gamma,            
                   const int32_t type, 
                   float * __restrict__ sig,
                   const uint32_t n_threads,
                   const uint32_t kernel_type,
                   const uint32_t n) {

     if(kernel_type==1) {
        uint32_t threadsBlock = 32;
        uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
        disk_rcs_kernel1<<<blocksGrid,threadsBlock>>>(Rd,theta,gamma,type,sig,n_threads);
     }
     else if(kernel_type==2) {
        uint32_t threadsBlock = 256;
        uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
        disk_rcs_kernel2<<<blocksGrid,threadsBlock>>>(Rd,theta,gamma,type,sig,n);
     }
}


__global__
void curved_edge_rcs_kernel(const float * __restrict__ Redg, //m, radius of edge contour, incidence edge-perpendicular
		            const float gamma,//m, wavelength
                            float * __restrict__ rcs ,
                            const uint32_t n_threads) { 

       uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
       if(tid < n_threads) {
          const float tRedg = Redg[tid];
          rcs[tid] = 0.5f*tRedg*gamma;
       }
}


void curved_edge_rcs_cuda(const float * __restrict__ Redg,
		          const float gamma,
                          float * __restrict__ rcs,
                          const uint32_t n_threads) {

      uint32_t threadsBlock = 32;
      uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
      curved_edge_rcs_kernel<<<blocksGrid,threadsBlock>>>(Redg,gamma,rcs,n_threads);
}


__global__
void straight_edge_rcs(const float * __restrict__ Ledg,
                       float * __restrict__ rcs,
                       const uint32_t n_threads) {

    uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
    if(tid < n_threads) {
       const float tLedg = Ledg[tid];
       rcs[tid] = (tLedg*tLedg)/PI;
    }
}


void straight_edge_rcs_cuda(const float * __restrict__ Ledg,
                            float * __restrict__ rcs,
                            const uint32_t n_threads) {

     uint32_t threadsBlock = 32;
     uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
     straight_edge_rcs_kernel<<<blocksGrid,threadsBlock>>>(Ledg,rcs,n_threads);
}


__global__
void ellipsoid_rcs_kernel(const float * __restrict__ a, //m, semimajor axis
                          const float * __restrict__ b, //m, semiminor axis
                          float * __restrict__ rcs,
                          const uint32_t n_threads) {

      uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
      if(tid < n_threads) {
         const float ta = a[tid];
         const float tb = b[tid];
         const float b4 = tb*tb*tb*tb;
         const float a2 = ta*ta;
         rcs[tid]       = (PI*b4)/a2;
      }   
}


void ellipsoid_rcs_cuda(  const float * __restrict__ a,
                          const float * __restrict__ b,
                          float * __restrict__ rcs,
                          const uint32_t n_threads) {

     uint32_t threadsBlock = 32;
     uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
     ellipsoid_rcs_kernel<<<blocksGrid,threadsBlock>>>(a,b,rcs,n_threads);
}


__global__
void plate_rcs_kernel1(const float * __restrict__ x, //m, plate length
                       const float * __restrict__ y, //m, plate width
                       const float theta,           //deg, angle relative to plate norm
                       const float phi,             //deg, angle between the plate normal and radar los
	               const float gamma,           //m, wavelength
                       float * __restrict__ sig,
                       const uint32_t n_threads) {         

     uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
     if(tid < n_threads) {
        const float tx     = x[tid];
        const float zth    = z*theta;
        const float ty     = y[tid];
        const float gamma2 = gamma*gamma;
        const float sinzt  = sinf(zth);
	const float xy2    = PI4*(x*y*x*y);
	const float coszt  = cosf(z*phi);
	const float zth2  = zth*zth;
	const float u     = wavenumber_r4(gamm)*x*sinzt*coszt;
	const float sinu  = sinf(u)/u;
	const float v     = wavenumber_r4(gamm)*y*sinzt*coszt;
	const float sinv  = sinf(v)/v;
	const float tleft = xy2/gamma2;
	const float tmid  = sinu*sinv*sinu*sinv;
	sig[tid]          = tleft*tmid*coszt;
     }
}


__global__
void plate_rcs_kernel2(const float * __restrict__ x, //m, plate length
                       const float * __restrict__ y, //m, plate width
                       const float theta,           //deg, angle relative to plate norm
                       const float phi,             //deg, angle between the plate normal and radar los
	               const float gamma,           //m, wavelength
                       float * __restrict__ sig,
                       const uint32_t n) {         

        uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
        uint32_t stride = blockDim.x*gridDim.x;
        for(uint32_t i = tid; i < n; i += stride) {
            const float tx     = x[i];
            const float zth    = z*theta;
            const float ty     = y[i];
            const float gamma2 = gamma*gamma;
            const float sinzt  = sinf(zth);
	    const float xy2    = PI4*(x*y*x*y);
	    const float coszt  = cosf(z*phi);
	    const float zth2  = zth*zth;
	    const float u     = wavenumber_r4(gamm)*x*sinzt*coszt;
	    const float sinu  = sinf(u)/u;
	    const float v     = wavenumber_r4(gamm)*y*sinzt*coszt;
	    const float sinv  = sinf(v)/v;
	    const float tleft = xy2/gamma2;
	    const float tmid  = sinu*sinv*sinu*sinv;
	    sig[i]            = tleft*tmid*coszt;
     }
}


void plate_rcs_cuda(   const float * __restrict__ x, 
                       const float * __restrict__ y, 
                       const float theta,           
                       const float phi,             
	               const float gamma,          
                       float * __restrict__ sig,
                       const uint32_t n_threads,
                       const uint32_t type,
                       const uint32_t n) {

    if(type==1) {
       uint32_t threadsBlock = 32;
       uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
       plate_rcs_kernel1<<<blocksGrid,threadsBlock>>>(x,y,theta,phi,gamma,sig,n_threads);
    }
    else if(type==2) {
       uint32_t threadsBlock = 256;
       uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
       plate_rcs_kernel2<<<blocksGrid,threadsBlock>>>(x,y,theta,phi,gamma,sig,n);
    }
}

