
#include <cuda_runtime.h>
#include <cuda.h>
#include "GMS_cpu_config.cuh"
#include "GMS_rcs.cuh"



static const float sig0 =   1.0f;
static const float PI   =   3.14159265358979323846264338328f;
static const float PI2  =   6.283185307179586476925286766559f;
static const float PI4  =   12.566370614359172953850573533118f;
static const float z    =   0.017453292519943295769236907685f; //coeff deg-to-rad conversion


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
    for(int32_t i = tid; i < n; i += stride) {
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
     for(int32_t i = tid; i < n; i += stride) {  
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
             uint32_t threadsBlock = 32;
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

