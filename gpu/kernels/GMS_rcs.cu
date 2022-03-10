
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
void empirical_K_kernel(const float * __restrict rcs, //// m^2, RCS monostatic target
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


void empirical_K_cuda(const float * __restrict rcs,
                      const float * __restrict A,
                      const float gam,
                      float * __restrict K,
                      const uint32_t n_threads) {

        uint32_t threadsBlock = 32;
        uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
        emprirical_K_kernel<<<blocksGrid,threadsBlock>>>(rcs,A,gam,K,n_threads);
}
