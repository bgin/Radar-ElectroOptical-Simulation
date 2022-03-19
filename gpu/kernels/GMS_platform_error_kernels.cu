
#include <cuda_runtime.h>
#include <cuda.h>
#include "GMS_cpu_config.cuh"
#include "GMS_platform_error_kernels.cuh"


#define z    0.0174532925199432957692f
#define PI   3.1415926535897932384626f
#define PI2  6.2831853071795864769253f
#define tb   0.0000000001f


__forceinline__
__device__ float
v_sum_pattern(const float u) {
   const float uu = u*u;
   return (expf(-1.3866f*uu));
}


__forceinline__
__device__ float
v_diff_pattern(const float u) {
    
   return (2.345f*u*v_sum_pattern(u));
}


__global__
void platform_orient_err_kernel1(const float * __restrict__ psi,
                                 const float * __restrict__ theta,
                                 const float * __restrict__ da1, 
                                 const float * __restrict__ da2,
                                 const float * __restrict__ da3,
                                 float * __restrict__ dpsi1,//ang,min, azimuth measurement error
                                 float * __restrict__ dth1, //ang,min, elevation measurement error
                                 const uint32_t n_threads) {

    uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
    if(tid < n_threads) {
       const float xpsi  = psi[tid];
       const float xth   = theta[tid];
       const float zpsi  = z*xpsi
       const float zth   = z*xth;
       const float xda1  = da1[tid];
       const float czpsi = cosf(zpsi);
       const float szpsi = sinf(zpsi);
       const float thc   = 90.0f-xth;
       const float czth  = cosf(zth);
       const float szth  = sinf(zth);
       const float arg1  = z*xda1/60.0f;
       const float xda2  = da2[tid];
       const float c1    = cosf(arg1);
       const float s1    = sinf(arg1);
       const float arg2  = z*xda2/60.0f;
       const float xda3  = da3[tid];
       const float c2    = cosf(arg2);
       const float s2    = sinf(arg2);
       const float arg3  = z*xda3/60.0f;
       const float c3    = cosf(arg3);
       const float s3    = sinf(arg3);
       const float A0    = czpsi*szth*c1*c2;
       const float B0    = -czpsi*szth*(s1*c3+c1*s2*s3);
       const float C0    = czpsi*szth*(s1*s3-c1*s2*c3);
       const float A1    = szpsi*szth*s1*c2;
       const float B1    = -szpsi*szth*(s1*s2*s3-c1*c3);
       const float C1    = -szpsi*szth*(c1*s3+s1*s2*c3);
       const float A2    = czth*s2;
       const float B2    = czth*c2*s3;
       const float C2    = czth*c2*c3;
       const float A012  = A0+A1+A2;
       const float Asqr  = A012*A012;
       const float B012  = B0+B1+B2;
       const float Bsqr  = B012*B012;
       const float C012  = C0+C1+C2;
       const float Csqr  = C012*C012;
       const float U     = A012/(Asqr+Bsqr);
       if(zpsi>=0.0f && z<=PI) {
          dpsi1[tid] = (acosf(U)-zpsi)/z*60.0f;
       }
       else {
          dpsi1[tid] = (PI2-acosf(U)-zpsi)/z*60.0f;
       }
       const float V = C012/(Asqr+Bsqr+Csqr);
       dth1[tid]     = -((acosf(V)-z*xth)/z*60.0f);
    }
}


__global__
void platform_orient_err_kernel2(const float * __restrict__ psi,
                                 const float * __restrict__ theta,
                                 const float * __restrict__ da1, 
                                 const float * __restrict__ da2,
                                 const float * __restrict__ da3,
                                 float * __restrict__ dpsi1,//ang,min, azimuth measurement error
                                 float * __restrict__ dth1, //ang,min, elevation measurement error
                                 const uint32_t n) {

     uint32_t tid    = blockIdx.x*blockDim.x+threadIdx.x;
     uint32_t stride = blockDim.x*gridDim.x;
     for(uint32_t i = tid; i < n; i += stride) {
         const float xpsi  = psi[i];
         const float xth   = theta[i];
         const float zpsi  = z*xpsi
         const float zth   = z*xth;
         const float xda1  = da1[i];
         const float czpsi = cosf(zpsi);
         const float szpsi = sinf(zpsi);
         const float thc   = 90.0f-xth;
         const float czth  = cosf(zth);
         const float szth  = sinf(zth);
         const float arg1  = z*xda1/60.0f;
         const float xda2  = da2[i];
         const float c1    = cosf(arg1);
         const float s1    = sinf(arg1);
         const float arg2  = z*xda2/60.0f;
         const float xda3  = da3[i];
         const float c2    = cosf(arg2);
         const float s2    = sinf(arg2);
         const float arg3  = z*xda3/60.0f;
         const float c3    = cosf(arg3);
         const float s3    = sinf(arg3);
         const float A0    = czpsi*szth*c1*c2;
         const float B0    = -czpsi*szth*(s1*c3+c1*s2*s3);
         const float C0    = czpsi*szth*(s1*s3-c1*s2*c3);
         const float A1    = szpsi*szth*s1*c2;
         const float B1    = -szpsi*szth*(s1*s2*s3-c1*c3);
         const float C1    = -szpsi*szth*(c1*s3+s1*s2*c3);
         const float A2    = czth*s2;
         const float B2    = czth*c2*s3;
         const float C2    = czth*c2*c3;
         const float A012  = A0+A1+A2;
         const float Asqr  = A012*A012;
         const float B012  = B0+B1+B2;
         const float Bsqr  = B012*B012;
         const float C012  = C0+C1+C2;
         const float Csqr  = C012*C012;
         const float U     = A012/(Asqr+Bsqr);
         if(zpsi>=0.0f && z<=PI) {
            dpsi1[i] = (acosf(U)-zpsi)/z*60.0f;
         }
         else {
            dpsi1[i] = (PI2-acosf(U)-zpsi)/z*60.0f;
         }
         const float V = C012/(Asqr+Bsqr+Csqr);
         dth1[i]     = -((acosf(V)-z*xth)/z*60.0f);
     }
}


void platform_orient_err_cuda(const float * __restrict__ psi,
                              const float * __restrict__ theta,
                              const float * __restrict__ da1, 
                              const float * __restrict__ da2,
                              const float * __restrict__ da3,
                              float * __restrict__ dpsi1,//ang,min, azimuth measurement error
                              float * __restrict__ dthi1, //ang,min, elevation measurement error
                              const uint32_t n_threads,
                              const uint32_t type,
                              const uint32_t n) {

     if(type==1) {
         uint32_t threadsBlock = 32;
         uint32_t blocksGrid  = (n_threads + threadsBlock - 1) / threadsBlock;
         platform_orient_err_kernel1<<<blocksGrid,threadsBlock>>>(psi,theta,da1,da2,da3,dpsi1,dthi1,n_threads);
     }
     else if(type==2) {
          uint32_t threadsBlock = 256;
          uint32_t blocksGrid  = (n + threadsBlock - 1) / threadsBlock;
          platform_orient_err_kernel2<<<blocksGrid,threadsBlock>>>(psi,theta,da1,da2,da3,dpsi1,dthi1,n);
     }
}


__global__
void platform_pos_err_kernel1(const float * __restrict__ R,
                              const float * __restrict__ psi,
                              const float * __restrict__ theta,
                              const float * __restrict__ dx1,
                              const float * __restrict__ dx2,
                              const float * __restrict__ dx3,
                              float * __restrict__ dpsi2, //min, azimuth error
                              float * __restrict__ dth2,  //min, elevation error
                              float * __restrict__ dR,    //m,   range error
                              const uint32_t n_threads) {

     uint32_t tid = blockDim.x*blockIdx.x+threadIdx.x;
     if(tid < n_threads) {
        const float xR    = R[tid];
	const float xpsi  = psi[tid];
	const float zpsi  = z*xpsi;
	const float xth   = theta[tid];
        const float zth   = z*xth;
	const float szth  = cephes_sinf(zth);
        const float czth  = cephes_cosf(zth);
        const float xdx1  = dx1[tid];
	const float xdx12 = xdx1*xdx1;
	const float xdx2  = dx2[tid];
	const float xdx22 = xdx2*xdx2;
	const float xdx3  = dx3[tid];
	const float xdx32 = xdx3*xdx3;
	const float eta1  = sqrtf(xdx12+xdx22+xdx32);
	const float eta1R = eta1/xR;
	const float eta1R2= eta1R*eta1R;
	const float eta2  = atanf(xdx2/(xdx1+tb));
	const float num   = sqrtf(xdx12+xdx22);
	const float eta3  = atanf(num/xdx3+tb);
	const float seta3 = sinf(eta3);
	const float ceta3 = cosf(eta3);
	const float lterm = cosf(zpsi)*cephes_sinf(zpsi)
	const float rterm = eta1R*cosf(eta2)*seta3;
	const float D0    = lterm-rterm;
        const float lterm2= sinf(zpsi)*szth;
	const float rterm2= eta1R*sinf(eta2)*seta3;
	const float D1    = lterm2-rterm2;
	const float D2    = czth-eta1R*ceta3;
	const float Dk01  = D0*D0+D1*D1;
	const float U     = D0/sqrtf(Dk01)
        if(zpsi>=0.0f && zpsi<=PI) {
              dpsi2[tid] = (acosf(U)-zpsi)/z*60.0f;
	}
	else {
              dpsi2[tid] = (acosf(U)-zpsi)/z*60.0f;
	}
	const float Dk012 = DK01+D2*D2;
	const float V     = D2/sqrtf(Dk012);
        dth2[tid] = -((acosf(V)-zth)/z*60.0f);
	const float lterm3= szth*cosf(zpsi-eta2)*seta3;
	const float rterm3= czth*ceta3;
	const float aR    = lterm3+rterm3;
	const float lterm4= 1.0f-2.0f*eta1R*aR+eta1R2;
	const float delR  = sqrtf(lterm4);
	dR[tid]           = xR*(delR-1.0f);
     }
}
                                 
