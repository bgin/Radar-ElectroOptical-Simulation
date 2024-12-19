#include "fdtd.hpp"
#include <algorithm>
#include <iostream>
#include <chrono>
#include <ctime>
#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <complex>
#include <fstream>
#include <iomanip>

#include <thread>
#include <vector>
#include <algorithm>
#include <string>

#include <atomic>
#include <mutex>
#include <condition_variable>
#include <numeric>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}

static std::mutex mutex_H, mutex_E;//, mutex_H2, mutex_E2;
static std::condition_variable cv_H, cv_E, cv_H2, cv_E2;
std::atomic<size_t> counter_H(0), counter_H2(0), counter_E(0), counter_E2(0);
int blockdim = 128;

void Barrier(std::atomic<size_t> &counter, std::mutex &mutex, std::condition_variable &cv, size_t nthreads)
{
    std::unique_lock<std::mutex> lock(mutex);
    if (++counter == nthreads){
        counter.store(0);
        cv.notify_all();
    }
    else
        cv.wait(lock, [&] { return counter == 0; });
}

template <typename T>
__global__ void initMat(T *matrix, size_t len, T val)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for(int i=idx;i<len;i+=gridDim.x*blockDim.x)
    {
        matrix[i]=val;
    }
}

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

__global__ void kernel_update_H_bulk(kernelpar *kpar)
{
//     1D grid of 1D blocks
   size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
//     1D grid of 2D blocks
//     size_t ind_global = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
//     2D grid of 1D blocks
//     size_t ind_global = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;

    // int i,j,k;
    // k = blockIdx.x * blockDim.x + threadIdx.x + kpar->k0;
    // j = blockIdx.y * blockDim.y + threadIdx.y + kpar->j0;
    // i = blockIdx.z + kpar->i0;
    // single GPU version
    // size_t ind_global = (i-kpar->i0) * kpar->Nx * kpar->Ny + (j-kpar->j0) * kpar->Nx + (k-kpar->k0);
    // multiple GPU version: decompose along X
    
    // index within each decomposed domain
    // size_t ind_global = (i-kpar->i0) * kpar->K * kpar->J + (j-kpar->j0) * kpar->K + (k-kpar->k0);
    if (ind_global >= kpar->size) return;

    int i, j, k; 
    i = ind_global/(kpar->J*kpar->K);
    j = (ind_global%(kpar->J*kpar->K))/kpar->K;
    k = (ind_global%(kpar->J*kpar->K))%kpar->K;

    if(kpar->K<kpar->Nx && kpar->k0!=0 && k==0) return;
    if(kpar->K<kpar->Nx && kpar->k0+kpar->K!=kpar->Nx && k==kpar->K-1) return;
    if(kpar->J<kpar->Ny && kpar->j0!=0 && j==0) return;
    if(kpar->J<kpar->Ny && kpar->j0+kpar->J!=kpar->Ny && j==kpar->J-1) return;
    if(kpar->I<kpar->Nz && kpar->i0!=0 && i==0) return;
    if(kpar->I<kpar->Nz && kpar->i0+kpar->I!=kpar->Nz && i==kpar->I-1) return;

    int ind_ijk, ind_ijp1k, ind_ip1jk, ind_ijkp1;
    double dExdy, dExdz, dEydx, dEydz, dEzdx, dEzdy;
    int ind_pml, ind_pml_param;
    double b, C, kappa;
    int srcind_src, srcind_size, srcind_size_offset;
    double src_t;

    ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2) + (j+1)*(kpar->K+2) + k + 1;
    ind_ijp1k = ind_ijk + kpar->K + 2;
    ind_ip1jk = ind_ijk + (kpar->J+2)*(kpar->K+2);
    ind_ijkp1 = ind_ijk + 1;

    // set up fields on the boundary
    if (kpar->bc[0] != 'P' && k == kpar->Nx - kpar->k0){
        kpar->Ey[ind_ijk] = 0.0;
        kpar->Ez[ind_ijk] = 0.0;
    }
    if (kpar->bc[1] != 'P' && j == kpar->Ny - kpar->j0){
        kpar->Ex[ind_ijk] = 0.0;
        kpar->Ez[ind_ijk] = 0.0;
    }
    if (kpar->bc[2] != 'P' && i == kpar->Nz - kpar->i0){
        kpar->Ex[ind_ijk] = 0.0;
        kpar->Ey[ind_ijk] = 0.0;
    }

    dEzdy = kpar->odx * (kpar->Ez[ind_ijp1k] - kpar->Ez[ind_ijk]);
    dEydz = kpar->odx * (kpar->Ey[ind_ip1jk] - kpar->Ey[ind_ijk]);
    dExdz = kpar->odx * (kpar->Ex[ind_ip1jk] - kpar->Ex[ind_ijk]);
    dEzdx = kpar->odx * (kpar->Ez[ind_ijkp1] - kpar->Ez[ind_ijk]);
    dEydx = kpar->odx * (kpar->Ey[ind_ijkp1] - kpar->Ey[ind_ijk]);
    dExdy = kpar->odx * (kpar->Ex[ind_ijp1k] - kpar->Ex[ind_ijk]);
    kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + kpar->dt * (dEydz - dEzdy);
    kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + kpar->dt * (dEzdx - dExdz);
    kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + kpar->dt * (dExdy - dEydx);

    // Update PML
    if (k + kpar->k0 < kpar->pml_xmin){
        ind_pml = i * kpar->J * (kpar->pml_xmin - kpar->k0) + j * (kpar->pml_xmin - kpar->k0) + k;
        ind_pml_param = kpar->pml_xmin - k - kpar->k0 - 1;
        kappa = kpar->kappa_H_x[ind_pml_param];
        b = kpar->bHx[ind_pml_param];
        C = kpar->cHx[ind_pml_param];
        kpar->pml_Eyx0[ind_pml] = C * dEydx + b * kpar->pml_Eyx0[ind_pml];
        kpar->pml_Ezx0[ind_pml] = C * dEzdx + b * kpar->pml_Ezx0[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] - kpar->dt * (kpar->pml_Eyx0[ind_pml] - dEydx + dEydx / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + kpar->dt * (kpar->pml_Ezx0[ind_pml] - dEzdx + dEzdx / kappa);
    }
    else if(k + kpar->k0 >= kpar->pml_xmax){
        ind_pml = i * kpar->J * (kpar->k0 + kpar->K - kpar->pml_xmax) + j * (kpar->k0 + kpar->K - kpar->pml_xmax)
                + k + kpar->k0 - kpar->pml_xmax;
        ind_pml_param = k + kpar->k0 - kpar->pml_xmax + kpar->w_pml_x0;
        kappa = kpar->kappa_H_x[ind_pml_param];
        b = kpar->bHx[ind_pml_param];
        C = kpar->cHx[ind_pml_param];
        kpar->pml_Eyx1[ind_pml] = C * dEydx + b * kpar->pml_Eyx1[ind_pml];
        kpar->pml_Ezx1[ind_pml] = C * dEzdx + b * kpar->pml_Ezx1[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] - kpar->dt * (kpar->pml_Eyx1[ind_pml] - dEydx + dEydx / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + kpar->dt * (kpar->pml_Ezx1[ind_pml] - dEzdx + dEzdx / kappa);
    }
    if (j + kpar->j0 < kpar->pml_ymin){
        ind_pml = i * (kpar->pml_ymin - kpar->j0) * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_ymin - j - kpar->j0 - 1;
        kappa = kpar->kappa_H_y[ind_pml_param];
        b = kpar->bHy[ind_pml_param];
        C = kpar->cHy[ind_pml_param];
        kpar->pml_Exy0[ind_pml] = C * dExdy + b * kpar->pml_Exy0[ind_pml];
        kpar->pml_Ezy0[ind_pml] = C * dEzdy + b * kpar->pml_Ezy0[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + kpar->dt * (kpar->pml_Exy0[ind_pml] - dExdy + dExdy / kappa);
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] - kpar->dt * (kpar->pml_Ezy0[ind_pml] - dEzdy + dEzdy / kappa);
    }
    else if(j + kpar->j0 >= kpar->pml_ymax){
        ind_pml = i * (kpar->j0 + kpar->J - kpar->pml_ymax) * kpar->K + (kpar->j0 + j - kpar->pml_ymax) * kpar->K + k;
        ind_pml_param = j + kpar->j0 - kpar->pml_ymax + kpar->w_pml_y0;
        kappa = kpar->kappa_H_y[ind_pml_param];
        b = kpar->bHy[ind_pml_param];
        C = kpar->cHy[ind_pml_param];
        kpar->pml_Exy1[ind_pml] = C * dExdy + b * kpar->pml_Exy1[ind_pml];
        kpar->pml_Ezy1[ind_pml] = C * dEzdy + b * kpar->pml_Ezy1[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + kpar->dt * (kpar->pml_Exy1[ind_pml] - dExdy + dExdy / kappa);
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] - kpar->dt * (kpar->pml_Ezy1[ind_pml] - dEzdy + dEzdy / kappa);
    }
    if (i + kpar->i0 < kpar->pml_zmin){
        ind_pml = i * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_zmin - i - kpar->i0 - 1;
        kappa = kpar->kappa_H_z[ind_pml_param];
        b = kpar->bHz[ind_pml_param];
        C = kpar->cHz[ind_pml_param];
        kpar->pml_Exz0[ind_pml] = C * dExdz + b * kpar->pml_Exz0[ind_pml];
        kpar->pml_Eyz0[ind_pml] = C * dEydz + b * kpar->pml_Eyz0[ind_pml];
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + kpar->dt * (kpar->pml_Eyz0[ind_pml] - dEydz + dEydz / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] - kpar->dt * (kpar->pml_Exz0[ind_pml] - dExdz + dExdz / kappa);
    }
    else if(i + kpar->i0 > kpar->pml_zmax){
        ind_pml = (kpar->i0 + i - kpar->pml_zmax) * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = i + kpar->i0 - kpar->pml_zmax + kpar->w_pml_z0;
        kappa = kpar->kappa_H_z[ind_pml_param];
        b = kpar->bHz[ind_pml_param];
        C = kpar->cHz[ind_pml_param];
        kpar->pml_Exz1[ind_pml] = C * dExdz + b * kpar->pml_Exz1[ind_pml];
        kpar->pml_Eyz1[ind_pml] = C * dEydz + b * kpar->pml_Eyz1[ind_pml];
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + kpar->dt * (kpar->pml_Eyz1[ind_pml] - dEydz + dEydz / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] - kpar->dt * (kpar->pml_Exz1[ind_pml] - dExdz + dExdz / kappa);
    }

    // update sources
    // kernel/domain's ind_global = i*_J*_K + j*_K + k;
    srcind_size_offset = 0;
    for(int ii = 0; ii < kpar->srclen; ii ++){
        srcind_size = kpar->Is[ii]*kpar->Js[ii]*kpar->Ks[ii];
        if( i+kpar->i0 >= kpar->i0s[ii] && j+kpar->j0 >= kpar->j0s[ii] && k+kpar->k0 >= kpar->k0s[ii]
            && i+kpar->i0 < kpar->i0s[ii]+kpar->Is[ii] && j+kpar->j0 < kpar->j0s[ii]+kpar->Js[ii] && k+kpar->k0 < kpar->k0s[ii]+kpar->Ks[ii]){
            srcind_src = (i+kpar->i0-kpar->i0s[ii])*kpar->Js[ii]*kpar->Ks[ii] + (j+kpar->j0-kpar->j0s[ii])*kpar->Ks[ii] + k+kpar->k0-kpar->k0s[ii];
            if (kpar->t <= kpar->src_T){
                src_t = sin(kpar->t + kpar->Mx[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + src_t*kpar->Mx[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->My[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + src_t*kpar->My[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->Mz[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + src_t*kpar->Mz[srcind_src+srcind_size_offset].real*kpar->dt;
            }
            else{
                src_t = sin(kpar->t + kpar->Mx[srcind_src+srcind_size_offset].imag);
                kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + src_t*kpar->Mx[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->My[srcind_src+srcind_size_offset].imag);
                kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + src_t*kpar->My[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->Mz[srcind_src+srcind_size_offset].imag);
                kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + src_t*kpar->Mz[srcind_src+srcind_size_offset].real*kpar->dt;
            }
        }
        srcind_size_offset += srcind_size;
    }

    // // store left ghost layers
    // if(kpar->K<kpar->Nx && kpar->k0!=0 && k==0){
    //     kpar->H_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Hx[ind_ijk];
    //     kpar->H_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Hy[ind_ijk];
    //     kpar->H_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Hz[ind_ijk];
    // }
    // if(kpar->J<kpar->Ny && kpar->j0!=0 && j==0){
    //     kpar->H_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Hx[ind_ijk];
    //     kpar->H_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Hy[ind_ijk];
    //     kpar->H_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Hz[ind_ijk];
    // }
    // if(kpar->I<kpar->Nz && kpar->i0!=0 && i==0){
    //     kpar->H_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Hx[ind_ijk];
    //     kpar->H_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Hy[ind_ijk];
    //     kpar->H_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Hz[ind_ijk];
    // }
    // // store right ghost layers
    // if(kpar->K<kpar->Nx && kpar->k0+kpar->K!=kpar->Nx && k==kpar->K-1){
    //     kpar->H_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Hx[ind_ijk];
    //     kpar->H_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Hy[ind_ijk];
    //     kpar->H_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Hz[ind_ijk];
    // }
    // if(kpar->J<kpar->Ny && kpar->j0+kpar->J!=kpar->Ny && j==kpar->J-1){
    //     kpar->H_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Hx[ind_ijk];
    //     kpar->H_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Hy[ind_ijk];
    //     kpar->H_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Hz[ind_ijk];
    // }
    // if(kpar->I<kpar->Nz && kpar->i0+kpar->I!=kpar->Nz && i==kpar->I-1){
    //     kpar->H_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Hx[ind_ijk];
    //     kpar->H_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Hy[ind_ijk];
    //     kpar->H_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Hz[ind_ijk];
    // }
}

__global__ void kernel_update_H_border(kernelpar *kpar)
{
//     1D grid of 1D blocks
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind_global >= kpar->buffer_size) return;

    int i, j, k; 
    if(kpar->K < kpar->Nx){
        i = (ind_global%(kpar->I*kpar->J))/kpar->J;
        j = (ind_global%(kpar->I*kpar->J))%kpar->J;
        if(kpar->k0==0){ k = kpar->K-1; }
        else if(kpar->k0+kpar->K==kpar->Nx){ k=0; }
        else{
            k = (ind_global < kpar->I*kpar->J)? 0 : kpar->K-1;
        }
    }
    else if(kpar->J < kpar->Ny){
        i = (ind_global%(kpar->I*kpar->K))/kpar->K;
        k = (ind_global%(kpar->I*kpar->K))%kpar->K;
        if(kpar->j0==0){ j = kpar->J-1; }
        else if(kpar->j0+kpar->J==kpar->Ny){ j=0; }
        else{
            j = (ind_global < kpar->I*kpar->K)? 0 : kpar->J-1;
        }
    }
    else if(kpar->I < kpar->Nz){
        j = (ind_global%(kpar->J*kpar->K))/kpar->K;
        k = (ind_global%(kpar->J*kpar->K))%kpar->K;
        if(kpar->i0==0){ i = kpar->I-1; }
        else if(kpar->i0+kpar->I==kpar->Nz){ i=0; }
        else{
            i = (ind_global < kpar->J*kpar->K)? 0 : kpar->I-1;
        }
    }
    else{
        // if(ind_global==0) printf("No H border to update.\n");
        return;
    }

    int ind_ijk, ind_ijp1k, ind_ip1jk, ind_ijkp1;
    double dExdy, dExdz, dEydx, dEydz, dEzdx, dEzdy;
    int ind_pml, ind_pml_param;
    double b, C, kappa;
    int srcind_src, srcind_size, srcind_size_offset;
    double src_t;

    ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2) + (j+1)*(kpar->K+2) + k + 1;
    ind_ijp1k = ind_ijk + kpar->K + 2;
    ind_ip1jk = ind_ijk + (kpar->J+2)*(kpar->K+2);
    ind_ijkp1 = ind_ijk + 1;

    // set up fields on the boundary
    if (kpar->bc[0] != 'P' && k == kpar->Nx - kpar->k0){
        kpar->Ey[ind_ijk] = 0.0;
        kpar->Ez[ind_ijk] = 0.0;
    }
    if (kpar->bc[1] != 'P' && j == kpar->Ny - kpar->j0){
        kpar->Ex[ind_ijk] = 0.0;
        kpar->Ez[ind_ijk] = 0.0;
    }
    if (kpar->bc[2] != 'P' && i == kpar->Nz - kpar->i0){
        kpar->Ex[ind_ijk] = 0.0;
        kpar->Ey[ind_ijk] = 0.0;
    }

    dEzdy = kpar->odx * (kpar->Ez[ind_ijp1k] - kpar->Ez[ind_ijk]);
    dEydz = kpar->odx * (kpar->Ey[ind_ip1jk] - kpar->Ey[ind_ijk]);
    dExdz = kpar->odx * (kpar->Ex[ind_ip1jk] - kpar->Ex[ind_ijk]);
    dEzdx = kpar->odx * (kpar->Ez[ind_ijkp1] - kpar->Ez[ind_ijk]);
    dEydx = kpar->odx * (kpar->Ey[ind_ijkp1] - kpar->Ey[ind_ijk]);
    dExdy = kpar->odx * (kpar->Ex[ind_ijp1k] - kpar->Ex[ind_ijk]);
    kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + kpar->dt * (dEydz - dEzdy);
    kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + kpar->dt * (dEzdx - dExdz);
    kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + kpar->dt * (dExdy - dEydx);

    // Update PML
    if (k + kpar->k0 < kpar->pml_xmin){
        ind_pml = i * kpar->J * (kpar->pml_xmin - kpar->k0) + j * (kpar->pml_xmin - kpar->k0) + k;
        ind_pml_param = kpar->pml_xmin - k - kpar->k0 - 1;
        kappa = kpar->kappa_H_x[ind_pml_param];
        b = kpar->bHx[ind_pml_param];
        C = kpar->cHx[ind_pml_param];
        kpar->pml_Eyx0[ind_pml] = C * dEydx + b * kpar->pml_Eyx0[ind_pml];
        kpar->pml_Ezx0[ind_pml] = C * dEzdx + b * kpar->pml_Ezx0[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] - kpar->dt * (kpar->pml_Eyx0[ind_pml] - dEydx + dEydx / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + kpar->dt * (kpar->pml_Ezx0[ind_pml] - dEzdx + dEzdx / kappa);
    }
    else if(k + kpar->k0 >= kpar->pml_xmax){
        ind_pml = i * kpar->J * (kpar->k0 + kpar->K - kpar->pml_xmax) + j * (kpar->k0 + kpar->K - kpar->pml_xmax)
                + k + kpar->k0 - kpar->pml_xmax;
        ind_pml_param = k + kpar->k0 - kpar->pml_xmax + kpar->w_pml_x0;
        kappa = kpar->kappa_H_x[ind_pml_param];
        b = kpar->bHx[ind_pml_param];
        C = kpar->cHx[ind_pml_param];
        kpar->pml_Eyx1[ind_pml] = C * dEydx + b * kpar->pml_Eyx1[ind_pml];
        kpar->pml_Ezx1[ind_pml] = C * dEzdx + b * kpar->pml_Ezx1[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] - kpar->dt * (kpar->pml_Eyx1[ind_pml] - dEydx + dEydx / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + kpar->dt * (kpar->pml_Ezx1[ind_pml] - dEzdx + dEzdx / kappa);
    }
    if (j + kpar->j0 < kpar->pml_ymin){
        ind_pml = i * (kpar->pml_ymin - kpar->j0) * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_ymin - j - kpar->j0 - 1;
        kappa = kpar->kappa_H_y[ind_pml_param];
        b = kpar->bHy[ind_pml_param];
        C = kpar->cHy[ind_pml_param];
        kpar->pml_Exy0[ind_pml] = C * dExdy + b * kpar->pml_Exy0[ind_pml];
        kpar->pml_Ezy0[ind_pml] = C * dEzdy + b * kpar->pml_Ezy0[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + kpar->dt * (kpar->pml_Exy0[ind_pml] - dExdy + dExdy / kappa);
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] - kpar->dt * (kpar->pml_Ezy0[ind_pml] - dEzdy + dEzdy / kappa);
    }
    else if(j + kpar->j0 >= kpar->pml_ymax){
        ind_pml = i * (kpar->j0 + kpar->J - kpar->pml_ymax) * kpar->K + (kpar->j0 + j - kpar->pml_ymax) * kpar->K + k;
        ind_pml_param = j + kpar->j0 - kpar->pml_ymax + kpar->w_pml_y0;
        kappa = kpar->kappa_H_y[ind_pml_param];
        b = kpar->bHy[ind_pml_param];
        C = kpar->cHy[ind_pml_param];
        kpar->pml_Exy1[ind_pml] = C * dExdy + b * kpar->pml_Exy1[ind_pml];
        kpar->pml_Ezy1[ind_pml] = C * dEzdy + b * kpar->pml_Ezy1[ind_pml];
        kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + kpar->dt * (kpar->pml_Exy1[ind_pml] - dExdy + dExdy / kappa);
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] - kpar->dt * (kpar->pml_Ezy1[ind_pml] - dEzdy + dEzdy / kappa);
    }
    if (i + kpar->i0 < kpar->pml_zmin){
        ind_pml = i * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_zmin - i - kpar->i0 - 1;
        kappa = kpar->kappa_H_z[ind_pml_param];
        b = kpar->bHz[ind_pml_param];
        C = kpar->cHz[ind_pml_param];
        kpar->pml_Exz0[ind_pml] = C * dExdz + b * kpar->pml_Exz0[ind_pml];
        kpar->pml_Eyz0[ind_pml] = C * dEydz + b * kpar->pml_Eyz0[ind_pml];
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + kpar->dt * (kpar->pml_Eyz0[ind_pml] - dEydz + dEydz / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] - kpar->dt * (kpar->pml_Exz0[ind_pml] - dExdz + dExdz / kappa);
    }
    else if(i + kpar->i0 > kpar->pml_zmax){
        ind_pml = (kpar->i0 + i - kpar->pml_zmax) * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = i + kpar->i0 - kpar->pml_zmax + kpar->w_pml_z0;
        kappa = kpar->kappa_H_z[ind_pml_param];
        b = kpar->bHz[ind_pml_param];
        C = kpar->cHz[ind_pml_param];
        kpar->pml_Exz1[ind_pml] = C * dExdz + b * kpar->pml_Exz1[ind_pml];
        kpar->pml_Eyz1[ind_pml] = C * dEydz + b * kpar->pml_Eyz1[ind_pml];
        kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + kpar->dt * (kpar->pml_Eyz1[ind_pml] - dEydz + dEydz / kappa);
        kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] - kpar->dt * (kpar->pml_Exz1[ind_pml] - dExdz + dExdz / kappa);
    }

    // update sources
    // kernel/domain's ind_global = i*_J*_K + j*_K + k;
    srcind_size_offset = 0;
    for(int ii = 0; ii < kpar->srclen; ii ++){
        srcind_size = kpar->Is[ii]*kpar->Js[ii]*kpar->Ks[ii];
        if( i+kpar->i0 >= kpar->i0s[ii] && j+kpar->j0 >= kpar->j0s[ii] && k+kpar->k0 >= kpar->k0s[ii]
            && i+kpar->i0 < kpar->i0s[ii]+kpar->Is[ii] && j+kpar->j0 < kpar->j0s[ii]+kpar->Js[ii] && k+kpar->k0 < kpar->k0s[ii]+kpar->Ks[ii]){
            srcind_src = (i+kpar->i0-kpar->i0s[ii])*kpar->Js[ii]*kpar->Ks[ii] + (j+kpar->j0-kpar->j0s[ii])*kpar->Ks[ii] + k+kpar->k0-kpar->k0s[ii];
            if (kpar->t <= kpar->src_T){
                src_t = sin(kpar->t + kpar->Mx[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + src_t*kpar->Mx[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->My[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + src_t*kpar->My[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->Mz[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + src_t*kpar->Mz[srcind_src+srcind_size_offset].real*kpar->dt;
            }
            else{
                src_t = sin(kpar->t + kpar->Mx[srcind_src+srcind_size_offset].imag);
                kpar->Hx[ind_ijk] = kpar->Hx[ind_ijk] + src_t*kpar->Mx[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->My[srcind_src+srcind_size_offset].imag);
                kpar->Hy[ind_ijk] = kpar->Hy[ind_ijk] + src_t*kpar->My[srcind_src+srcind_size_offset].real*kpar->dt;
                src_t = sin(kpar->t + kpar->Mz[srcind_src+srcind_size_offset].imag);
                kpar->Hz[ind_ijk] = kpar->Hz[ind_ijk] + src_t*kpar->Mz[srcind_src+srcind_size_offset].real*kpar->dt;
            }
        }
        srcind_size_offset += srcind_size;
    }

    // store left ghost layers
    if(kpar->K<kpar->Nx && kpar->k0!=0 && k==0){
        kpar->H_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Hx[ind_ijk];
        kpar->H_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Hy[ind_ijk];
        kpar->H_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Hz[ind_ijk];
    }
    if(kpar->J<kpar->Ny && kpar->j0!=0 && j==0){
        kpar->H_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Hx[ind_ijk];
        kpar->H_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Hy[ind_ijk];
        kpar->H_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Hz[ind_ijk];
    }
    if(kpar->I<kpar->Nz && kpar->i0!=0 && i==0){
        kpar->H_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Hx[ind_ijk];
        kpar->H_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Hy[ind_ijk];
        kpar->H_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Hz[ind_ijk];
    }
    // store right ghost layers
    if(kpar->K<kpar->Nx && kpar->k0+kpar->K!=kpar->Nx && k==kpar->K-1){
        kpar->H_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Hx[ind_ijk];
        kpar->H_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Hy[ind_ijk];
        kpar->H_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Hz[ind_ijk];
    }
    if(kpar->J<kpar->Ny && kpar->j0+kpar->J!=kpar->Ny && j==kpar->J-1){
        kpar->H_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Hx[ind_ijk];
        kpar->H_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Hy[ind_ijk];
        kpar->H_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Hz[ind_ijk];
    }
    if(kpar->I<kpar->Nz && kpar->i0+kpar->I!=kpar->Nz && i==kpar->I-1){
        kpar->H_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Hx[ind_ijk];
        kpar->H_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Hy[ind_ijk];
        kpar->H_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Hz[ind_ijk];
    }
}

__global__ void kernel_load_ghostlayerH(kernelpar *kpar)
{
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind_global >= kpar->buffer_size) return;
    
    int i, j, k, ind_ijk;
    if(kpar->K<kpar->Nx){
        i = (ind_global%(kpar->I*kpar->J))/kpar->J;
        j = (ind_global%(kpar->I*kpar->J))%kpar->J;
        if(kpar->k0==0){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+kpar->K+1;
            kpar->Hx[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
            kpar->Hy[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
            kpar->Hz[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
        }
        else if(kpar->k0+kpar->K==kpar->Nx){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+0;
            kpar->Hx[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
            kpar->Hy[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
            kpar->Hz[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
        }
        else{
            if(ind_global < kpar->I*kpar->J){
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+kpar->K+1;
                kpar->Hx[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
                kpar->Hy[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
                kpar->Hz[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
            }
            else{
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+0;
                kpar->Hx[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
                kpar->Hy[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
                kpar->Hz[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
            }
        }
    }
    if(kpar->J<kpar->Ny){
        i = (ind_global%(kpar->I*kpar->K))/kpar->K;
        k = (ind_global%(kpar->I*kpar->K))%kpar->K;
        if(kpar->j0==0){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(kpar->J+1)*(kpar->K+2)+k+1;
            kpar->Hx[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
            kpar->Hy[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
            kpar->Hz[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
        }
        else if(kpar->j0+kpar->J==kpar->Ny){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(0)*(kpar->K+2)+k+1;
            kpar->Hx[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
            kpar->Hy[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
            kpar->Hz[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
        }
        else{
            if(ind_global < kpar->I*kpar->K){
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(kpar->J+1)*(kpar->K+2)+k+1;
                kpar->Hx[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
                kpar->Hy[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
                kpar->Hz[ind_ijk] = kpar->H_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
            }
            else{
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(0)*(kpar->K+2)+k+1;
                kpar->Hx[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
                kpar->Hy[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
                kpar->Hz[ind_ijk] = kpar->H_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
            }
        }
    }
    if(kpar->I<kpar->Nz){
        j = (ind_global%(kpar->J*kpar->K))/kpar->K;
        k = (ind_global%(kpar->J*kpar->K))%kpar->K;
        if(kpar->i0==0){
            ind_ijk = (kpar->I+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
            kpar->Hx[ind_ijk] = kpar->H_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
            kpar->Hy[ind_ijk] = kpar->H_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
            kpar->Hz[ind_ijk] = kpar->H_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
        }
        else if(kpar->i0+kpar->I==kpar->Nz){
            ind_ijk = (0)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
            kpar->Hx[ind_ijk] = kpar->H_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
            kpar->Hy[ind_ijk] = kpar->H_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
            kpar->Hz[ind_ijk] = kpar->H_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
        }
        else{
            if(ind_global < kpar->J*kpar->K){
                ind_ijk = (kpar->I+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
                kpar->Hx[ind_ijk] = kpar->H_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
                kpar->Hy[ind_ijk] = kpar->H_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
                kpar->Hz[ind_ijk] = kpar->H_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
            }
            else{
                ind_ijk = (0)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
                kpar->Hx[ind_ijk] = kpar->H_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
                kpar->Hy[ind_ijk] = kpar->H_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
                kpar->Hz[ind_ijk] = kpar->H_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
            }
        }
    }
}

__global__ void kernel_load_ghostlayerE(kernelpar *kpar)
{
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind_global >= kpar->buffer_size) return;

    int i, j, k, ind_ijk;
    if(kpar->K<kpar->Nx){
        i = (ind_global%(kpar->I*kpar->J))/kpar->J;
        j = (ind_global%(kpar->I*kpar->J))%kpar->J;
        if(kpar->k0==0){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+kpar->K+1;
            kpar->Ex[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
            kpar->Ey[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
            kpar->Ez[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
        }
        else if(kpar->k0+kpar->K==kpar->Nx){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+0;
            kpar->Ex[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
            kpar->Ey[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
            kpar->Ez[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
        }
        else{
            if(ind_global < kpar->I*kpar->J){
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+0;
                kpar->Ex[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
                kpar->Ey[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
                kpar->Ez[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
            }
            else{
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+kpar->K+1;
                kpar->Ex[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0];
                kpar->Ey[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1];
                kpar->Ez[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2];
            }
        }
    }
    if(kpar->J<kpar->Ny){
        i = (ind_global%(kpar->I*kpar->K))/kpar->K;
        k = (ind_global%(kpar->I*kpar->K))%kpar->K;
        if(kpar->j0==0){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(kpar->J+1)*(kpar->K+2)+k+1;
            kpar->Ex[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
            kpar->Ey[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
            kpar->Ez[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
        }
        else if(kpar->j0+kpar->J==kpar->Ny){
            ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(0)*(kpar->K+2)+k+1;
            kpar->Ex[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
            kpar->Ey[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
            kpar->Ez[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
        }
        else{
            if(ind_global < kpar->I*kpar->K){
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(0)*(kpar->K+2)+k+1;
                kpar->Ex[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
                kpar->Ey[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
                kpar->Ez[ind_ijk] = kpar->E_bufferleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
            }
            else{
                ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2)+(kpar->J+1)*(kpar->K+2)+k+1;
                kpar->Ex[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0];
                kpar->Ey[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1];
                kpar->Ez[ind_ijk] = kpar->E_bufferright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2];
            }
        }
    }
    if(kpar->I<kpar->Nz){
        j = (ind_global%(kpar->J*kpar->K))/kpar->K;
        k = (ind_global%(kpar->J*kpar->K))%kpar->K;
        if(kpar->i0==0){
            ind_ijk = (kpar->I+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
            kpar->Ex[ind_ijk] = kpar->E_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
            kpar->Ey[ind_ijk] = kpar->E_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
            kpar->Ez[ind_ijk] = kpar->E_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
        }
        else if(kpar->i0+kpar->I==kpar->Nz){
            ind_ijk = (0)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
            kpar->Ex[ind_ijk] = kpar->E_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
            kpar->Ey[ind_ijk] = kpar->E_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
            kpar->Ez[ind_ijk] = kpar->E_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
        }
        else{
            if(ind_global < kpar->J*kpar->K){
                ind_ijk = (0)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
                kpar->Ex[ind_ijk] = kpar->E_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
                kpar->Ey[ind_ijk] = kpar->E_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
                kpar->Ez[ind_ijk] = kpar->E_bufferleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
            }
            else{
                ind_ijk = (kpar->I+1)*(kpar->J+2)*(kpar->K+2)+(j+1)*(kpar->K+2)+k+1;
                kpar->Ex[ind_ijk] = kpar->E_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0];
                kpar->Ey[ind_ijk] = kpar->E_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1];
                kpar->Ez[ind_ijk] = kpar->E_bufferright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2];
            }
        }
    }
}

__global__ void kernel_update_E_bulk(kernelpar *kpar)
{
//     1D grid of 1D blocks
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    // multiple GPU version: decompose along X
    // size_t ind_global = (i-kpar->i0) * kpar->K * kpar->J + (j-kpar->j0) * kpar->K + (k-kpar->k0);
//     1D grid of 2D blocks
//     size_t ind_global = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
//     2D grid of 1D blocks
//     size_t ind_global = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;

    if (ind_global >= kpar->size) return;

    int i, j, k, ind_ijk, ind_ijm1k, ind_im1jk, ind_ijkm1;
    double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy, b_x, b_y, b_z;
    int ind_pml, ind_pml_param;
    double b, C, kappa;
    int srcind_src, srcind_size, srcind_size_offset;
    double src_t;

    i = ind_global/(kpar->J*kpar->K);
    j = (ind_global%(kpar->J*kpar->K))/kpar->K;
    k = (ind_global%(kpar->J*kpar->K))%kpar->K;

    if(kpar->K<kpar->Nx && kpar->k0!=0 && k==0) return;
    if(kpar->K<kpar->Nx && kpar->k0+kpar->K!=kpar->Nx && k==kpar->K-1) return;
    if(kpar->J<kpar->Ny && kpar->j0!=0 && j==0) return;
    if(kpar->J<kpar->Ny && kpar->j0+kpar->J!=kpar->Ny && j==kpar->J-1) return;
    if(kpar->I<kpar->Nz && kpar->i0!=0 && i==0) return;
    if(kpar->I<kpar->Nz && kpar->i0+kpar->I!=kpar->Nz && i==kpar->I-1) return;

    ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2) + (j+1)*(kpar->K+2) + k + 1;
    ind_ijm1k = ind_ijk - kpar->K - 2;
    ind_im1jk = ind_ijk - (kpar->J+2)*(kpar->K+2);
    ind_ijkm1 = ind_ijk - 1;

    // set up fields on the boundary
    if(kpar->k0==0 && k==0){
        if(kpar->bc[0]=='0'){
            kpar->Hy[ind_ijk-1] = 0.0;
            kpar->Hz[ind_ijk-1] = 0.0;
        }
        else if(kpar->bc[0]=='E'){
            kpar->Hy[ind_ijk-1] = -1.0*kpar->Hy[ind_ijk];
            kpar->Hz[ind_ijk-1] = -1.0*kpar->Hz[ind_ijk];
        }
        else if(kpar->bc[0]=='H'){
            kpar->Hy[ind_ijk-1] = kpar->Hy[ind_ijk];
            kpar->Hz[ind_ijk-1] = kpar->Hz[ind_ijk];
        }
    }
    if(kpar->j0==0 && j==0){
        if(kpar->bc[1]=='0'){
            kpar->Hx[ind_ijk-kpar->K-2] = 0.0;
            kpar->Hz[ind_ijk-kpar->K-2] = 0.0;
        }
        else if(kpar->bc[1]=='E'){
            kpar->Hx[ind_ijk-kpar->K-2] = -1.0*kpar->Hx[ind_ijk];
            kpar->Hz[ind_ijk-kpar->K-2] = -1.0*kpar->Hz[ind_ijk];
        }
        else if(kpar->bc[1]=='H'){
            kpar->Hx[ind_ijk-kpar->K-2] = kpar->Hx[ind_ijk];
            kpar->Hz[ind_ijk-kpar->K-2] = kpar->Hz[ind_ijk];
        }
    }
    if(kpar->i0==0 && i==0){
        if(kpar->bc[2]=='0'){
            kpar->Hx[ind_ijk-(kpar->J+2)*(kpar->K+2)] = 0.0;
            kpar->Hy[ind_ijk-(kpar->J+2)*(kpar->K+2)] = 0.0;
        }
        else if(kpar->bc[2]=='E'){
            kpar->Hx[ind_ijk-(kpar->J+2)*(kpar->K+2)] = -1.0*kpar->Hx[ind_ijk];
            kpar->Hy[ind_ijk-(kpar->J+2)*(kpar->K+2)] = -1.0*kpar->Hy[ind_ijk];
        }
        else if(kpar->bc[2]=='H'){
            kpar->Hx[ind_ijk-(kpar->J+2)*(kpar->K+2)] = kpar->Hx[ind_ijk];
            kpar->Hy[ind_ijk-(kpar->J+2)*(kpar->K+2)] = kpar->Hy[ind_ijk];
        }
    }

    // update fields
    b_x = kpar->dt/kpar->epsx[ind_global].real;
    b_y = kpar->dt/kpar->epsy[ind_global].real;
    b_z = kpar->dt/kpar->epsz[ind_global].real;
    dHzdy = kpar->odx * (kpar->Hz[ind_ijk] - kpar->Hz[ind_ijm1k]);
    dHydz = kpar->odx * (kpar->Hy[ind_ijk] - kpar->Hy[ind_im1jk]);
    dHxdz = kpar->odx * (kpar->Hx[ind_ijk] - kpar->Hx[ind_im1jk]);
    dHzdx = kpar->odx * (kpar->Hz[ind_ijk] - kpar->Hz[ind_ijkm1]);
    dHydx = kpar->odx * (kpar->Hy[ind_ijk] - kpar->Hy[ind_ijkm1]);
    dHxdy = kpar->odx * (kpar->Hx[ind_ijk] - kpar->Hx[ind_ijm1k]);
    kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] + b_x * (dHzdy - dHydz);
    kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] + b_y * (dHxdz - dHzdx);
    kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] + b_z * (dHydx - dHxdy);

    // Update PML
    if (k + kpar->k0 < kpar->pml_xmin){
        ind_pml = i * kpar->J * (kpar->pml_xmin - kpar->k0) + j * (kpar->pml_xmin - kpar->k0) + k;
        ind_pml_param = kpar->pml_xmin - k - kpar->k0 - 1;
        kappa = kpar->kappa_E_x[ind_pml_param];
        b = kpar->bEx[ind_pml_param];
        C = kpar->cEx[ind_pml_param];
        kpar->pml_Hyx0[ind_pml] = C * dHydx + b * kpar->pml_Hyx0[ind_pml];
        kpar->pml_Hzx0[ind_pml] = C * dHzdx + b * kpar->pml_Hzx0[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] + b_z * (kpar->pml_Hyx0[ind_pml] - dHydx + dHydx / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - b_y * (kpar->pml_Hzx0[ind_pml] - dHzdx + dHzdx / kappa);
    }
    else if(k + kpar->k0 >= kpar->pml_xmax){
        ind_pml = i * kpar->J * (kpar->k0 + kpar->K - kpar->pml_xmax) + j * (kpar->k0 + kpar->K - kpar->pml_xmax)
                + k + kpar->k0 - kpar->pml_xmax;
        ind_pml_param = k + kpar->k0 - kpar->pml_xmax + kpar->w_pml_x0;
        kappa = kpar->kappa_E_x[ind_pml_param];
        b = kpar->bEx[ind_pml_param];
        C = kpar->cEx[ind_pml_param];
        kpar->pml_Hyx1[ind_pml] = C * dHydx + b * kpar->pml_Hyx1[ind_pml];
        kpar->pml_Hzx1[ind_pml] = C * dHzdx + b * kpar->pml_Hzx1[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] + b_z * (kpar->pml_Hyx1[ind_pml] - dHydx + dHydx / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - b_y * (kpar->pml_Hzx1[ind_pml] - dHzdx + dHzdx / kappa);
    }

    if (j + kpar->j0 < kpar->pml_ymin){
        ind_pml = i * (kpar->pml_ymin - kpar->j0) * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_ymin - j - kpar->j0 - 1;
        kappa = kpar->kappa_E_y[ind_pml_param];
        b = kpar->bEy[ind_pml_param];
        C = kpar->cEy[ind_pml_param];
        kpar->pml_Hxy0[ind_pml] = C * dHxdy + b * kpar->pml_Hxy0[ind_pml];
        kpar->pml_Hzy0[ind_pml] = C * dHzdy + b * kpar->pml_Hzy0[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - b_z * (kpar->pml_Hxy0[ind_pml] - dHxdy + dHxdy / kappa);
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] + b_x * (kpar->pml_Hzy0[ind_pml] - dHzdy + dHzdy / kappa);
    }
    else if(j + kpar->j0 >= kpar->pml_ymax){
        ind_pml = i * (kpar->j0 + kpar->J - kpar->pml_ymax) * kpar->K + (kpar->j0 + j - kpar->pml_ymax) * kpar->K + k;
        ind_pml_param = j + kpar->j0 - kpar->pml_ymax + kpar->w_pml_y0;
        kappa = kpar->kappa_E_y[ind_pml_param];
        b = kpar->bEy[ind_pml_param];
        C = kpar->cEy[ind_pml_param];
        kpar->pml_Hxy1[ind_pml] = C * dHxdy + b * kpar->pml_Hxy1[ind_pml];
        kpar->pml_Hzy1[ind_pml] = C * dHzdy + b * kpar->pml_Hzy1[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - b_z * (kpar->pml_Hxy1[ind_pml] - dHxdy + dHxdy / kappa);
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] + b_x * (kpar->pml_Hzy1[ind_pml] - dHzdy + dHzdy / kappa);
    }

    if (i + kpar->i0 < kpar->pml_zmin){
        ind_pml = i * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_zmin - i - kpar->i0 - 1;
        kappa = kpar->kappa_E_z[ind_pml_param];
        b = kpar->bEz[ind_pml_param];
        C = kpar->cEz[ind_pml_param];
        kpar->pml_Hxz0[ind_pml] = C * dHxdz + b * kpar->pml_Hxz0[ind_pml];
        kpar->pml_Hyz0[ind_pml] = C * dHydz + b * kpar->pml_Hyz0[ind_pml];
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - b_x * (kpar->pml_Hyz0[ind_pml] - dHydz + dHydz / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] + b_y * (kpar->pml_Hxz0[ind_pml] - dHxdz + dHxdz / kappa);
    }
    else if(i + kpar->i0 > kpar->pml_zmax){
        ind_pml = (kpar->i0 + i - kpar->pml_zmax) * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = i + kpar->i0 - kpar->pml_zmax + kpar->w_pml_z0;
        kappa = kpar->kappa_E_z[ind_pml_param];
        b = kpar->bEz[ind_pml_param];
        C = kpar->cEz[ind_pml_param];
        kpar->pml_Hxz1[ind_pml] = C * dHxdz + b * kpar->pml_Hxz1[ind_pml];
        kpar->pml_Hyz1[ind_pml] = C * dHydz + b * kpar->pml_Hyz1[ind_pml];
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - b_x * (kpar->pml_Hyz1[ind_pml] - dHydz + dHydz / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] + b_y * (kpar->pml_Hxz1[ind_pml] - dHxdz + dHxdz / kappa);
    }

    // update sources
    srcind_size_offset = 0;
    for(int ii = 0; ii < kpar->srclen; ii++){
        srcind_size = kpar->Is[ii] * kpar->Js[ii] * kpar->Ks[ii];
        if(i+kpar->i0 >= kpar->i0s[ii] && j+kpar->j0 >= kpar->j0s[ii] && k+kpar->k0 >= kpar->k0s[ii]
           && i+kpar->i0 < kpar->i0s[ii]+kpar->Is[ii] && j+kpar->j0 < kpar->j0s[ii]+kpar->Js[ii] && k+kpar->k0 < kpar->k0s[ii]+kpar->Ks[ii]){
            srcind_src = (i+kpar->i0-kpar->i0s[ii])*kpar->Js[ii]*kpar->Ks[ii] + (j+kpar->j0-kpar->j0s[ii])*kpar->Ks[ii] + k+kpar->k0- kpar->k0s[ii];
            if (kpar->t <= kpar->src_T){
                src_t = sin(kpar->t + kpar->Jx[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - src_t*kpar->Jx[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsx[ind_global].real;
                src_t = sin(kpar->t + kpar->Jy[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - src_t*kpar->Jy[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsy[ind_global].real;
                src_t = sin(kpar->t + kpar->Jz[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - src_t*kpar->Jz[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsz[ind_global].real;
            }
            else{
                src_t = sin(kpar->t + kpar->Jx[srcind_src+srcind_size_offset].imag);
                kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - src_t*kpar->Jx[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsx[ind_global].real;
                src_t = sin(kpar->t + kpar->Jy[srcind_src+srcind_size_offset].imag);
                kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - src_t*kpar->Jy[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsy[ind_global].real;
                src_t = sin(kpar->t + kpar->Jz[srcind_src+srcind_size_offset].imag);
                kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - src_t*kpar->Jz[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsz[ind_global].real;
           }
        }
        srcind_size_offset += srcind_size;
    }

    // // store left ghost layers
    // if(kpar->K<kpar->Nx && kpar->k0!=0 && k==0){
    //     kpar->E_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Ex[ind_ijk];
    //     kpar->E_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Ey[ind_ijk];
    //     kpar->E_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Ez[ind_ijk];
    // }
    // if(kpar->J<kpar->Ny && kpar->j0!=0 && j==0){
    //     kpar->E_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Ex[ind_ijk];
    //     kpar->E_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Ey[ind_ijk];
    //     kpar->E_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Ez[ind_ijk];
    // }
    // if(kpar->I<kpar->Nz && kpar->i0!=0 && i==0){
    //     kpar->E_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Ex[ind_ijk];
    //     kpar->E_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Ey[ind_ijk];
    //     kpar->E_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Ez[ind_ijk];
    // }
    // // store right ghost layers
    // if(kpar->K<kpar->Nx && kpar->k0+kpar->K!=kpar->Nx && k==kpar->K-1){
    //     kpar->E_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Ex[ind_ijk];
    //     kpar->E_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Ey[ind_ijk];
    //     kpar->E_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Ez[ind_ijk];
    // }
    // if(kpar->J<kpar->Ny && kpar->j0+kpar->J!=kpar->Ny && j==kpar->J-1){
    //     kpar->E_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Ex[ind_ijk];
    //     kpar->E_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Ey[ind_ijk];
    //     kpar->E_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Ez[ind_ijk];
    // }
    // if(kpar->I<kpar->Nz && kpar->i0+kpar->I!=kpar->Nz && i==kpar->I-1){
    //     kpar->E_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Ex[ind_ijk];
    //     kpar->E_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Ey[ind_ijk];
    //     kpar->E_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Ez[ind_ijk];
    // }
}

__global__ void kernel_update_E_border(kernelpar *kpar)
{
//     1D grid of 1D blocks
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind_global >= kpar->buffer_size) return;

    int i, j, k;
    if(kpar->K < kpar->Nx){
        i = (ind_global%(kpar->I*kpar->J))/kpar->J;
        j = (ind_global%(kpar->I*kpar->J))%kpar->J;
        if(kpar->k0==0){ k = kpar->K-1; }
        else if(kpar->k0+kpar->K==kpar->Nx){ k=0; }
        else{
            k = (ind_global < kpar->I*kpar->J) ? 0 : kpar->K-1;
        }
    }
    else if(kpar->J < kpar->Ny){
        i = (ind_global%(kpar->I*kpar->K))/kpar->K;
        k = (ind_global%(kpar->I*kpar->K))%kpar->K;
        if(kpar->j0==0){ j = kpar->J-1; }
        else if(kpar->j0+kpar->J==kpar->Ny){ j=0; }
        else{
            j = (ind_global < kpar->I*kpar->K) ? 0 : kpar->J-1;
        }
    }
    else if(kpar->I < kpar->Nz){
        j = (ind_global%(kpar->J*kpar->K))/kpar->K;
        k = (ind_global%(kpar->J*kpar->K))%kpar->K;
        if(kpar->i0==0){ i = kpar->I-1; }
        else if(kpar->i0+kpar->I==kpar->Nz){ i=0; }
        else{
            i = (ind_global < kpar->J*kpar->K) ? 0 : kpar->I-1;
        }
    }
    else{
        // if(ind_global==0) printf("No E border to update.\n");
        return;
    }
    
    int ind_ijk, ind_ijm1k, ind_im1jk, ind_ijkm1;
    double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy, b_x, b_y, b_z;
    int ind_pml, ind_pml_param;
    double b, C, kappa;
    int srcind_src, srcind_size, srcind_size_offset;
    double src_t;

    ind_ijk = (i+1)*(kpar->J+2)*(kpar->K+2) + (j+1)*(kpar->K+2) + k + 1;
    ind_ijm1k = ind_ijk - kpar->K - 2;
    ind_im1jk = ind_ijk - (kpar->J+2)*(kpar->K+2);
    ind_ijkm1 = ind_ijk - 1;

    // set up fields on the boundary
    if(kpar->k0==0 && k==0){
        if(kpar->bc[0]=='0'){
            kpar->Hy[ind_ijk-1] = 0.0;
            kpar->Hz[ind_ijk-1] = 0.0;
        }
        else if(kpar->bc[0]=='E'){
            kpar->Hy[ind_ijk-1] = -1.0*kpar->Hy[ind_ijk];
            kpar->Hz[ind_ijk-1] = -1.0*kpar->Hz[ind_ijk];
        }
        else if(kpar->bc[0]=='H'){
            kpar->Hy[ind_ijk-1] = kpar->Hy[ind_ijk];
            kpar->Hz[ind_ijk-1] = kpar->Hz[ind_ijk];
        }
    }
    if(kpar->j0==0 && j==0){
        if(kpar->bc[1]=='0'){
            kpar->Hx[ind_ijk-kpar->K-2] = 0.0;
            kpar->Hz[ind_ijk-kpar->K-2] = 0.0;
        }
        else if(kpar->bc[1]=='E'){
            kpar->Hx[ind_ijk-kpar->K-2] = -1.0*kpar->Hx[ind_ijk];
            kpar->Hz[ind_ijk-kpar->K-2] = -1.0*kpar->Hz[ind_ijk];
        }
        else if(kpar->bc[1]=='H'){
            kpar->Hx[ind_ijk-kpar->K-2] = kpar->Hx[ind_ijk];
            kpar->Hz[ind_ijk-kpar->K-2] = kpar->Hz[ind_ijk];
        }
    }
    if(kpar->i0==0 && i==0){
        if(kpar->bc[2]=='0'){
            kpar->Hx[ind_ijk-(kpar->J+2)*(kpar->K+2)] = 0.0;
            kpar->Hy[ind_ijk-(kpar->J+2)*(kpar->K+2)] = 0.0;
        }
        else if(kpar->bc[2]=='E'){
            kpar->Hx[ind_ijk-(kpar->J+2)*(kpar->K+2)] = -1.0*kpar->Hx[ind_ijk];
            kpar->Hy[ind_ijk-(kpar->J+2)*(kpar->K+2)] = -1.0*kpar->Hy[ind_ijk];
        }
        else if(kpar->bc[2]=='H'){
            kpar->Hx[ind_ijk-(kpar->J+2)*(kpar->K+2)] = kpar->Hx[ind_ijk];
            kpar->Hy[ind_ijk-(kpar->J+2)*(kpar->K+2)] = kpar->Hy[ind_ijk];
        }
    }

    // update fields
    b_x = kpar->dt/kpar->epsx[ind_global].real;
    b_y = kpar->dt/kpar->epsy[ind_global].real;
    b_z = kpar->dt/kpar->epsz[ind_global].real;
    dHzdy = kpar->odx * (kpar->Hz[ind_ijk] - kpar->Hz[ind_ijm1k]);
    dHydz = kpar->odx * (kpar->Hy[ind_ijk] - kpar->Hy[ind_im1jk]);
    dHxdz = kpar->odx * (kpar->Hx[ind_ijk] - kpar->Hx[ind_im1jk]);
    dHzdx = kpar->odx * (kpar->Hz[ind_ijk] - kpar->Hz[ind_ijkm1]);
    dHydx = kpar->odx * (kpar->Hy[ind_ijk] - kpar->Hy[ind_ijkm1]);
    dHxdy = kpar->odx * (kpar->Hx[ind_ijk] - kpar->Hx[ind_ijm1k]);
    kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] + b_x * (dHzdy - dHydz);
    kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] + b_y * (dHxdz - dHzdx);
    kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] + b_z * (dHydx - dHxdy);

    // Update PML
    if (k + kpar->k0 < kpar->pml_xmin){
        ind_pml = i * kpar->J * (kpar->pml_xmin - kpar->k0) + j * (kpar->pml_xmin - kpar->k0) + k;
        ind_pml_param = kpar->pml_xmin - k - kpar->k0 - 1;
        kappa = kpar->kappa_E_x[ind_pml_param];
        b = kpar->bEx[ind_pml_param];
        C = kpar->cEx[ind_pml_param];
        kpar->pml_Hyx0[ind_pml] = C * dHydx + b * kpar->pml_Hyx0[ind_pml];
        kpar->pml_Hzx0[ind_pml] = C * dHzdx + b * kpar->pml_Hzx0[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] + b_z * (kpar->pml_Hyx0[ind_pml] - dHydx + dHydx / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - b_y * (kpar->pml_Hzx0[ind_pml] - dHzdx + dHzdx / kappa);
    }
    else if(k + kpar->k0 >= kpar->pml_xmax){
        ind_pml = i * kpar->J * (kpar->k0 + kpar->K - kpar->pml_xmax) + j * (kpar->k0 + kpar->K - kpar->pml_xmax)
                + k + kpar->k0 - kpar->pml_xmax;
        ind_pml_param = k + kpar->k0 - kpar->pml_xmax + kpar->w_pml_x0;
        kappa = kpar->kappa_E_x[ind_pml_param];
        b = kpar->bEx[ind_pml_param];
        C = kpar->cEx[ind_pml_param];
        kpar->pml_Hyx1[ind_pml] = C * dHydx + b * kpar->pml_Hyx1[ind_pml];
        kpar->pml_Hzx1[ind_pml] = C * dHzdx + b * kpar->pml_Hzx1[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] + b_z * (kpar->pml_Hyx1[ind_pml] - dHydx + dHydx / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - b_y * (kpar->pml_Hzx1[ind_pml] - dHzdx + dHzdx / kappa);
    }

    if (j + kpar->j0 < kpar->pml_ymin){
        ind_pml = i * (kpar->pml_ymin - kpar->j0) * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_ymin - j - kpar->j0 - 1;
        kappa = kpar->kappa_E_y[ind_pml_param];
        b = kpar->bEy[ind_pml_param];
        C = kpar->cEy[ind_pml_param];
        kpar->pml_Hxy0[ind_pml] = C * dHxdy + b * kpar->pml_Hxy0[ind_pml];
        kpar->pml_Hzy0[ind_pml] = C * dHzdy + b * kpar->pml_Hzy0[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - b_z * (kpar->pml_Hxy0[ind_pml] - dHxdy + dHxdy / kappa);
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] + b_x * (kpar->pml_Hzy0[ind_pml] - dHzdy + dHzdy / kappa);
    }
    else if(j + kpar->j0 >= kpar->pml_ymax){
        ind_pml = i * (kpar->j0 + kpar->J - kpar->pml_ymax) * kpar->K + (kpar->j0 + j - kpar->pml_ymax) * kpar->K + k;
        ind_pml_param = j + kpar->j0 - kpar->pml_ymax + kpar->w_pml_y0;
        kappa = kpar->kappa_E_y[ind_pml_param];
        b = kpar->bEy[ind_pml_param];
        C = kpar->cEy[ind_pml_param];
        kpar->pml_Hxy1[ind_pml] = C * dHxdy + b * kpar->pml_Hxy1[ind_pml];
        kpar->pml_Hzy1[ind_pml] = C * dHzdy + b * kpar->pml_Hzy1[ind_pml];
        kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - b_z * (kpar->pml_Hxy1[ind_pml] - dHxdy + dHxdy / kappa);
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] + b_x * (kpar->pml_Hzy1[ind_pml] - dHzdy + dHzdy / kappa);
    }

    if (i + kpar->i0 < kpar->pml_zmin){
        ind_pml = i * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = kpar->pml_zmin - i - kpar->i0 - 1;
        kappa = kpar->kappa_E_z[ind_pml_param];
        b = kpar->bEz[ind_pml_param];
        C = kpar->cEz[ind_pml_param];
        kpar->pml_Hxz0[ind_pml] = C * dHxdz + b * kpar->pml_Hxz0[ind_pml];
        kpar->pml_Hyz0[ind_pml] = C * dHydz + b * kpar->pml_Hyz0[ind_pml];
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - b_x * (kpar->pml_Hyz0[ind_pml] - dHydz + dHydz / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] + b_y * (kpar->pml_Hxz0[ind_pml] - dHxdz + dHxdz / kappa);
    }
    else if(i + kpar->i0 > kpar->pml_zmax){
        ind_pml = (kpar->i0 + i - kpar->pml_zmax) * kpar->J * kpar->K + j * kpar->K + k;
        ind_pml_param = i + kpar->i0 - kpar->pml_zmax + kpar->w_pml_z0;
        kappa = kpar->kappa_E_z[ind_pml_param];
        b = kpar->bEz[ind_pml_param];
        C = kpar->cEz[ind_pml_param];
        kpar->pml_Hxz1[ind_pml] = C * dHxdz + b * kpar->pml_Hxz1[ind_pml];
        kpar->pml_Hyz1[ind_pml] = C * dHydz + b * kpar->pml_Hyz1[ind_pml];
        kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - b_x * (kpar->pml_Hyz1[ind_pml] - dHydz + dHydz / kappa);
        kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] + b_y * (kpar->pml_Hxz1[ind_pml] - dHxdz + dHxdz / kappa);
    }

    // update sources
    srcind_size_offset = 0;
    for(int ii = 0; ii < kpar->srclen; ii++){
        srcind_size = kpar->Is[ii] * kpar->Js[ii] * kpar->Ks[ii];
        if(i+kpar->i0 >= kpar->i0s[ii] && j+kpar->j0 >= kpar->j0s[ii] && k+kpar->k0 >= kpar->k0s[ii]
           && i+kpar->i0 < kpar->i0s[ii]+kpar->Is[ii] && j+kpar->j0 < kpar->j0s[ii]+kpar->Js[ii] && k+kpar->k0 < kpar->k0s[ii]+kpar->Ks[ii]){
            srcind_src = (i+kpar->i0-kpar->i0s[ii])*kpar->Js[ii]*kpar->Ks[ii] + (j+kpar->j0-kpar->j0s[ii])*kpar->Ks[ii] + k+kpar->k0- kpar->k0s[ii];
            if (kpar->t <= kpar->src_T){
                src_t = sin(kpar->t + kpar->Jx[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - src_t*kpar->Jx[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsx[ind_global].real;
                src_t = sin(kpar->t + kpar->Jy[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - src_t*kpar->Jy[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsy[ind_global].real;
                src_t = sin(kpar->t + kpar->Jz[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - src_t*kpar->Jz[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsz[ind_global].real;
            }
            else{
                src_t = sin(kpar->t + kpar->Jx[srcind_src+srcind_size_offset].imag);
                kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - src_t*kpar->Jx[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsx[ind_global].real;
                src_t = sin(kpar->t + kpar->Jy[srcind_src+srcind_size_offset].imag);
                kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - src_t*kpar->Jy[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsy[ind_global].real;
                src_t = sin(kpar->t + kpar->Jz[srcind_src+srcind_size_offset].imag);
                kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - src_t*kpar->Jz[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsz[ind_global].real;
           }
        }
        srcind_size_offset += srcind_size;
    }

    // store left ghost layers
    if(kpar->K<kpar->Nx && kpar->k0!=0 && k==0){
        kpar->E_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Ex[ind_ijk];
        kpar->E_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Ey[ind_ijk];
        kpar->E_ghostleft[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Ez[ind_ijk];
    }
    if(kpar->J<kpar->Ny && kpar->j0!=0 && j==0){
        kpar->E_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Ex[ind_ijk];
        kpar->E_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Ey[ind_ijk];
        kpar->E_ghostleft[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Ez[ind_ijk];
    }
    if(kpar->I<kpar->Nz && kpar->i0!=0 && i==0){
        kpar->E_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Ex[ind_ijk];
        kpar->E_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Ey[ind_ijk];
        kpar->E_ghostleft[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Ez[ind_ijk];
    }
    // store right ghost layers
    if(kpar->K<kpar->Nx && kpar->k0+kpar->K!=kpar->Nx && k==kpar->K-1){
        kpar->E_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*0] = kpar->Ex[ind_ijk];
        kpar->E_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*1] = kpar->Ey[ind_ijk];
        kpar->E_ghostright[(i-kpar->i0)*kpar->J+j-kpar->j0+kpar->I*kpar->J*2] = kpar->Ez[ind_ijk];
    }
    if(kpar->J<kpar->Ny && kpar->j0+kpar->J!=kpar->Ny && j==kpar->J-1){
        kpar->E_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*0] = kpar->Ex[ind_ijk];
        kpar->E_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*1] = kpar->Ey[ind_ijk];
        kpar->E_ghostright[(i-kpar->i0)*kpar->K+k-kpar->k0+kpar->I*kpar->K*2] = kpar->Ez[ind_ijk];
    }
    if(kpar->I<kpar->Nz && kpar->i0+kpar->I!=kpar->Nz && i==kpar->I-1){
        kpar->E_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*0] = kpar->Ex[ind_ijk];
        kpar->E_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*1] = kpar->Ey[ind_ijk];
        kpar->E_ghostright[(j-kpar->j0)*kpar->K+k-kpar->k0+kpar->J*kpar->K*2] = kpar->Ez[ind_ijk];
    }
}

void prepare_nvlink (int device_id, int gpus_count)
{
	// throw_on_error (cudaSetDevice (thread_info.thread_id), __FILE__, __LINE__);
    gpuErrchk(cudaSetDevice(device_id));
	bool nvlink_enabled = true;
	for (int other_device_id=0; other_device_id<gpus_count; other_device_id++)
	{
		if (other_device_id != device_id)
		{
			int can_access_other_device {};
			gpuErrchk(cudaDeviceCanAccessPeer (&can_access_other_device, device_id, other_device_id));
			if (can_access_other_device)
			{
				gpuErrchk(cudaDeviceEnablePeerAccess (other_device_id, 0));
			}
			else
			{
				std::cerr << "Warning in thread " << device_id <<":device" << device_id
					<< " can't access device " << other_device_id << " memory."
					<< " Fall back to normal copy through the host." << std::endl;
				nvlink_enabled = false;
			}
		}
	}
	for (int tid=0; tid<gpus_count; tid++)
	{
		if(tid==device_id)
		{
			if (nvlink_enabled)
				std::cout << "NVLINK enabled on thread " << device_id << std::endl;
		}
	}
}


void ghost_layer_copyH(kernelpar* left, int leftdev, kernelpar* right, int rightdev, size_t ghost_size)
{   
    // kernelpar *right, *left;
    // right = (kernelpar *)malloc(sizeof(kernelpar));
    // left = (kernelpar *)malloc(sizeof(kernelpar));
    // gpuErrchk(cudaMemcpy(right, right_kpar, sizeof(kernelpar), cudaMemcpyDeviceToHost));
    // gpuErrchk(cudaMemcpy(left, left_kpar, sizeof(kernelpar), cudaMemcpyDeviceToHost));

    // int right_I, right_J;
    // right_I = right->I; right_J = right->J;

    // int canAccessPeer{1};
    // gpuErrchk(cudaDeviceCanAccessPeer(&canAccessPeer, rightdev, leftdev));
    // std::cout << "Right can access left:" << canAccessPeer << std::endl;
    // gpuErrchk(cudaDeviceCanAccessPeer(&canAccessPeer, leftdev, rightdev));
    // std::cout << "Left can access right:" << canAccessPeer << std::endl;
    
    // cudaSetDevice(rightdev);
    // // if(canAccessPeer){ gpuErrchk(cudaDeviceEnablePeerAccess(leftdev, 0)); }
    // // else{ 
    //     cudaDeviceDisablePeerAccess(leftdev); 
    //     // }
    // cudaSetDevice(leftdev);
    // // if(canAccessPeer){ gpuErrchk(cudaDeviceEnablePeerAccess(rightdev, 0)); }
    // // else{ 
    //     cudaDeviceDisablePeerAccess(rightdev); 
    //     // }

    if(leftdev<rightdev){
        // copy ghost layers of right blocks to left blocks
        gpuErrchk(cudaMemcpyPeer(left->H_bufferright, leftdev, right->H_ghostleft, rightdev, 3*ghost_size*sizeof(double)));

        // // copy ghost layers of left blocks to right blocks
        // gpuErrchk(cudaMemcpyPeer(right->Hx_bufferleft, rightdev, left->Hx_ghostright, leftdev, right_I*right_J*sizeof(double)));
        // gpuErrchk(cudaMemcpyPeer(right->Hy_bufferleft, rightdev, left->Hy_ghostright, leftdev, right_I*right_J*sizeof(double)));
        // gpuErrchk(cudaMemcpyPeer(right->Hz_bufferleft, rightdev, left->Hz_ghostright, leftdev, right_I*right_J*sizeof(double)));
        gpuErrchk(cudaMemcpyPeer(right->H_bufferleft, rightdev, left->H_ghostright, leftdev, 3*ghost_size*sizeof(double)));
    }
    else{
        // gpuErrchk(cudaMemcpyPeer(left->Hx_bufferleft, leftdev, right->Hx_ghostright, rightdev, right_I*right_J*sizeof(double)));
        // gpuErrchk(cudaMemcpyPeer(left->Hy_bufferleft, leftdev, right->Hy_ghostright, rightdev, right_I*right_J*sizeof(double)));
        // gpuErrchk(cudaMemcpyPeer(left->Hz_bufferleft, leftdev, right->Hz_ghostright, rightdev, right_I*right_J*sizeof(double)));
        gpuErrchk(cudaMemcpyPeer(left->H_bufferleft, leftdev, right->H_ghostright, rightdev, 3*ghost_size*sizeof(double)));
    }
    // delete[] right; delete[] left;
}

void ghost_layer_copyE(kernelpar* left, int leftdev, kernelpar* right, int rightdev, size_t ghost_size)
{   
    // kernelpar *right, *left;
    // right = (kernelpar *)malloc(sizeof(kernelpar));
    // left = (kernelpar *)malloc(sizeof(kernelpar));
    // gpuErrchk(cudaMemcpy(right, right_kpar, sizeof(kernelpar), cudaMemcpyDeviceToHost));
    // gpuErrchk(cudaMemcpy(left, left_kpar, sizeof(kernelpar), cudaMemcpyDeviceToHost));

    // int right_I, right_J;
    // right_I = right->I; right_J = right->J;

    // int canAccessPeer{1};
    // gpuErrchk(cudaDeviceCanAccessPeer(&canAccessPeer, rightdev, leftdev));
    // std::cout << "Right can access left:" << canAccessPeer << std::endl;
    // gpuErrchk(cudaDeviceCanAccessPeer(&canAccessPeer, leftdev, rightdev));
    // std::cout << "Left can access right:" << canAccessPeer << std::endl;
    
    // cudaSetDevice(rightdev);
    // // if(canAccessPeer){ gpuErrchk(cudaDeviceEnablePeerAccess(leftdev, 0)); }
    // // else{ 
    //     cudaDeviceDisablePeerAccess(leftdev); 
    //     // }
    // cudaSetDevice(leftdev);
    // // if(canAccessPeer){ gpuErrchk(cudaDeviceEnablePeerAccess(rightdev, 0)); }
    // // else{ 
    //     cudaDeviceDisablePeerAccess(rightdev); 
    //     // }

    if(leftdev<rightdev){
        // copy ghost layers of right blocks to left blocks
        gpuErrchk(cudaMemcpyPeer(left->E_bufferright, leftdev, right->E_ghostleft, rightdev, 3*ghost_size*sizeof(double)));
        // copy ghost layers of left blocks to right blocks
        // gpuErrchk(cudaMemcpyPeer(right->Ex_bufferleft, rightdev, left->Ex_ghostright, leftdev, right_I*right_J*sizeof(double)));
        // gpuErrchk(cudaMemcpyPeer(right->Ey_bufferleft, rightdev, left->Ey_ghostright, leftdev, right_I*right_J*sizeof(double)));
        // gpuErrchk(cudaMemcpyPeer(right->Ez_bufferleft, rightdev, left->Ez_ghostright, leftdev, right_I*right_J*sizeof(double)));
        // gpuErrchk(cudaMemcpyPeer(right->E_bufferleft, rightdev, left->E_ghostright, leftdev, 3*ghost_size*sizeof(double)));
    }
    else{
        gpuErrchk(cudaMemcpyPeer(left->E_bufferleft, leftdev, right->E_ghostright, rightdev, 3*ghost_size*sizeof(double)));
    }
    // delete[] right; delete[] left;
}

int main(int argc, char* argv[])
{
    double t{};
    int i0{}, j0{}, k0{}, I{}, J{}, K{}, Nx{}, Ny{}, Nz{};
    int pml_xmin{}, pml_xmax{}, pml_ymin{}, pml_ymax{}, pml_zmin{}, pml_zmax{};
    int w_pml_x0{}, w_pml_x1{}, w_pml_y0{}, w_pml_y1{}, w_pml_z0{}, w_pml_z1{};
    size_t size{};
    char bc[]="000"; //int bc0{}, bc1{}, bc2{};
    double odx{}, dt{};
    //double src_t{};
    int *i0s, *j0s, *k0s, *Is, *Js, *Ks, srclen{};
    double src_T{}, src_min{}, src_k{};
    complex128 *epsx, *epsy, *epsz;
    complex128 *Jx, *Jy, *Jz, *Mx, *My, *Mz;

    double *pml_Exy0, *pml_Exy1, *pml_Exz0, *pml_Exz1,
           *pml_Eyx0, *pml_Eyx1, *pml_Eyz0, *pml_Eyz1,
           *pml_Ezx0, *pml_Ezx1, *pml_Ezy0, *pml_Ezy1,
           *pml_Hxy0, *pml_Hxy1, *pml_Hxz0, *pml_Hxz1,
           *pml_Hyx0, *pml_Hyx1, *pml_Hyz0, *pml_Hyz1,
           *pml_Hzx0, *pml_Hzx1, *pml_Hzy0, *pml_Hzy1;
    double *kappa_H_x, *kappa_H_y, *kappa_H_z,
           *kappa_E_x, *kappa_E_y, *kappa_E_z,
           *bHx, *bHy, *bHz, *bEx, *bEy, *bEz,
           *cHx, *cHy, *cHz, *cEx, *cEy, *cEz;

    kernelpar *_kpar_host, *_kpar_device;

    // std::ifstream pfin;
    // pfin.open("//home/pesun//emopt//examples//Silicon_Grating_Coupler//fielddata.dat");
    // if(pfin.fail()){ std::cerr << "open failure:" << strerror(errno) << "\n";}
    // else{ std::cerr << "data file open successful\n"; }
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> t;
    // pfin >> i0 >> j0 >> k0 >> I >> J >> K >> Nx >> Ny >> Nz;
    // pfin >> pml_xmin >> pml_xmax >> pml_ymin >> pml_ymax >> pml_zmin >> pml_zmax;
    // pfin >> w_pml_x0 >> w_pml_x1 >> w_pml_y0 >> w_pml_y1 >> w_pml_z0 >> w_pml_z1;
    // pfin >> size;
    // pfin >> bc0 >> bc1 >> bc2;
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> odx;
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> dt;
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> src_t;
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> srclen;

    i0=0; j0=0; k0=0; I=87; J=572; K=585; Nx=585; Ny=572; Nz=87;
    pml_xmin=16; pml_xmax=569; pml_ymin=16; pml_ymax=556; pml_zmin=16; pml_zmax=71;
    w_pml_x0=16; w_pml_x1=16; w_pml_y0=16; w_pml_y1=16; w_pml_z0=16; w_pml_z1=16;
    size=29111940; odx=6.9497658483460967; dt=0.1139619832606030; srclen=1;

    i0s = (int *)malloc(srclen * sizeof(int));
    j0s = (int *)malloc(srclen * sizeof(int));
    k0s = (int *)malloc(srclen * sizeof(int));
    Is = (int *)malloc(srclen * sizeof(int));
    Js = (int *)malloc(srclen * sizeof(int));
    Ks = (int *)malloc(srclen * sizeof(int));
    size_t sizeall=0;
    i0s[0] = 33; j0s[0] = 18; k0s[0] = 18; Is[0]=30; Js[0]=536; Ks[0]=1;
    sizeall += Is[0]*Js[0]+Ks[0];
    src_T =604.8014773776548054; src_min=1e-4; src_k = 39714.1518172891592258;
    // for(int i=0; i< srclen; i++){
    //     pfin >> i0s[i];
    //     pfin >> j0s[i];
    //     pfin >> k0s[i];
    //     pfin >> Is[i];
    //     pfin >> Js[i];
    //     pfin >> Ks[i];
    //     sizeall += Is[i] * Js[i] * Ks[i];
    // }
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> src_T;
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> src_min;
    // pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> src_k;

    /*
        Initialize materials
    */
    epsx = (complex128 *)malloc(size * sizeof(complex128));
    epsy = (complex128 *)malloc(size * sizeof(complex128));
    epsz = (complex128 *)malloc(size * sizeof(complex128));
    complex128 matfill{};
    matfill.real=12.85;//matfill.imag=0.0;
    std::fill(epsx, epsx+size, matfill);
    std::fill(epsy, epsy+size, matfill);
    std::fill(epsz, epsz+size, matfill);

    /*
        Initialize sources
    */
    // for(int i=0; i< size; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> epsx[i].real;
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> epsx[i].imag;
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> epsy[i].real;
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> epsy[i].imag;
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> epsz[i].real;
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> epsz[i].imag;
    // }
    Jx = (complex128 *)malloc(sizeall * sizeof(complex128));
    Jy = (complex128 *)malloc(sizeall * sizeof(complex128));
    Jz = (complex128 *)malloc(sizeall * sizeof(complex128));
    Mx = (complex128 *)malloc(sizeall * sizeof(complex128));
    My = (complex128 *)malloc(sizeall * sizeof(complex128));
    Mz = (complex128 *)malloc(sizeall * sizeof(complex128));
    complex128 srcfill{};
    srcfill.real=0.7e-3;
    srcfill.imag=1.0e-3;
    std::fill(Jx, Jx+sizeall, srcfill);
    std::fill(Jy, Jy+sizeall, srcfill);
    std::fill(Jz, Jz+sizeall, srcfill);
    std::fill(Mx, Mx+sizeall, srcfill);
    std::fill(My, My+sizeall, srcfill);
    std::fill(Mz, Mz+sizeall, srcfill);

    // size_t src_offset=0;
    size_t srcsize=0;
    // for(int i=0; i<srclen; i++){
    //     srcsize = Is[i] * Js[i] * Ks[i];
    //     for(int j=0; j<srcsize; j++){
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Jx[src_offset+j].real;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Jx[src_offset+j].imag;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Jy[src_offset+j].real;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Jy[src_offset+j].imag;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Jz[src_offset+j].real;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Jz[src_offset+j].imag;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Mx[src_offset+j].real;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Mx[src_offset+j].imag;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> My[src_offset+j].real;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> My[src_offset+j].imag;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Mz[src_offset+j].real;
    //         pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> Mz[src_offset+j].imag;
    //     }
    // }

    /*
        Initialize PML on host
    */
    int Npmlx = w_pml_x0 + w_pml_x1,
        Npmly = w_pml_y0 + w_pml_y1,
        Npmlz = w_pml_z0 + w_pml_z1;
    int N=0;
    N = I * J * (w_pml_x0 - k0);
    pml_Eyx0 = (double *)malloc(N * sizeof(double));
    pml_Ezx0 = (double *)malloc(N * sizeof(double));
    pml_Hyx0 = (double *)malloc(N * sizeof(double));
    pml_Hzx0 = (double *)malloc(N * sizeof(double));
    // for(int i=0; i<N; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Eyx0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Ezx0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hyx0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hzx0[i];
    // }
    N = I * J * (k0 + K -(Nx-w_pml_x1));
    pml_Eyx1 = (double *)malloc(N * sizeof(double));
    pml_Ezx1 = (double *)malloc(N * sizeof(double));
    pml_Hyx1 = (double *)malloc(N * sizeof(double));
    pml_Hzx1 = (double *)malloc(N * sizeof(double));
    // for(int i=0; i<N; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Eyx1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Ezx1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hyx1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hzx1[i];
    // }
    N = I * K * (w_pml_y0 - j0);
    pml_Exy0 = (double *)malloc(N * sizeof(double));
    pml_Ezy0 = (double *)malloc(N * sizeof(double));
    pml_Hxy0 = (double *)malloc(N * sizeof(double));
    pml_Hzy0 = (double *)malloc(N * sizeof(double));
    // for(int i=0; i<N; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Exy0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Ezy0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hxy0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hzy0[i];
    // }
    N = I * K * (j0 + J - (Ny-w_pml_y1));
    pml_Exy1 = (double *)malloc(N * sizeof(double));
    pml_Ezy1 = (double *)malloc(N * sizeof(double));
    pml_Hxy1 = (double *)malloc(N * sizeof(double));
    pml_Hzy1 = (double *)malloc(N * sizeof(double));
    // for(int i=0; i<N; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Exy1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Ezy1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hxy1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hzy1[i];
    // }
    N = J * K * (w_pml_z0 - i0);
    pml_Exz0 = (double *)malloc(N * sizeof(double));
    pml_Eyz0 = (double *)malloc(N * sizeof(double));
    pml_Hxz0 = (double *)malloc(N * sizeof(double));
    pml_Hyz0 = (double *)malloc(N * sizeof(double));
    // for(int i=0; i<N; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Exz0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Eyz0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hxz0[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hyz0[i];
    // }
    N = J * K * (i0 + I -(Nz-w_pml_z1));
    pml_Exz1 = (double *)malloc(N * sizeof(double));
    pml_Eyz1 = (double *)malloc(N * sizeof(double));
    pml_Hxz1 = (double *)malloc(N * sizeof(double));
    pml_Hyz1 = (double *)malloc(N * sizeof(double));
    // for(int i=0; i<N; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Exz1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Eyz1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hxz1[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> pml_Hyz1[i];
    // }

    Npmlx = w_pml_x0 + w_pml_x1;
    kappa_H_x = (double *)malloc(Npmlx * sizeof(double));
    kappa_E_x = (double *)malloc(Npmlx * sizeof(double));
    bHx = (double *)malloc(Npmlx * sizeof(double));
    bEx = (double *)malloc(Npmlx * sizeof(double));
    cHx = (double *)malloc(Npmlx * sizeof(double));
    cEx = (double *)malloc(Npmlx * sizeof(double));
    // for(int i=0; i<Npmlx; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> kappa_H_x[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> kappa_E_x[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> bHx[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> bEx[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> cHx[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> cEx[i];
    // }
    Npmly = w_pml_y0 + w_pml_y1;
    kappa_H_y = (double *)malloc(Npmly * sizeof(double));
    kappa_E_y = (double *)malloc(Npmly * sizeof(double));
    bHy = (double *)malloc(Npmly * sizeof(double));
    bEy = (double *)malloc(Npmly * sizeof(double));
    cHy = (double *)malloc(Npmly * sizeof(double));
    cEy = (double *)malloc(Npmly * sizeof(double));
    // for(int i=0; i<Npmly; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> kappa_H_y[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> kappa_E_y[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> bHy[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> bEy[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> cHy[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> cEy[i];
    // }
    Npmlz = w_pml_z0 + w_pml_z1;
    kappa_H_z = (double *)malloc(Npmlz * sizeof(double));
    kappa_E_z = (double *)malloc(Npmlz * sizeof(double));
    bHz = (double *)malloc(Npmlz * sizeof(double));
    bEz = (double *)malloc(Npmlz * sizeof(double));
    cHz = (double *)malloc(Npmlz * sizeof(double));
    cEz = (double *)malloc(Npmlz * sizeof(double));
    // for(int i=0; i<Npmlz; i++){
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> kappa_H_z[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> kappa_E_z[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> bHz[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> bEz[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> cHz[i];
    //     pfin >> std::fixed >> std::setprecision(std::numeric_limits<double>::digits10+1) >> cEz[i];
    // }
    double pml_dist, pml_factor, sigma, alpha, kappa, b, c;
    double pow;
    double sigma0, alpha0, kappa0;
    sigma0=3.0; alpha0=0.0; kappa0=1.0; pow=3.0;
    for(int k = 0; k < w_pml_x0; k++) {
        pml_dist = double(k - 0.5)/w_pml_x0; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);
        if(pml_factor < 0) pml_factor = 0;

        // compute H coefficients
        sigma = sigma0 * pml_factor;
        alpha = alpha0 * (1-pml_factor);
        kappa = (kappa0-1.0) * pml_factor+1.0;
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_H_x[k] = kappa;
        bHx[k] = b;
        cHx[k] = c;

        pml_dist = double(k)/w_pml_x0; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        // compute E coefficients
        sigma = sigma0 * pml_factor;
        alpha = alpha0 * (1-pml_factor);
        kappa = (kappa0-1) * pml_factor+1;
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_E_x[k] = kappa;
        bEx[k] = b;
        cEx[k] = c;

    }
    for(int k = 0; k < w_pml_x1; k++) {
        // compute H coefficients
        pml_dist = double(k + 0.5)/w_pml_x1; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_H_x[w_pml_x0 + k] = kappa;
        bHx[w_pml_x0 + k] = b;
        cHx[w_pml_x0 + k] = c;

        //compute E coefficients
        pml_dist = double(k)/w_pml_x1; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_E_x[w_pml_x0 + k] = kappa;
        bEx[w_pml_x0 + k] = b;
        cEx[w_pml_x0 + k] = c;
    }
    for(int j = 0; j < w_pml_y0; j++) {
        // calc H coefficients
        pml_dist = double(j - 0.5)/w_pml_y0; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);
        if(pml_factor < 0) pml_factor = 0;

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_H_y[j] = kappa;
        bHy[j] = b;
        cHy[j] = c;

        // calc E coefficients
        pml_dist = double(j)/w_pml_y0; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_E_y[j] = kappa;
        bEy[j] = b;
        cEy[j] = c;
    
    }
    for(int j = 0; j < w_pml_y1; j++) {
        // calc H coeffs
        pml_dist = double(j + 0.5)/w_pml_y1; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_H_y[w_pml_y0 + j] = kappa;
        bHy[w_pml_y0 + j] = b;
        cHy[w_pml_y0 + j] = c;

        // compute E coefficients
        pml_dist = double(j)/w_pml_y1; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha); 

        kappa_E_y[w_pml_y0 + j] = kappa;
        bEy[w_pml_y0 + j] = b;
        cEy[w_pml_y0 + j] = c;
    }
    for(int i = 0; i < w_pml_z0; i++) {
        // calc H coeffs
        pml_dist = double(i)/w_pml_z0; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c= 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_H_z[i] = kappa;
        bHz[i] = b;
        cHz[i] = c;

        // calc E coeffs
        pml_dist = double(i+0.5)/w_pml_z0; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        // compute coefficients
        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_E_z[i] = kappa;
        bEz[i] = b;
        cEz[i] = c;
    }
    for(int i = 0; i < w_pml_z1; i++) {
        // calc H coeffs
        pml_dist = double(i)/w_pml_z1; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);

        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_H_z[w_pml_z0 + i] = kappa;
        bHz[w_pml_z0 + i] = b;
        cHz[w_pml_z0 + i] = c;

        // calc E coeffs
        pml_dist = double(i - 0.5)/w_pml_z1; // distance from pml edge
        pml_factor = std::pow(pml_dist,pow);
        if(pml_factor < 0) pml_factor = 0;

        // compute coefficients
        sigma = sigma0 * pml_factor;
        kappa = (kappa0-1) * pml_factor+1;
        alpha = alpha0 * (1-pml_factor);
        b = exp(-dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        kappa_E_z[w_pml_z0 + i] = kappa;
        bEz[w_pml_z0 + i] = b;
        cEz[w_pml_z0 + i] = c;
    }

    // pfin.close();

    /*
        Initialize device arrays
    */
    complex128 *dJx, *dJy, *dJz, *dMx, *dMy, *dMz;
    int *di0s, *dj0s, *dk0s, *dIs, *dJs, *dKs;
    size_t size_offset = 0;
    char *dbc;
    double *dEx, *dEy, *dEz, *dHx, *dHy, *dHz;
    complex128 *depsx, *depsy, *depsz;
    double *dpml_Exy0, *dpml_Exy1, *dpml_Exz0, *dpml_Exz1,
        *dpml_Eyx0, *dpml_Eyx1, *dpml_Eyz0, *dpml_Eyz1,
        *dpml_Ezx0, *dpml_Ezx1, *dpml_Ezy0, *dpml_Ezy1,
        *dpml_Hxy0, *dpml_Hxy1, *dpml_Hxz0, *dpml_Hxz1,
        *dpml_Hyx0, *dpml_Hyx1, *dpml_Hyz0, *dpml_Hyz1,
        *dpml_Hzx0, *dpml_Hzx1, *dpml_Hzy0, *dpml_Hzy1;
    double *dkappa_H_x, *dkappa_H_y, *dkappa_H_z,
        *dkappa_E_x, *dkappa_E_y, *dkappa_E_z,
        *dbHx, *dbHy, *dbHz,
        *dbEx, *dbEy, *dbEz,
        *dcHx, *dcHy, *dcHz,
        *dcEx, *dcEy, *dcEz;

    double *dE_bufferleft, *dE_bufferright, *dH_bufferleft, *dH_bufferright,
        *dE_ghostleft, *dE_ghostright, *dH_ghostleft, *dH_ghostright;

    int gpus_count{}, Ngpus{2};
    int FDTD_steps{3000};
    char domain_decomp{'x'};
    cudaGetDeviceCount (&Ngpus);

    /*
        Allocate and copy data to GPUs
    */
    int It{}, Jt{}, Kt{}, i0t{}, j0t{}, k0t{};

    //     /*
    //         Domain decomposition: along X for development
    //     */
    gpus_count=Ngpus;
    if(argc>1)
        gpus_count = std::stoi(argv[1]);
    if(argc>2){
        gpus_count = std::stoi(argv[1]);
        FDTD_steps = std::stoi(argv[2]);
    }
    if(argc>3){
        gpus_count = std::stoi(argv[1]);
        FDTD_steps = std::stoi(argv[2]);
        domain_decomp = argv[3][0];
    }
    std::cout << "Use " << gpus_count << " of " << Ngpus << " available GPUs" << std::endl;
    std::cout << "Total No. of FDTD steps: " << FDTD_steps << std::endl;

    std::vector<kernelpar *> kpar_list;
    std::vector<dim3> griddim_list;
    std::vector<dim3> griddimghost_list;

    kpar_list.reserve(gpus_count);
    kernelpar** kpar_list_host = new kernelpar*[gpus_count];

    griddim_list.reserve(gpus_count);
    griddimghost_list.reserve(gpus_count);

    size_t ghost_size{};
    size_t buffer_size{};
    std::vector<size_t> ghostsize_list;
    std::vector<size_t> buffersize_list;
    ghostsize_list.reserve(gpus_count);
    buffersize_list.reserve(gpus_count);

    for(int device_id=0; device_id<gpus_count; device_id++)
    {
        cudaSetDevice(device_id);
        _kpar_host = (kernelpar *)malloc(sizeof(kernelpar));
        gpuErrchk(cudaMalloc((void **)&_kpar_device, sizeof(kernelpar)));
        kpar_list.push_back(_kpar_device);

        /*
            Domain decomposition: along X for development
        */
        // i0t = i0; j0t = j0; k0t = device_id*Nx/gpus_count;
        // It = I; Jt = J; Kt = (device_id==gpus_count-1)?Nx-device_id*Nx/gpus_count:(device_id+1)*Nx/gpus_count-device_id*Nx/gpus_count;
        switch(domain_decomp){
            case 'x':
                i0t = i0; j0t = j0; k0t = device_id*Nx/gpus_count;
                It = I; Jt = J; Kt = (device_id==gpus_count-1)?Nx-device_id*Nx/gpus_count:(device_id+1)*Nx/gpus_count-device_id*Nx/gpus_count;
                ghost_size = I*J;
                buffer_size = (k0t==0 || k0t+Kt==Nx) ? I*J : 2*I*J;
                break;
            case 'y':
                i0t = i0; j0t = device_id*Ny/gpus_count; k0t = k0;
                It = I; Jt = (device_id==gpus_count-1)?Ny-device_id*Ny/gpus_count:(device_id+1)*Ny/gpus_count-device_id*Ny/gpus_count; Kt = K;
                ghost_size = I*K;
                buffer_size = (j0t==0 || j0t+Jt==Ny) ? I*K : 2*I*K;
                break;
            case 'z':
                i0t = device_id*Nz/gpus_count; j0t = j0; k0t = k0;
                It = (device_id==gpus_count-1)?Nz-device_id*Nz/gpus_count:(device_id+1)*Nz/gpus_count-device_id*Nz/gpus_count; Jt = J; Kt = K;
                ghost_size = J*K;
                buffer_size = (i0t==0 || i0t+It==Nz) ? J*K : 2*J*K;
                break;        
        }
        ghostsize_list.push_back(ghost_size);
        buffersize_list.push_back(buffer_size);

        _kpar_host->I = It; _kpar_host->J = Jt; _kpar_host->K = Kt;
        _kpar_host->i0 = i0t; _kpar_host->j0 = j0t; _kpar_host->k0 = k0t;
        // _kpar_host->I = I; _kpar_host->J = J; _kpar_host->K = K;
        // _kpar_host->i0 = i0; _kpar_host->j0 = j0; _kpar_host->k0 = k0;
 
        _kpar_host->size = _kpar_host->I * _kpar_host->J * _kpar_host->K; //It*Jt*Kt;
        _kpar_host->ghost_size = ghost_size;
        _kpar_host->buffer_size = buffer_size;

        std::cout << "Thread " << device_id << ", (I, i0)=(" << _kpar_host->I << ", " << _kpar_host->i0 << "), (J,j0)=" 
                << _kpar_host->J << ", " << _kpar_host->j0 << "), (K, k0)=(" << _kpar_host->K << ", " << _kpar_host->k0 << ")" << std::endl;

        griddim_list.push_back(ceil(_kpar_host->I*_kpar_host->J*_kpar_host->K/(double)(blockdim)));
        griddimghost_list.push_back(ceil(buffer_size/(double)(blockdim)));

        _kpar_host->Nx = Nx; _kpar_host->Ny = Ny; _kpar_host->Nz = Nz;
        _kpar_host->dt = dt; _kpar_host->odx = odx; _kpar_host->t = t;
        _kpar_host->src_T = src_T; _kpar_host->src_min = src_min; _kpar_host->src_k = src_k;

        /*
            Fields
        */
       // size in X is always padded to Kt+2: ghost layers for inner edges and boundary condition for outer edges
        gpuErrchk(cudaMalloc((void **)&dEx, (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2) * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dEy, (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2) * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dEz, (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2) * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dHx, (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2) * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dHy, (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2) * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dHz, (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2) * sizeof(double)));

        gpuErrchk(cudaMalloc((void **)&dE_bufferleft, 3*ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_bufferleft, 3*ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dE_bufferright, 3*ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_bufferright, 3*ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dE_ghostleft, 3*ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_ghostleft, 3*ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dE_ghostright, 3*ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_ghostright, 3*ghost_size*sizeof(double)));

        _kpar_host->E_bufferleft = dE_bufferleft;
        _kpar_host->H_bufferleft = dH_bufferleft;
        _kpar_host->E_bufferright = dE_bufferright;
        _kpar_host->H_bufferright = dH_bufferright;
        _kpar_host->E_ghostleft = dE_ghostleft;
        _kpar_host->H_ghostleft = dH_ghostleft;
        _kpar_host->E_ghostright = dE_ghostright;
        _kpar_host->H_ghostright = dH_ghostright;

        size_t Exlen = (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2);
        size_t blockSize = 256;
        size_t numBlocks = (Exlen+blockSize-1)/Exlen;
        initMat<<<numBlocks, blockSize>>>(dEx, Exlen, 0.0);
        initMat<<<numBlocks, blockSize>>>(dEy, Exlen, 0.0);
        initMat<<<numBlocks, blockSize>>>(dEz, Exlen, 0.0);
        initMat<<<numBlocks, blockSize>>>(dHx, Exlen, 0.0);
        initMat<<<numBlocks, blockSize>>>(dHy, Exlen, 0.0);
        initMat<<<numBlocks, blockSize>>>(dHz, Exlen, 0.0);

        // Exlen = (_kpar_host->I)*(_kpar_host->J);
        // numBlocks = (Exlen+blockSize-1)/Exlen;
        // initMat<<<numBlocks, blockSize>>>(dEx_ghostleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEy_ghostleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEz_ghostleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHx_ghostleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHy_ghostleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHz_ghostleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEx_bufferleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEy_bufferleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEz_bufferleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHx_bufferleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHy_bufferleft, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHz_bufferleft, Exlen, 0.0);

        // initMat<<<numBlocks, blockSize>>>(dEx_ghostright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEy_ghostright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEz_ghostright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHx_ghostright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHy_ghostright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHz_ghostright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEx_bufferright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEy_bufferright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dEz_bufferright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHx_bufferright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHy_bufferright, Exlen, 0.0);
        // initMat<<<numBlocks, blockSize>>>(dHz_bufferright, Exlen, 0.0);

        _kpar_host->Ex = dEx; _kpar_host->Ey = dEy; _kpar_host->Ez = dEz;
        _kpar_host->Hx = dHx; _kpar_host->Hy = dHy; _kpar_host->Hz = dHz;

        /*
            Materials
        */
        gpuErrchk(cudaMalloc((void **)&depsx, _kpar_host->I*_kpar_host->J*_kpar_host->K * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&depsy, _kpar_host->I*_kpar_host->J*_kpar_host->K * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&depsz, _kpar_host->I*_kpar_host->J*_kpar_host->K * sizeof(complex128)));
        gpuErrchk(cudaMemcpy(depsx, epsx, _kpar_host->I*_kpar_host->J*_kpar_host->K * sizeof(complex128), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(depsy, epsy, _kpar_host->I*_kpar_host->J*_kpar_host->K * sizeof(complex128), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(depsz, epsz, _kpar_host->I*_kpar_host->J*_kpar_host->K * sizeof(complex128), cudaMemcpyHostToDevice));

        _kpar_host->epsx = depsx; _kpar_host->epsy = depsy; _kpar_host->epsz = depsz;

        /*
            PML
        */
        if(k0t< w_pml_x0){
            N = _kpar_host->I*_kpar_host->J * (w_pml_x0 - _kpar_host->k0);
            gpuErrchk(cudaMalloc((void **)&dpml_Eyx0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezx0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyx0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzx0, N * sizeof(double)));
        }
        if(k0t+Kt>Nx-w_pml_x1){
            N = _kpar_host->I*_kpar_host->J * (_kpar_host->k0 + _kpar_host->K -(Nx-w_pml_x1));
            gpuErrchk(cudaMalloc((void **)&dpml_Eyx1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezx1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyx1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzx1, N * sizeof(double)));
        }
        if(j0t<w_pml_y0){
            N = _kpar_host->I*_kpar_host->K * (w_pml_y0 - _kpar_host->j0);
            gpuErrchk(cudaMalloc((void **)&dpml_Exy0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezy0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxy0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzy0, N * sizeof(double)));
        }
        if(j0t+Jt>Ny-w_pml_y1){
            N = _kpar_host->I*_kpar_host->K * (_kpar_host->j0 + _kpar_host->J - (Ny-w_pml_y1));
            gpuErrchk(cudaMalloc((void **)&dpml_Exy1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezy1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxy1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzy1, N * sizeof(double)));
        }
        if(i0t<w_pml_z0){
            N = _kpar_host->J*_kpar_host->K * (w_pml_z0 - _kpar_host->i0);
            gpuErrchk(cudaMalloc((void **)&dpml_Exz0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Eyz0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxz0, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyz0, N * sizeof(double)));
        }
        if(i0t+It>Nz-w_pml_z1){
            N = _kpar_host->J*_kpar_host->K * (_kpar_host->i0 + _kpar_host->I -(Nz-w_pml_z1));
            gpuErrchk(cudaMalloc((void **)&dpml_Exz1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Eyz1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxz1, N * sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyz1, N * sizeof(double)));
        }

        gpuErrchk(cudaMalloc((void **)&dkappa_H_x, Npmlx * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_H_y, Npmly * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_H_z, Npmlz * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_E_x, Npmlx * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_E_y, Npmly * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_E_z, Npmlz * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbHx, Npmlx * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbHy, Npmly * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbHz, Npmlz * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbEx, Npmlx * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbEy, Npmly * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbEz, Npmlz * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcHx, Npmlx * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcHy, Npmly * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcHz, Npmlz * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcEx, Npmlx * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcEy, Npmly * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcEz, Npmlz * sizeof(double)));

        gpuErrchk(cudaMemcpy(dkappa_H_x, kappa_H_x, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_H_y, kappa_H_y, Npmly * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_H_z, kappa_H_z, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_E_x, kappa_E_x, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_E_y, kappa_E_y, Npmly * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_E_z, kappa_E_z, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbHx, bHx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbHy, bHy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbHz, bHz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbEx, bEx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbEy, bEy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbEz, bEz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcHx, cHx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcHy, cHy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcHz, cHz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcEx, cEx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcEy, cEy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcEz, cEz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));

        _kpar_host->pml_Exy0 = dpml_Exy0;
        _kpar_host->pml_Exy1 = dpml_Exy1;
        _kpar_host->pml_Exz0 = dpml_Exz0;
        _kpar_host->pml_Exz1 = dpml_Exz1;
        _kpar_host->pml_Eyx0 = dpml_Eyx0;
        _kpar_host->pml_Eyx1 = dpml_Eyx1;
        _kpar_host->pml_Eyz0 = dpml_Eyz0;
        _kpar_host->pml_Eyz1 = dpml_Eyz1;
        _kpar_host->pml_Ezx0 = dpml_Ezx0;
        _kpar_host->pml_Ezx1 = dpml_Ezx1;
        _kpar_host->pml_Ezy0 = dpml_Ezy0;
        _kpar_host->pml_Ezy1 = dpml_Ezy1;

        _kpar_host->pml_Hxy0 = dpml_Hxy0;
        _kpar_host->pml_Hxy1 = dpml_Hxy1;
        _kpar_host->pml_Hxz0 = dpml_Hxz0;
        _kpar_host->pml_Hxz1 = dpml_Hxz1;
        _kpar_host->pml_Hyx0 = dpml_Hyx0;
        _kpar_host->pml_Hyx1 = dpml_Hyx1;
        _kpar_host->pml_Hyz0 = dpml_Hyz0;
        _kpar_host->pml_Hyz1 = dpml_Hyz1;
        _kpar_host->pml_Hzx0 = dpml_Hzx0;
        _kpar_host->pml_Hzx1 = dpml_Hzx1;
        _kpar_host->pml_Hzy0 = dpml_Hzy0;
        _kpar_host->pml_Hzy1 = dpml_Hzy1;

        _kpar_host->kappa_H_x = dkappa_H_x;
        _kpar_host->kappa_H_y = dkappa_H_y;
        _kpar_host->kappa_H_z = dkappa_H_z;
        _kpar_host->kappa_E_x = dkappa_E_x;
        _kpar_host->kappa_E_y = dkappa_E_y;
        _kpar_host->kappa_E_z = dkappa_E_z;
        _kpar_host->bHx = dbHx;
        _kpar_host->bHy = dbHy;
        _kpar_host->bHz = dbHz;
        _kpar_host->bEx = dbEx;
        _kpar_host->bEy = dbEy;
        _kpar_host->bEz = dbEz;
        _kpar_host->cHx = dcHx;
        _kpar_host->cHy = dcHy;
        _kpar_host->cHz = dcHz;
        _kpar_host->cEx = dcEx;
        _kpar_host->cEy = dcEy;
        _kpar_host->cEz = dcEz;

        _kpar_host->w_pml_x0 = w_pml_x0;
        _kpar_host->w_pml_x1 = w_pml_x1;
        _kpar_host->w_pml_y0 = w_pml_y0;
        _kpar_host->w_pml_y1 = w_pml_y1;
        _kpar_host->w_pml_z0 = w_pml_z0;
        _kpar_host->w_pml_z1 = w_pml_z1;
        _kpar_host->pml_xmin = pml_xmin;
        _kpar_host->pml_xmax = pml_xmax;
        _kpar_host->pml_ymin = pml_ymin;
        _kpar_host->pml_ymax = pml_ymax;
        _kpar_host->pml_zmin = pml_zmin;
        _kpar_host->pml_zmax = pml_zmax;

        /*
            Sources
        */
        _kpar_host->srclen = srclen;
        gpuErrchk(cudaMalloc((void **)&dJx, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dJy, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dJz, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dMx, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dMy, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dMz, sizeall * sizeof(complex128)));

        gpuErrchk(cudaMalloc((void **)&di0s, srclen * sizeof(int)));
        gpuErrchk(cudaMalloc((void **)&dj0s, srclen * sizeof(int)));
        gpuErrchk(cudaMalloc((void **)&dk0s, srclen * sizeof(int)));
        gpuErrchk(cudaMalloc((void **)&dIs, srclen * sizeof(int)));
        gpuErrchk(cudaMalloc((void **)&dJs, srclen * sizeof(int)));
        gpuErrchk(cudaMalloc((void **)&dKs, srclen * sizeof(int)));

        gpuErrchk(cudaMemcpy(di0s, i0s, srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dj0s, j0s, srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dk0s, k0s, srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dIs, Is, srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dJs, Js, srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dKs, Ks, srclen * sizeof(int), cudaMemcpyHostToDevice));

        size_offset=0;
        for(int i=0; i < srclen; i++){
            srcsize = Is[i] * Js[i] * Ks[i];
            gpuErrchk(cudaMemcpy(dJx+size_offset, Jx+size_offset, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dJy+size_offset, Jy+size_offset, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dJz+size_offset, Jz+size_offset, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dMx+size_offset, Mx+size_offset, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dMy+size_offset, My+size_offset, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dMz+size_offset, Mz+size_offset, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            size_offset += srcsize;
        }

        _kpar_host->i0s = di0s;
        _kpar_host->j0s = dj0s;
        _kpar_host->k0s = dk0s;
        _kpar_host->Is = dIs;
        _kpar_host->Js = dJs;
        _kpar_host->Ks = dKs;
        _kpar_host->Jx = dJx;
        _kpar_host->Jy = dJy;
        _kpar_host->Jz = dJz;
        _kpar_host->Mx = dMx;
        _kpar_host->My = dMy;
        _kpar_host->Mz = dMz;

        // std::cout << "Thread " << device_id << ", (Is, i0s)=(" << Is[0] << ", " << i0s[0] << "), (Js,j0s)=" 
        //         << Js[0] << ", " << j0s[0] << "), (K, k0)=(" << Ks[0] << ", " << k0s[0] << ",)" << std::endl;

        /*
            Boundary conditions
        */
        gpuErrchk(cudaMalloc((void **)&dbc, 3 * sizeof(char)));
        gpuErrchk(cudaMemcpy(dbc, bc, 3 * sizeof(char), cudaMemcpyHostToDevice));
        _kpar_host->bc = dbc;

        gpuErrchk(cudaMemcpy(kpar_list[device_id], _kpar_host, sizeof(kernelpar), cudaMemcpyHostToDevice));
        kpar_list_host[device_id]= _kpar_host;
    }

    std::vector<std::thread> threads_nvlink, threads_kernels;
    threads_nvlink.reserve(gpus_count);

    /*
        Enable NVlink
    */
    bool Enable_NVLINK{true};
    // if(Enable_NVLINK==true)
    // {
    //     for(int device_id=0; device_id<gpus_count; device_id++)
    //     {
    //         thread_info_class thread_info (device_id, gpus_count, barrier);
    //         threads_nvlink.emplace_back([thread_info] () {
    //             try{
    //                 prepare_nvlink (thread_info);
    //             }
    //             catch (std::runtime_error &error) {
    //                 std::cerr << "Error in thread " << thread_info.thread_id << ":" << error.what() <<std::endl;
    //             }
    //         });
    //     }
    //     for (auto &thread: threads_nvlink)
    //         thread.join();
    // }

    if(Enable_NVLINK==true)
    {
        for(int device_id=0; device_id<gpus_count; device_id++)
        {
            threads_nvlink.emplace_back([=](){
                try{ prepare_nvlink (device_id, gpus_count); }
                catch(std::runtime_error &error){ std::cerr << "Error in thread " << device_id << ":" << error.what() <<std::endl; }
            });
        }
        for (auto &thread: threads_nvlink)
            thread.join();
    }

    /*
        Prepare streams
    */
    std::vector<cudaStream_t> compute_stream_list, sync_stream_list;
    std::vector<cudaEvent_t> e_border_event_list, h_border_event_list, e_bulk_event_list, h_bulk_event_list, e_time_event_list, h_time_event_list;
    compute_stream_list.reserve(gpus_count);
    sync_stream_list.reserve(gpus_count);
    e_border_event_list.reserve(gpus_count);
    h_border_event_list.reserve(gpus_count);
    e_bulk_event_list.reserve(gpus_count);
    h_bulk_event_list.reserve(gpus_count);
    e_time_event_list.reserve(gpus_count);
    h_time_event_list.reserve(gpus_count);
    int least_priority{}, greatest_priority{};
    cudaDeviceGetStreamPriorityRange(&least_priority, &greatest_priority);

    for(int device_id=0; device_id<gpus_count; device_id++)
    {
        cudaSetDevice(device_id%gpus_count);
        cudaStream_t compute_stream, sync_stream;
        gpuErrchk(cudaStreamCreateWithPriority(&compute_stream, cudaStreamDefault, least_priority));
        compute_stream_list.push_back(compute_stream);
        gpuErrchk(cudaStreamCreateWithPriority(&sync_stream, cudaStreamDefault, greatest_priority));
        sync_stream_list.push_back(sync_stream);

        cudaEvent_t h_bulk_computed, e_bulk_computed, h_border_computed, e_border_computed, h_time_computed, e_time_computed;
        cudaEventCreateWithFlags(&h_bulk_computed, cudaEventDisableTiming);
        cudaEventCreateWithFlags(&e_bulk_computed, cudaEventDisableTiming);
        cudaEventCreateWithFlags(&h_border_computed, cudaEventDisableTiming);
        cudaEventCreateWithFlags(&e_border_computed, cudaEventDisableTiming);
        cudaEventCreateWithFlags(&h_time_computed, cudaEventDisableTiming);
        cudaEventCreateWithFlags(&e_time_computed, cudaEventDisableTiming);
        h_bulk_event_list.push_back(h_bulk_computed);
        e_bulk_event_list.push_back(e_bulk_computed);
        h_border_event_list.push_back(h_border_computed);
        e_border_event_list.push_back(e_border_computed);
        h_time_event_list.push_back(h_time_computed);
        e_time_event_list.push_back(e_time_computed);
    }

    /*
        Kernel launch
    */
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    threads_kernels.reserve(gpus_count);
    
    for(int device_id=0; device_id<gpus_count; device_id++)
    {
        threads_kernels.emplace_back([=] () {
            try{
                cudaSetDevice(device_id%gpus_count);
                for(int step=0;step<FDTD_steps;step++)
                {
                    kpar_list_host[device_id]->t = step*dt;

                    gpuErrchk(cudaMemcpyAsync(kpar_list[device_id], kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(h_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], h_time_event_list[device_id], 0);
                    cudaStreamWaitEvent(compute_stream_list[device_id], e_border_event_list[device_id], 0);
                    kernel_update_H_bulk <<< griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (kpar_list[device_id]);
                    cudaEventRecord(h_bulk_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], h_time_event_list[device_id], 0);
                    cudaStreamWaitEvent(sync_stream_list[device_id], e_bulk_event_list[device_id], 0);
                    cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id], 0);
                    kernel_update_H_border <<< griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (kpar_list[device_id]);

                    cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);

                    if(gpus_count>1){
                        Barrier(counter_H, mutex_H, cv_H, gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->H_bufferright, device_id, kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            kernel_load_ghostlayerH <<< griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (kpar_list[device_id]);
                        }
                        else if(device_id==gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->H_bufferleft, device_id, kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            kernel_load_ghostlayerH <<< griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (kpar_list[device_id]);
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->H_bufferright, device_id, kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->H_bufferleft, device_id, kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            kernel_load_ghostlayerH <<< griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (kpar_list[device_id]);
                        }
                        kernel_load_ghostlayerH <<< griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (kpar_list[device_id]);
                        cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);
                    }

                    kpar_list_host[device_id]->t = (step+0.5)*dt;
                    gpuErrchk(cudaMemcpyAsync(kpar_list[device_id], kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(e_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], e_time_event_list[device_id], 0);
                    cudaStreamWaitEvent(compute_stream_list[device_id], h_border_event_list[device_id], 0);
                    kernel_update_E_bulk <<< griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (kpar_list[device_id]);
                    cudaEventRecord(e_bulk_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], e_time_event_list[device_id], 0);
                    cudaStreamWaitEvent(sync_stream_list[device_id], h_bulk_event_list[device_id], 0);
                    cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id], 0);
                    kernel_update_E_border <<< griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (kpar_list[device_id]);
                    cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    
                    if(gpus_count>1){
                        Barrier(counter_E, mutex_E, cv_E, gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->E_bufferright, device_id, kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                        }
                        else if(device_id==gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->E_bufferleft, device_id, kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->E_bufferright, device_id, kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            gpuErrchk(cudaMemcpyPeerAsync(kpar_list_host[device_id]->E_bufferleft, device_id, kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                        }
                        kernel_load_ghostlayerE <<< griddimghost_list[device_id], blockdim,0, sync_stream_list[device_id] >>> (kpar_list[device_id]);
                        cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    }

                }
            }
            catch (std::runtime_error &error){
                std::cerr << "FDTD error in thread " << device_id << ":" << error.what() << std::endl;
                // std::cerr << "FDTD error in thread " << thread_info.thread_id << ":" << error.what() << std::endl;
            }
        });
    }
    for (auto &thread: threads_kernels)
        thread.join();

    // cudaError_t cudaerr;
    // cudaerr = cudaDeviceSynchronize();
    // if (cudaerr != cudaSuccess)
    //     printf("FDTD status \"%s\".\n", cudaGetErrorString(cudaerr));

    // cudaError_t cudaerr{};
    // for (int step=0; step<1000; step++)
    // {
    //     kernel_update_H <<< grid_dim, block_dim >>> (_kpar_device);
        // cudaerr = cudaDeviceSynchronize();
        // if (cudaerr != cudaSuccess)
        // {
        //    printf("update H: kernel launch status \"%s\".\n", cudaGetErrorString(cudaerr));
        //    break;
        // }             
        // kernel_update_E <<< grid_dim, block_dim >>> (_kpar_device);
        // cudaerr = cudaDeviceSynchronize();
        // if (cudaerr != cudaSuccess)
        // {
        //    printf("update E: kernel launch status \"%s\".\n", cudaGetErrorString(cudaerr));
        //    break;
        // }
    // }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()/1e6 << "s" << std::endl;

    /*
        Assemble Ex
    */
    if(false){   
        double max{1.0f}, min{-1.0f};
        double *Ex, *Ext;
        Ex = (double *)malloc(I*J*K*sizeof(double));

        // int idmax =0;
        for(int device_id=0; device_id<gpus_count; device_id++)
        {
            cudaSetDevice(device_id);
            gpuErrchk(cudaMemcpy(_kpar_host, kpar_list[device_id], sizeof(kernelpar), cudaMemcpyDeviceToHost));

            // i0t = i0; j0t = j0; k0t = device_id*Nx/gpus_count;
            // It = I; Jt = J; Kt = (device_id==gpus_count-1)?Nx-device_id*Nx/gpus_count:(device_id+1)*Nx/gpus_count-device_id*Nx/gpus_count;
            switch(domain_decomp){
                case 'x':
                    i0t = i0; j0t = j0; k0t = device_id*Nx/gpus_count;
                    It = I; Jt = J; Kt = (device_id==gpus_count-1)?Nx-device_id*Nx/gpus_count:(device_id+1)*Nx/gpus_count-device_id*Nx/gpus_count;
                    break;
                case 'y':
                    i0t = i0; j0t = device_id*Ny/gpus_count; k0t = k0;
                    It = I; Jt = (device_id==gpus_count-1)?Ny-device_id*Ny/gpus_count:(device_id+1)*Ny/gpus_count-device_id*Ny/gpus_count; Kt = K;
                    break;
                case 'z':
                    i0t = device_id*Nz/gpus_count; j0t = j0; k0t = k0;
                    It = (device_id==gpus_count-1)?Nz-device_id*Nz/gpus_count:(device_id+1)*Nz/gpus_count-device_id*Nz/gpus_count; Jt = J; Kt = K;
                    break;        
            }

            Ext = (double *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(double));
            gpuErrchk(cudaMemcpy(Ext, _kpar_host->Ex, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
            gpuErrchk(cudaDeviceSynchronize());
            for(int ii=0; ii<It; ii++){
                for(int jj=0; jj<Jt; jj++){
                    for(int kk=0; kk<Kt; kk++){
                        Ex[(ii+i0t)*J*K + (jj+j0t)*K + kk+k0t] = Ext[(ii+1)*(Jt+2)*(Kt+2) + (jj+1)*(Kt+2) + kk+1];
                    }
                }
            }
            delete[] Ext;
        }

        // double max{1.0f}, min{-1.0f};
        // double* Ex;
        // // thrust::device_ptr<double> d_ptr = thrust::device_pointer_cast(_kpar_host->Ex);
        // // max = *(thrust::max_element(d_ptr, d_ptr+((I+2)*(J+2)*(K*2))));
        // // min = *(thrust::min_element(d_ptr, d_ptr+((I+2)*(J+2)*(K*2))));
        max = *std::max_element(Ex, Ex+I*J*K);
        min = *std::min_element(Ex, Ex+I*J*K);
        double sum = 0;
        sum= std::accumulate(Ex, Ex+I*J*K, sum);
        std::cout << "Max Ex:" << std::setprecision(15) << max << std::endl;
        std::cout << "Min Ex:" << std::setprecision(15) << min << std::endl;
        std::cout << "Sum Ex:" << std::setprecision(15) << sum << std::endl;

        std::ofstream out("field_slice.csv");
        out << "Ex, i, k" << '\n';
        for(int ii=55; ii<55+1; ii++){
            for(int jj=0; jj<J; jj++){
                for(int kk=0; kk<K; kk++){
                    out << Ex[ii*J*K+jj*K+kk] << ',' << jj << ',' << kk << '\n';
                }
            }
        }

        delete[] Ex;
    }

    /*
        Clean up
    */
    delete[] epsx; delete[] epsy; delete[] epsz;
    delete[] Jx; delete[] Jy; delete[] Jz; delete[] Mx; delete[] My; delete[] Mz;

    delete[] pml_Exy0; delete[] pml_Exy1; delete[] pml_Exz0; delete[] pml_Exz1;
    delete[] pml_Eyx0; delete[] pml_Eyx1; delete[] pml_Eyz0; delete[] pml_Eyz1;
    delete[] pml_Ezx0; delete[] pml_Ezx1; delete[] pml_Ezy0; delete[] pml_Ezy1;
    delete[] pml_Hxy0; delete[] pml_Hxy1; delete[] pml_Hxz0; delete[] pml_Hxz1;
    delete[] pml_Hyx0; delete[] pml_Hyx1; delete[] pml_Hyz0; delete[] pml_Hyz1;
    delete[] pml_Hzx0; delete[] pml_Hzx1; delete[] pml_Hzy0; delete[] pml_Hzy1;

    delete[] kappa_H_x; delete[] kappa_H_y; delete[] kappa_H_z; delete[] kappa_E_x; delete[] kappa_E_y; delete[] kappa_E_z;

    delete[] bHx; delete[] bHy; delete[] bHz; delete[] bEx; delete[] bEy; delete[] bEz;
    delete[] cHx; delete[] cHy; delete[] cHz; delete[] cEx; delete[] cEy; delete[] cEz;

    // fields
    gpuErrchk(cudaFree(dEx));gpuErrchk(cudaFree(dEy));gpuErrchk(cudaFree(dEz));
    gpuErrchk(cudaFree(dHx));gpuErrchk(cudaFree(dHy));gpuErrchk(cudaFree(dHz));
    // materials
    gpuErrchk(cudaFree(depsx));gpuErrchk(cudaFree(depsy));gpuErrchk(cudaFree(depsz));
    // PML
    gpuErrchk(cudaFree(dpml_Exy0));cudaFree(dpml_Exy1);cudaFree(dpml_Exz0);cudaFree(dpml_Exz1);
    cudaFree(dpml_Eyx0);cudaFree(dpml_Eyx1);cudaFree(dpml_Eyz0);cudaFree(dpml_Eyz1);
    cudaFree(dpml_Ezx0);cudaFree(dpml_Ezx1);cudaFree(dpml_Ezy0);cudaFree(dpml_Ezy1);
    cudaFree(dpml_Hxy0);cudaFree(dpml_Hxy1);cudaFree(dpml_Hxz0);cudaFree(dpml_Hxz1);
    cudaFree(dpml_Hyx0);cudaFree(dpml_Hyx1);cudaFree(dpml_Hyz0);cudaFree(dpml_Hyz1);
    cudaFree(dpml_Hzx0);cudaFree(dpml_Hzx1);cudaFree(dpml_Hzy0);cudaFree(dpml_Hzy1);

    gpuErrchk(cudaFree(dkappa_H_x));cudaFree(dkappa_H_y);cudaFree(dkappa_H_z);
    cudaFree(dkappa_E_x);cudaFree(dkappa_E_y);cudaFree(dkappa_E_z);
    gpuErrchk(cudaFree(dbHx));cudaFree(dbHy);cudaFree(dbHz);
    cudaFree(dbEx);cudaFree(dbEy);cudaFree(dbEz);
    cudaFree(dcHx);cudaFree(dcHy);cudaFree(dcHz);
    cudaFree(dcEx);cudaFree(dcEy);cudaFree(dcEz);

    // source index arrays
    gpuErrchk(cudaFree(di0s)); cudaFree(dj0s); cudaFree(dk0s);
    cudaFree(dIs); cudaFree(dJs); cudaFree(dKs);

    // source arrays
    gpuErrchk(cudaFree(dJx)); cudaFree(dJy); cudaFree(dJz);
    cudaFree(dMx); cudaFree(dMy); cudaFree(dMz);

    return 0;
}
