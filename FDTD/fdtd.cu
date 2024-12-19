/******************************************************************************
 * Copyright (c) 2023, Andrew Michaels.  All rights reserved.
 * Copyright (c) 2023, NVIDIA CORPORATION.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

#include "fdtd.hpp"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <ctime>
#include <cstdio>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <thrust/complex.h>

#include <cuda_runtime.h>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <fstream>
#include <signal.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

static std::mutex mutex_H, mutex_E, mutex_conv, mutex_conv2, mutex_t1, mutex_t2;
static std::condition_variable cv_H, cv_E, cv_conv, cv_conv2, cv_t1, cv_t2;
std::atomic<size_t> counter_H(0), counter_E(0), counter_conv(0), counter_conv2(0), counter_t1(0), counter_t2(0);
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

void signal_callback_handler(int signum)
{
    std::cout << "Caught signal " << signum << std::endl;
    exit(signum);
}

inline double norm2(double *mat, int len)
{
    double sum = 0.0;
    for(int i=0;i<len;i++)
    {
        sum += mat[i]*mat[i];
    }
    return sqrt(sum);
}

inline double norm2(double *mat, double *mat2, int len)
{
    double sum = 0.0;
    for(int i=0;i<len;i++)
    {
        sum += (mat[i]-mat2[i])*(mat[i]-mat2[i]);
    }
    return sqrt(sum);
}

__host__ __device__ inline double calc_phase(double t0, double t1, double t2, double f0, double f1, double f2)
{
    double f10=f1-f0, f21=f2-f1;
    double ret{0.0f};
    if(f10==0 && f21==0) return 0.0;
    else{
        ret = atan2(f10*(sin(t2)-sin(t1))-f21*(sin(t1)-sin(t0)), f21*(cos(t1)-cos(t0))-f10*(cos(t2)-cos(t1)));
        if(ret!=ret) return M_PI/2;
        else return ret;
    } 
}

__host__ __device__ inline double calc_amplitude(double t0, double t1, double t2, double f0, double f1, double f2, double phase)
{
    double f10=f1-f0, f21=f2-f1;
    double ret;
    if(f10==0 && f21==0) return 0.0;
    else if(f21*f21 >= f10*f10){
        ret = f21/(cos(phase)*(sin(t2)-sin(t1))+sin(phase)*(cos(t2)-cos(t1)));
        if(ret!=ret) return 0.0;
        else return ret;
    }else{
        ret = f10/(cos(phase)*(sin(t1)-sin(t0))+sin(phase)*(cos(t1)-cos(t0)));
        if(ret!=ret) return 0.0;
        else return ret;
    }
}

void prepare_nvlink (int device_id, int gpus_count)
{
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
                // No gpuErrchk here to avoid unnecessary "peer access already enabled" error/exit
				cudaDeviceEnablePeerAccess (other_device_id, 0);
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
				std::cout << "P2P capable on device " << device_id << std::endl;
		}
	}
}

void cudaP2PAssert_pair(int* P2Pworking_pair, int src_device_id, int dst_device_id)
{
    int Ngpus = 2;
    double *host_data = new double[Ngpus]{2.6, 7.1};
    double *host_buffer = new double[Ngpus]{-1.0, -1.0};
    int *device_id_list = new int[Ngpus]{src_device_id, dst_device_id};
    std::vector<std::thread> threads_list;
    threads_list.reserve(Ngpus);
    double** d_sender_list = new double*[Ngpus];
    double** d_receiver_list = new double*[Ngpus];
    double *d_sender, *d_receiver;
    for(int device_id=0; device_id<Ngpus; device_id++){
        gpuErrchk(cudaSetDevice(device_id_list[device_id]));
        gpuErrchk(cudaMalloc((void **)&d_sender, sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&d_receiver, sizeof(double)));
        gpuErrchk(cudaMemcpy(d_sender, host_data+device_id, sizeof(double), cudaMemcpyHostToDevice));
        d_sender_list[device_id] = d_sender;
        d_receiver_list[device_id] = d_receiver;
    }    
    for(int device_id=0; device_id<Ngpus; device_id++){
        threads_list.emplace_back([=](){
            gpuErrchk(cudaSetDevice(device_id_list[device_id]));
            if(device_id==0) {
                gpuErrchk(cudaMemcpyPeer(d_receiver_list[device_id], device_id_list[device_id], d_sender_list[device_id+1], device_id_list[device_id+1], sizeof(double)));
                gpuErrchk(cudaMemcpy(host_buffer+device_id, d_receiver_list[device_id], sizeof(double), cudaMemcpyDeviceToHost));
            }
            else {
                gpuErrchk(cudaMemcpyPeer(d_receiver_list[device_id], device_id_list[device_id], d_sender_list[device_id-1], device_id_list[device_id-1], sizeof(double)));
                gpuErrchk(cudaMemcpy(host_buffer+device_id, d_receiver_list[device_id], sizeof(double), cudaMemcpyDeviceToHost));
            }
        });
    }
    for(auto &thread: threads_list)
        thread.join();

    if(host_buffer[0]==host_data[1] && host_buffer[1]==host_data[0])
        *P2Pworking_pair = 1;
    else
        *P2Pworking_pair = 0;

    for(int device_id=0; device_id<Ngpus; device_id++){
        gpuErrchk(cudaFree(d_sender_list[device_id]));
        gpuErrchk(cudaFree(d_receiver_list[device_id]));
    }
}

void cudaP2PAssert(int* _P2Pworking, int Ngpus)
{
    int P2Pworking_pair{0};
    *_P2Pworking = 1;
    for(int src=0; src<Ngpus; src++){
        for(int dst=0; dst<Ngpus; dst++){
            if(src!=dst){
                cudaP2PAssert_pair(&P2Pworking_pair, src, dst);
                if(!P2Pworking_pair){
                    std::cout << "P2P not working between devices " << src << " and " << dst << std::endl;
                    *_P2Pworking = 0;
                }
            }
        }
    }
}

template <typename T>
__global__ void init_mat(T *mat, size_t len, T val)
{
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind_global >= len) return; 
    mat[ind_global] = val;
}

__global__ void kernel_copy_eps(kernelpar *kpar)
{
//     1D grid of 1D blocks
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind_global >= kpar->I*kpar->J*kpar->K) return;

    int i, j, k;
    i = ind_global/(kpar->J*kpar->K);
    j = (ind_global%(kpar->J*kpar->K))/kpar->K;
    k = (ind_global%(kpar->J*kpar->K))%kpar->K;

    kpar->epsx[ind_global] = kpar->epsx_full[(i+kpar->i0)*kpar->Ny*kpar->Nx + (j+kpar->j0)*kpar->Nx + k + kpar->k0];
    kpar->epsy[ind_global] = kpar->epsy_full[(i+kpar->i0)*kpar->Ny*kpar->Nx + (j+kpar->j0)*kpar->Nx + k + kpar->k0];
    kpar->epsz[ind_global] = kpar->epsz_full[(i+kpar->i0)*kpar->Ny*kpar->Nx + (j+kpar->j0)*kpar->Nx + k + kpar->k0];
}

__global__ void kernel_calc_complex_fields(double t0, double t1, double t2, double *F_t0, double *F_t1, double *F_t2, size_t len, complex128 *F_out)
{
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind_global >= len) return;
    double f0, f1, f2, phi, A;
    f0 = F_t0[ind_global];
    f1 = F_t1[ind_global];
    f2 = F_t2[ind_global];
    phi = calc_phase(t0, t1, t2, f0, f1, f2);
    A = calc_amplitude(t0, t1, t2, f0, f1, f2, phi);
    if(A<0) {A*=-1; phi+=M_PI;}
    F_out[ind_global].real = A*cos(phi);
    F_out[ind_global].imag = -A*sin(phi);

}

__global__ void kernel_update_H(kernelpar *kpar)
{
//     1D grid of 1D blocks
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
//     1D grid of 2D blocks
//     size_t ind_global = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
//     2D grid of 1D blocks
//     size_t ind_global = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;
    if (ind_global >= kpar->size){
        return;
    }
    int i, j, k, ind_ijk, ind_ijp1k, ind_ip1jk, ind_ijkp1;
    double dExdy, dExdz, dEydx, dEydz, dEzdx, dEzdy;
    int ind_pml, ind_pml_param;
    double b, C, kappa;
    int //srcind_ijk,
        srcind_src, srcind_size, srcind_size_offset;
    double src_t;

    // ind_global = i*_J*_K + j*_K + k;
    i = ind_global/(kpar->J*kpar->K);
    j = (ind_global%(kpar->J*kpar->K))/kpar->K;
    k = (ind_global%(kpar->J*kpar->K))%kpar->K;

    //blockDim.x = 128; gridDim.x=585
//     k = ind_global/(kpar->I*kpar->J);
//     j = (ind_global%(kpar->I*kpar->J))/kpar->I;
//     i = (ind_global%(kpar->I*kpar->J))%kpar->I;

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
        if( i >= kpar->i0s[ii] && j >= kpar->j0s[ii] && k >= kpar->k0s[ii]
            && i < kpar->i0s[ii]+kpar->Is[ii] && j < kpar->j0s[ii]+kpar->Js[ii] && k < kpar->k0s[ii]+kpar->Ks[ii]){
            srcind_src = (i-kpar->i0s[ii])*kpar->Js[ii]*kpar->Ks[ii] + (j-kpar->j0s[ii])*kpar->Ks[ii] + k-kpar->k0s[ii];
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
}

__global__ void kernel_update_E(kernelpar *kpar)
{
//     1D grid of 1D blocks
    size_t ind_global = blockIdx.x * blockDim.x + threadIdx.x;
//     1D grid of 2D blocks
//     size_t ind_global = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
//     2D grid of 1D blocks
//     size_t ind_global = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;

    if (ind_global >= kpar->size){
        return;
    }
    int i, j, k, ind_ijk, ind_ijm1k, ind_im1jk, ind_ijkm1;
    double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy,
           b_x, b_y, b_z;
    int ind_pml, ind_pml_param;
    double b, C, kappa;
    int srcind_src,
        //srcind_ijk, srcind_global,
        srcind_size, srcind_size_offset;
    double src_t;

    // ind_global = i*_J*_K + j*_K + k;
    i = ind_global/(kpar->J*kpar->K);
    j = (ind_global%(kpar->J*kpar->K))/kpar->K;
    k = (ind_global%(kpar->J*kpar->K))%kpar->K;

    //blockDim.x = 128; gridDim.x=585
    // K=585, J=572, I=87
//     k = ind_global/(kpar->I*kpar->J);
//     j = (ind_global%(kpar->I*kpar->J))/kpar->I;
//     i = (ind_global%(kpar->I*kpar->J))%kpar->I;

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
    // kernel/domain's ind_global = i*_J*_K + j*_K + k;
    srcind_size_offset = 0;
    for(int ii = 0; ii < kpar->srclen; ii++){
        srcind_size = kpar->Is[ii] * kpar->Js[ii] * kpar->Ks[ii];
        if(i >= kpar->i0s[ii] && j >= kpar->j0s[ii] && k >= kpar->k0s[ii]
           && i < kpar->i0s[ii]+kpar->Is[ii] && j < kpar->j0s[ii]+kpar->Js[ii] && k < kpar->k0s[ii]+kpar->Ks[ii]){
            srcind_src = (i-kpar->i0s[ii])*kpar->Js[ii]*kpar->Ks[ii] + (j-kpar->j0s[ii])*kpar->Ks[ii] + k - kpar->k0s[ii];
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
    int ind_border = i*kpar->J*kpar->K + j*kpar->K + k;
    b_x = kpar->dt/kpar->epsx[ind_border].real;
    b_y = kpar->dt/kpar->epsy[ind_border].real;
    b_z = kpar->dt/kpar->epsz[ind_border].real;
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
                kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - src_t*kpar->Jx[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsx[ind_border].real;
                src_t = sin(kpar->t + kpar->Jy[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - src_t*kpar->Jy[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsy[ind_border].real;
                src_t = sin(kpar->t + kpar->Jz[srcind_src+srcind_size_offset].imag)*((1+kpar->src_min)*exp(-(kpar->t-kpar->src_T)*(kpar->t-kpar->src_T)/kpar->src_k)-kpar->src_min);
                kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - src_t*kpar->Jz[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsz[ind_border].real;
            }
            else{
                src_t = sin(kpar->t + kpar->Jx[srcind_src+srcind_size_offset].imag);
                kpar->Ex[ind_ijk] = kpar->Ex[ind_ijk] - src_t*kpar->Jx[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsx[ind_border].real;
                src_t = sin(kpar->t + kpar->Jy[srcind_src+srcind_size_offset].imag);
                kpar->Ey[ind_ijk] = kpar->Ey[ind_ijk] - src_t*kpar->Jy[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsy[ind_border].real;
                src_t = sin(kpar->t + kpar->Jz[srcind_src+srcind_size_offset].imag);
                kpar->Ez[ind_ijk] = kpar->Ez[ind_ijk] - src_t*kpar->Jz[srcind_src+srcind_size_offset].real*kpar->dt/kpar->epsz[ind_border].real;
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

__global__ void kernel_calc_ydAx(size_t size, size_t Nx, size_t Ny, size_t Nz, size_t i0, size_t i1,
                thrust::complex<double> *ydAx,
                thrust::complex<double> *Ex_adj, thrust::complex<double> *Ey_adj, thrust::complex<double> *Ez_adj,
                thrust::complex<double> *Ex_fwd, thrust::complex<double> *Ey_fwd, thrust::complex<double> *Ez_fwd,
                thrust::complex<double> *epsx0, thrust::complex<double> *epsy0, thrust::complex<double> *epsz0,
                thrust::complex<double> *epsxp, thrust::complex<double> *epsyp, thrust::complex<double> *epszp)
{
    int i,j,k;

    k = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    i = blockIdx.z;
    if( k >=Nx || j>=Ny || i>=Nz ) { return; }
    size_t ind = (i1-i0)*Nx*Ny + j *Nx+k;
    size_t offset = Nx * Ny * i0;
    if(ind >= size) { return; }

    ydAx[ind] = Ex_adj[ind] * Ex_fwd[ind] * (epsxp[ind+offset]-epsx0[ind+offset]) +
                Ey_adj[ind] * Ey_fwd[ind] * (epsyp[ind+offset]-epsy0[ind+offset]) +
                Ez_adj[ind] * Ez_fwd[ind] * (epszp[ind+offset]-epsz0[ind+offset]);
}

void fdtd::FDTD::calc_ydAx(size_t size, size_t Nx, size_t Ny, size_t Nz, size_t i0, size_t i1, size_t i2,
                std::complex<double> *ydAx,
                std::complex<double> *Ex_adj, std::complex<double> *Ey_adj, std::complex<double> *Ez_adj,
                std::complex<double> *Ex_fwd, std::complex<double> *Ey_fwd, std::complex<double> *Ez_fwd,
                std::complex<double> *epsx0, std::complex<double> *epsy0, std::complex<double> *epsz0,
                std::complex<double> *epsxp, std::complex<double> *epsyp, std::complex<double> *epszp)
{
    thrust::complex<double> *dydAx, *dEx_adj, *dEy_adj, *dEz_adj, *dEx_fwd, *dEy_fwd, *dEz_fwd,
                            *depsx0, *depsy0, *depsz0, *depsxp, *depsyp, *depszp;
    size_t sizefull = Nx * Ny * Nz;
    gpuErrchk(cudaMalloc((void **)&dydAx, size * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dEx_adj, size * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dEy_adj, size * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dEz_adj, size * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dEx_fwd, size * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dEy_fwd, size * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dEz_fwd, size * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&depsx0, sizefull * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&depsy0, sizefull * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&depsz0, sizefull * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&depsxp, sizefull * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&depsyp, sizefull * sizeof(std::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&depszp, sizefull * sizeof(std::complex<double>)));

    gpuErrchk(cudaMemcpy(dEx_adj, Ex_adj, size * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dEy_adj, Ey_adj, size * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dEz_adj, Ez_adj, size * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dEx_fwd, Ex_fwd, size * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dEy_fwd, Ey_fwd, size * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dEz_fwd, Ez_fwd, size * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(depsx0, epsx0, sizefull * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(depsy0, epsy0, sizefull * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(depsz0, epsz0, sizefull * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(depsxp, epsxp, sizefull * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(depsyp, epsyp, sizefull * sizeof(std::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(depszp, epszp, sizefull * sizeof(std::complex<double>), cudaMemcpyHostToDevice));

//    kernel_calc_ydAx <<< ceil(size/128.0), 128 >>> (size, NxNy, i0, i1, dydAx, dEx_adj, dEy_adj, dEz_adj, dEx_fwd, dEy_fwd, dEz_fwd,
//                depsx0, depsy0, depsz0, depsxp, depsyp, depszp);

    dim3 block_dim(128,1);
    dim3 grid_dim((int)ceil(Nx/128.0), (int)ceil(Ny/1.0));
    for(int i=i1; i<i2; i++){
        kernel_calc_ydAx <<< grid_dim, block_dim >>> (size, Nx, Ny, Nz, i1, i, dydAx, dEx_adj, dEy_adj, dEz_adj, dEx_fwd, dEy_fwd, dEz_fwd,
                depsx0, depsy0, depsz0, depsxp, depsyp, depszp);
    }

    cudaError_t cudaerr = cudaDeviceSynchronize();
    if(cudaerr != cudaSuccess){ printf("fdtd::calc_ydAx: kernel launch status \"%s\".\n", cudaGetErrorString(cudaerr)); }

    gpuErrchk(cudaMemcpy(ydAx, dydAx, size * sizeof(std::complex<double>), cudaMemcpyDeviceToHost));

    gpuErrchk(cudaFree(dydAx));
    gpuErrchk(cudaFree(dEx_adj)); gpuErrchk(cudaFree(dEy_adj)); gpuErrchk(cudaFree(dEz_adj));
    gpuErrchk(cudaFree(dEx_fwd)); gpuErrchk(cudaFree(dEy_fwd)); gpuErrchk(cudaFree(dEz_fwd));
    gpuErrchk(cudaFree(depsx0)); gpuErrchk(cudaFree(depsy0)); gpuErrchk(cudaFree(depsz0));
    gpuErrchk(cudaFree(depsxp)); gpuErrchk(cudaFree(depsyp)); gpuErrchk(cudaFree(depszp));
}

fdtd::FDTD::FDTD()
{
    // make sure all of our PML arrays start NULL
    _pml_Exy0 = NULL; _pml_Exy1 = NULL; _pml_Exz0 = NULL; _pml_Exz1 = NULL;
    _pml_Eyx0 = NULL; _pml_Eyx1 = NULL; _pml_Eyz0 = NULL; _pml_Eyz1 = NULL;
    _pml_Ezx0 = NULL; _pml_Ezx1 = NULL; _pml_Ezy0 = NULL; _pml_Ezy1 = NULL;
    _pml_Hxy0 = NULL; _pml_Hxy1 = NULL; _pml_Hxz0 = NULL; _pml_Hxz1 = NULL;
    _pml_Hyx0 = NULL; _pml_Hyx1 = NULL; _pml_Hyz0 = NULL; _pml_Hyz1 = NULL;
    _pml_Hzx0 = NULL; _pml_Hzx1 = NULL; _pml_Hzy0 = NULL; _pml_Hzy1 = NULL;
    
    _kappa_H_x = NULL; _kappa_H_y = NULL; _kappa_H_z = NULL;
    _kappa_E_x = NULL; _kappa_E_y = NULL; _kappa_E_z = NULL;

    _bHx = NULL; _bHy = NULL; _bHz = NULL;
    _bEx = NULL; _bEy = NULL; _bEz = NULL;

    _cHx = NULL; _cHy = NULL; _cHz = NULL;
    _cEx = NULL; _cEy = NULL; _cEz = NULL;

    _w_pml_x0 = 0; _w_pml_x1 = 0;
    _w_pml_y0 = 0; _w_pml_y1 = 0;
    _w_pml_z0 = 0; _w_pml_z1 = 0;

    _complex_eps = false;

    _kpar_host = (kernelpar *)malloc(sizeof(kernelpar));
    _kpar_host->srclen = 0;
    _srclen = 0;
    // kernel parameter structures
    gpuErrchk(cudaMalloc((void **)&_kpar_device, sizeof(kernelpar)));

}

fdtd::FDTD::~FDTD()
{
    // Clean up PML arrays
    delete[] _pml_Exy0; delete[] _pml_Exy1; delete[] _pml_Exz0; delete[] _pml_Exz1;
    delete[] _pml_Eyx0; delete[] _pml_Eyx1; delete[] _pml_Eyz0; delete[] _pml_Eyz1;
    delete[] _pml_Ezx0; delete[] _pml_Ezx1; delete[] _pml_Ezy0; delete[] _pml_Ezy1;
    delete[] _pml_Hxy0; delete[] _pml_Hxy1; delete[] _pml_Hxz0; delete[] _pml_Hxz1;
    delete[] _pml_Hyx0; delete[] _pml_Hyx1; delete[] _pml_Hyz0; delete[] _pml_Hyz1;
    delete[] _pml_Hzx0; delete[] _pml_Hzx1; delete[] _pml_Hzy0; delete[] _pml_Hzy1;

    delete [] _kappa_H_x;
    delete [] _kappa_H_y;
    delete [] _kappa_H_z;

    delete [] _kappa_E_x;
    delete [] _kappa_E_y;
    delete [] _kappa_E_z;

    delete [] _bHx;
    delete [] _bHy;
    delete [] _bHz;

    delete [] _bEx;
    delete [] _bEy;
    delete [] _bEz;

    delete [] _cHx;
    delete [] _cHy;
    delete [] _cHz;

    delete [] _cEx;
    delete [] _cEy;
    delete [] _cEz;

    delete [] _kpar_host;

    gpuErrchk(cudaFree(_kpar_device));

}

void fdtd::FDTD::set_physical_dims(double X, double Y, double Z,
                                         double dx, double dy, double dz)
{
    _X = X; _Y = Y; _Z = Z;
    _dx = dx; _dy = dy; _dz = dz;
}

void fdtd::FDTD::set_grid_dims(int Nx, int Ny, int Nz)
{
    _Nx = Nx;
    _Ny = Ny;
    _Nz = Nz;
}


void fdtd::FDTD::set_local_grid(int k0, int j0, int i0, int K, int J, int I)
{
    _i0 = i0; _j0 = j0; _k0 = k0;
    _I = I; _J = J; _K = K;

}

void fdtd::FDTD::set_local_grid_perturb(int i1, int i2)
{
    _i1 = i1; _i2 = i2;
}

void fdtd::FDTD::set_wavelength(double wavelength)
{
    _wavelength = wavelength;
    _R = _wavelength/(2*M_PI);
}

void fdtd::FDTD::set_rtol(double rtol)
{
    _rtol = rtol;
}

void fdtd::FDTD::set_domain_decomposition(char domain_decomp)
{
    _domain_decomp = domain_decomp;
}

void fdtd::FDTD::set_dt(double dt)
{
    _dt = dt;
    _odt = 1.0/_dt;
}

void fdtd::FDTD::set_complex_eps(bool complex_eps)
{
    _complex_eps = complex_eps;
}

void fdtd::FDTD::set_gpus_count(int gpus_count)
{
    gpuErrchk(cudaGetDeviceCount(&_Ngpus));
    _gpus_count = (gpus_count<_Ngpus)?gpus_count:_Ngpus;
    std::cout << "Use " << _gpus_count << " of " << _Ngpus << " available GPUs." << std::endl;
}

void fdtd::FDTD::set_Ncycle(double Ncycle)
{
    _Ncycle = Ncycle;
}

void fdtd::FDTD::copyCUDA_field_arrays()
{
    size_t size = (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2);
    gpuErrchk(cudaMemcpy(_Ex, _kpar_host->Ex, size * sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(_Ey, _kpar_host->Ey, size * sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(_Ez, _kpar_host->Ez, size * sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(_Hx, _kpar_host->Hx, size * sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(_Hy, _kpar_host->Hy, size * sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(_Hz, _kpar_host->Hz, size * sizeof(double), cudaMemcpyDeviceToHost));
}

void fdtd::FDTD::block_CUDA_free()
{
    // fields
    cudaFree(dEx);cudaFree(dEy);cudaFree(dEz);cudaFree(dHx);cudaFree(dHy);cudaFree(dHz);
    // materials
    cudaFree(depsx);cudaFree(depsy);cudaFree(depsz);
    // PML
    cudaFree(dpml_Exy0);cudaFree(dpml_Exy1);cudaFree(dpml_Exz0);cudaFree(dpml_Exz1);
    cudaFree(dpml_Eyx0);cudaFree(dpml_Eyx1);cudaFree(dpml_Eyz0);cudaFree(dpml_Eyz1);
    cudaFree(dpml_Ezx0);cudaFree(dpml_Ezx1);cudaFree(dpml_Ezy0);cudaFree(dpml_Ezy1);
    cudaFree(dpml_Hxy0);cudaFree(dpml_Hxy1);cudaFree(dpml_Hxz0);cudaFree(dpml_Hxz1);
    cudaFree(dpml_Hyx0);cudaFree(dpml_Hyx1);cudaFree(dpml_Hyz0);cudaFree(dpml_Hyz1);
    cudaFree(dpml_Hzx0);cudaFree(dpml_Hzx1);cudaFree(dpml_Hzy0);cudaFree(dpml_Hzy1);

    cudaFree(dkappa_H_x);cudaFree(dkappa_H_y);cudaFree(dkappa_H_z);
    cudaFree(dkappa_E_x);cudaFree(dkappa_E_y);cudaFree(dkappa_E_z);
    cudaFree(dbHx);cudaFree(dbHy);cudaFree(dbHz);
    cudaFree(dbEx);cudaFree(dbEy);cudaFree(dbEz);
    cudaFree(dcHx);cudaFree(dcHy);cudaFree(dcHz);
    cudaFree(dcEx);cudaFree(dcEy);cudaFree(dcEz);

}

void fdtd::FDTD::block_CUDA_src_free()
{
    // source index arrays
    gpuErrchk(cudaFree(di0s));
    cudaFree(dj0s);
    cudaFree(dk0s);
    cudaFree(dIs);
    cudaFree(dJs);
    cudaFree(dKs);

    // source arrays
    gpuErrchk(cudaFree(dJx));
    cudaFree(dJy);
    cudaFree(dJz);
    cudaFree(dMx);
    cudaFree(dMy);
    cudaFree(dMz);
}

void fdtd::FDTD::block_CUDA_src_malloc_memcpy()
{
    // extract the list of tuples of the sources
    _i0s = (int *)malloc(_srclen * sizeof(int));
    _j0s = (int *)malloc(_srclen * sizeof(int));
    _k0s = (int *)malloc(_srclen * sizeof(int));
    _Is = (int *)malloc(_srclen * sizeof(int));
    _Js = (int *)malloc(_srclen * sizeof(int));
    _Ks = (int *)malloc(_srclen * sizeof(int));
    size_t size = 0,
           size_offset = 0,
           sizeall = 0;
    for(int i=0; i < _srclen; i++){
        _i0s[i] = _sources[i].i0;
        _j0s[i] = _sources[i].j0;
        _k0s[i] = _sources[i].k0;
        _Is[i] = _sources[i].I;
        _Js[i] = _sources[i].J;
        _Ks[i] = _sources[i].K;
        sizeall += _Is[i] * _Js[i] * _Ks[i];
    }

    // initialize GPU memory to store source arrays
    gpuErrchk(cudaMalloc((void **)&dJx, sizeall * sizeof(complex128)));
    cudaMalloc((void **)&dJy, sizeall * sizeof(complex128));
    cudaMalloc((void **)&dJz, sizeall * sizeof(complex128));
    cudaMalloc((void **)&dMx, sizeall * sizeof(complex128));
    cudaMalloc((void **)&dMy, sizeall * sizeof(complex128));
    cudaMalloc((void **)&dMz, sizeall * sizeof(complex128));

    gpuErrchk(cudaMalloc((void **)&di0s, _srclen * sizeof(int)));
    cudaMalloc((void **)&dj0s, _srclen * sizeof(int));
    cudaMalloc((void **)&dk0s, _srclen * sizeof(int));
    cudaMalloc((void **)&dIs, _srclen * sizeof(int));
    cudaMalloc((void **)&dJs, _srclen * sizeof(int));
    cudaMalloc((void **)&dKs, _srclen * sizeof(int));

    gpuErrchk(cudaMemcpy(di0s, _i0s, _srclen * sizeof(int), cudaMemcpyHostToDevice));
    cudaMemcpy(dj0s, _j0s, _srclen * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dk0s, _k0s, _srclen * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dIs, _Is, _srclen * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dJs, _Js, _srclen * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dKs, _Ks, _srclen * sizeof(int), cudaMemcpyHostToDevice);

    size_offset = 0;
    for(int i=0; i < _srclen; i++){
        size = _Is[i] * _Js[i] * _Ks[i];
        gpuErrchk(cudaMemcpy(dJx+size_offset, _sources[i].Jx, size * sizeof(complex128), cudaMemcpyHostToDevice));
        cudaMemcpy(dJy+size_offset, _sources[i].Jy, size * sizeof(complex128), cudaMemcpyHostToDevice);
        cudaMemcpy(dJz+size_offset, _sources[i].Jz, size * sizeof(complex128), cudaMemcpyHostToDevice);
        cudaMemcpy(dMx+size_offset, _sources[i].Mx, size * sizeof(complex128), cudaMemcpyHostToDevice);
        cudaMemcpy(dMy+size_offset, _sources[i].My, size * sizeof(complex128), cudaMemcpyHostToDevice);
        cudaMemcpy(dMz+size_offset, _sources[i].Mz, size * sizeof(complex128), cudaMemcpyHostToDevice);
        size_offset += size;
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

}

void fdtd::FDTD::block_CUDA_multigpu_init()
{
    // size and indices of subdomains
    int It{}, Jt{}, Kt{}, i0t{}, j0t{}, k0t{};
    int N{0};
    int Npmlx = _w_pml_x0 + _w_pml_x1,
        Npmly = _w_pml_y0 + _w_pml_y1,
        Npmlz = _w_pml_z0 + _w_pml_z1;
    int pml_xmin = _w_pml_x0, pml_xmax = _Nx-_w_pml_x1,
        pml_ymin = _w_pml_y0, pml_ymax = _Ny-_w_pml_y1,
        pml_zmin = _w_pml_z0, pml_zmax = _Nz-_w_pml_z1;

    _i0s = (int *)malloc(_srclen*sizeof(int));
    _j0s = (int *)malloc(_srclen*sizeof(int));
    _k0s = (int *)malloc(_srclen*sizeof(int));
    _Is = (int *)malloc(_srclen*sizeof(int));
    _Js = (int *)malloc(_srclen*sizeof(int));
    _Ks = (int *)malloc(_srclen*sizeof(int));
    size_t srcsize=0, size_offset=0, sizeall=0;
    for(int i=0; i<_srclen; i++){
        _i0s[i] = _sources[i].i0;
        _j0s[i] = _sources[i].j0;
        _k0s[i] = _sources[i].k0;
        _Is[i] = _sources[i].I;
        _Js[i] = _sources[i].J;
        _Ks[i] = _sources[i].K;
        sizeall += _Is[i] * _Js[i] * _Ks[i];
    }

    _kpar_list.reserve(_gpus_count);
    _kpar_list_host = new kernelpar*[_gpus_count];
    for(int device_id=0; device_id<_gpus_count; device_id++)
    {
        // create kpar for each device
        gpuErrchk(cudaSetDevice(device_id));
        _kpar_host = (kernelpar *)malloc(sizeof(kernelpar));
        gpuErrchk(cudaMalloc((void **)&_kpar_device, sizeof(kernelpar)));
        _kpar_list.push_back(_kpar_device);

        // BCs
        gpuErrchk(cudaMalloc((void **)&dbc, 3 * sizeof(char)));
        gpuErrchk(cudaMemcpy(dbc, _bc, 3 * sizeof(char), cudaMemcpyHostToDevice));
        _kpar_host->bc = dbc;

        // Domain decomposition
        switch(_domain_decomp){
            case 'x':
                i0t = _i0; j0t = _j0; k0t = device_id*_Nx/_gpus_count;
                It = _I; Jt = _J; Kt = (device_id==_gpus_count-1)?_Nx-device_id*_Nx/_gpus_count:(device_id+1)*_Nx/_gpus_count-device_id*_Nx/_gpus_count;
                _ghost_size = _I*_J; _buffer_size = (k0t==0 || k0t+Kt==_Nx)?_I*_J:2*_I*_J;
                break;
            case 'y':
                i0t = _i0; j0t = device_id*_Ny/_gpus_count; k0t = _k0;
                It = _I; Jt = (device_id==_gpus_count-1)?_Ny-device_id*_Ny/_gpus_count:(device_id+1)*_Ny/_gpus_count-device_id*_Ny/_gpus_count; Kt = _K;
                _ghost_size = _I*_K; _buffer_size = (j0t==0 || j0t+Jt==_Ny)?_I*_K:2*_I*_K;
                break;
            case 'z':
                i0t = device_id*_Nz/_gpus_count; j0t = _j0; k0t = _k0;
                It = (device_id==_gpus_count-1)?_Nz-device_id*_Nz/_gpus_count:(device_id+1)*_Nz/_gpus_count-device_id*_Nz/_gpus_count; Jt = _J; Kt = _K;
                _ghost_size = _J*_K; _buffer_size = (i0t==0 || i0t+It==_Nz)?_J*_K:2*_J*_K;
                break;
            default:
                std::cout << "Unsupported domain decomposition method.  Default to X-cut." << std::endl;
                i0t = _i0; j0t = _j0; k0t = device_id*_Nx/_gpus_count;
                It = _I; Jt = _J; Kt = (device_id==_gpus_count-1)?_Nx-device_id*_Nx/_gpus_count:(device_id+1)*_Nx/_gpus_count-device_id*_Nx/_gpus_count;
                _ghost_size = _I*_J; _buffer_size = (k0t==0 || k0t+Kt==_Nx)?_I*_J:2*_I*_J;
                break;
        }
        _ghostsize_list.push_back(_ghost_size);
        _buffersize_list.push_back(_buffer_size);

        _kpar_host->I = It; _kpar_host->J = Jt; _kpar_host->K = Kt;
        _kpar_host->i0 = i0t; _kpar_host->j0 = j0t; _kpar_host->k0 = k0t;
        _kpar_host->size = _kpar_host->I * _kpar_host->J * _kpar_host->K;
        _kpar_host->odx = _R/_dx;
        _kpar_host->ghost_size = _ghost_size;
        _kpar_host->buffer_size = _buffer_size;

        std::cout << "Device " << device_id << ", (I, i0)=(" << _kpar_host->I << ", " << _kpar_host->i0 << "), (J, j0)=(" 
                << _kpar_host->J << ", " << _kpar_host->j0 << "), (K, k0)=(" << _kpar_host->K << ", " << _kpar_host->k0 << ")" << std::endl;

        _griddim_list.push_back(ceil(_kpar_host->I*_kpar_host->J*_kpar_host->K/(double)(blockdim)));
        _griddimghost_list.push_back(ceil(_kpar_host->buffer_size/(double)(blockdim)));

        _kpar_host->Nx = _Nx; _kpar_host->Ny = _Ny; _kpar_host->Nz = _Nz;
        _kpar_host->dt = _dt; 

        // fields
        size_t field_len = (_kpar_host->I+2)*(_kpar_host->J+2)*(_kpar_host->K+2);
        gpuErrchk(cudaMalloc((void **)&dEx, field_len * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dEy, field_len * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dEz, field_len * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dHx, field_len * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dHy, field_len * sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dHz, field_len * sizeof(double)));

        gpuErrchk(cudaMalloc((void **)&dE_bufferleft, 3*_ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dE_bufferright, 3*_ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_bufferleft, 3*_ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_bufferright, 3*_ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dE_ghostleft, 3*_ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dE_ghostright, 3*_ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_ghostleft, 3*_ghost_size*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dH_ghostright, 3*_ghost_size*sizeof(double)));

        size_t numBlocks = ceil(field_len/(double)blockdim);
        init_mat <<< numBlocks, blockdim >>> (dEx, field_len, 0.0);
        init_mat <<< numBlocks, blockdim >>> (dEy, field_len, 0.0);
        init_mat <<< numBlocks, blockdim >>> (dEz, field_len, 0.0);
        init_mat <<< numBlocks, blockdim >>> (dHx, field_len, 0.0);
        init_mat <<< numBlocks, blockdim >>> (dHy, field_len, 0.0);
        init_mat <<< numBlocks, blockdim >>> (dHz, field_len, 0.0);

        _kpar_host->E_bufferleft = dE_bufferleft; 
        _kpar_host->E_bufferright = dE_bufferright;
        _kpar_host->H_bufferleft = dH_bufferleft;
        _kpar_host->H_bufferright = dH_bufferright;
        _kpar_host->E_ghostleft = dE_ghostleft;
        _kpar_host->E_ghostright = dE_ghostright;
        _kpar_host->H_ghostleft = dH_ghostleft;
        _kpar_host->H_ghostright = dH_ghostright;

        _kpar_host->Ex = dEx; _kpar_host->Ey = dEy; _kpar_host->Ez = dEz;
        _kpar_host->Hx = dHx; _kpar_host->Hy = dHy; _kpar_host->Hz = dHz;

        // materials
        size_t mat_len = _kpar_host->I*_kpar_host->J*_kpar_host->K;
        gpuErrchk(cudaMalloc((void **)&depsx, mat_len*sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&depsy, mat_len*sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&depsz, mat_len*sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&depsx_full, _Nx*_Ny*_Nz*sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&depsy_full, _Nx*_Ny*_Nz*sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&depsz_full, _Nx*_Ny*_Nz*sizeof(complex128)));
        gpuErrchk(cudaMemcpy(depsx_full, _eps_x, _Nx*_Ny*_Nz*sizeof(complex128), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(depsy_full, _eps_y, _Nx*_Ny*_Nz*sizeof(complex128), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(depsz_full, _eps_z, _Nx*_Ny*_Nz*sizeof(complex128), cudaMemcpyHostToDevice));

        _kpar_host->epsx = depsx; _kpar_host->epsy = depsy; _kpar_host->epsz = depsz;
        _kpar_host->epsx_full = depsx_full; _kpar_host->epsy_full = depsy_full; _kpar_host->epsz_full = depsz_full;
        
        gpuErrchk(cudaMemcpy(_kpar_list[device_id], _kpar_host, sizeof(kernelpar), cudaMemcpyHostToDevice));
        numBlocks = ceil(mat_len/(double)blockdim);
        kernel_copy_eps <<< numBlocks, blockdim >>> (_kpar_list[device_id]);

        gpuErrchk(cudaFree(depsx_full)); gpuErrchk(cudaFree(depsy_full)); gpuErrchk(cudaFree(depsz_full));

        // PML
        if(k0t < _w_pml_x0){
            N = _kpar_host->I * _kpar_host->J * (_w_pml_x0 - _kpar_host->k0);
            gpuErrchk(cudaMalloc((void **)&dpml_Eyx0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezx0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyx0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzx0, N*sizeof(double)));
            numBlocks = ceil(N/(double)blockdim);
            init_mat <<< numBlocks, blockdim >>> (dpml_Eyx0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Ezx0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hyx0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hzx0, N, 0.0);
        }
        if(k0t+Kt>_Nx-_w_pml_x1){
            N = _kpar_host->I * _kpar_host->J * (_kpar_host->k0 + _kpar_host->K - (_Nx - _w_pml_x1));
            gpuErrchk(cudaMalloc((void **)&dpml_Eyx1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezx1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyx1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzx1, N*sizeof(double)));
            numBlocks = ceil(N/(double)blockdim);
            init_mat <<< numBlocks, blockdim >>> (dpml_Eyx1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Ezx1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hyx1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hzx1, N, 0.0);
        }
        if(j0t<_w_pml_y0){
            N = _kpar_host->I * _kpar_host->K * (_w_pml_y0 - _kpar_host->j0);
            gpuErrchk(cudaMalloc((void **)&dpml_Exy0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezy0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxy0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzy0, N*sizeof(double)));
            numBlocks = ceil(N/(double)blockdim);
            init_mat <<< numBlocks, blockdim >>> (dpml_Exy0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Ezy0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hxy0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hzy0, N, 0.0);
        }
        if(j0t+Jt>_Ny-_w_pml_y1){
            N = _kpar_host->I * _kpar_host->K * (_kpar_host->j0 + _kpar_host->J - (_Ny - _w_pml_y1));
            gpuErrchk(cudaMalloc((void **)&dpml_Exy1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Ezy1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxy1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hzy1, N*sizeof(double)));
            numBlocks = ceil(N/(double)blockdim);
            init_mat <<< numBlocks, blockdim >>> (dpml_Exy1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Ezy1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hxy1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hzy1, N, 0.0);
        }
        if(i0t<_w_pml_z0){
            N = _kpar_host->J * _kpar_host->K * (_w_pml_z0 - _kpar_host->i0);
            gpuErrchk(cudaMalloc((void **)&dpml_Exz0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Eyz0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxz0, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyz0, N*sizeof(double)));
            numBlocks = ceil(N/(double)blockdim);
            init_mat <<< numBlocks, blockdim >>> (dpml_Exz0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Eyz0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hxz0, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hyz0, N, 0.0);
        }
        if(i0t+It>_Nz-_w_pml_z1){
            N = _kpar_host->J * _kpar_host->K * (_kpar_host->i0 + _kpar_host->I - (_Nz - _w_pml_z1));
            gpuErrchk(cudaMalloc((void **)&dpml_Exz1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Eyz1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hxz1, N*sizeof(double)));
            gpuErrchk(cudaMalloc((void **)&dpml_Hyz1, N*sizeof(double)));
            numBlocks = ceil(N/(double)blockdim);
            init_mat <<< numBlocks, blockdim >>> (dpml_Exz1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Eyz1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hxz1, N, 0.0);
            init_mat <<< numBlocks, blockdim >>> (dpml_Hyz1, N, 0.0);
        }

        gpuErrchk(cudaMalloc((void **)&dkappa_H_x, Npmlx*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_H_y, Npmly*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_H_z, Npmlz*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_E_x, Npmlx*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_E_y, Npmly*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dkappa_E_z, Npmlz*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbHx, Npmlx*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbHy, Npmly*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbHz, Npmlz*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbEx, Npmlx*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbEy, Npmly*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dbEz, Npmlz*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcHx, Npmlx*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcHy, Npmly*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcHz, Npmlz*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcEx, Npmlx*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcEy, Npmly*sizeof(double)));
        gpuErrchk(cudaMalloc((void **)&dcEz, Npmlz*sizeof(double)));

        gpuErrchk(cudaMemcpy(dkappa_H_x, _kappa_H_x, Npmlx*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_H_y, _kappa_H_y, Npmly*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_H_z, _kappa_H_z, Npmlz*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_E_x, _kappa_E_x, Npmlx*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_E_y, _kappa_E_y, Npmly*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dkappa_E_z, _kappa_E_z, Npmlz*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbHx, _bHx, Npmlx*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbHy, _bHy, Npmly*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbHz, _bHz, Npmlz*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbEx, _bEx, Npmlx*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbEy, _bEy, Npmly*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dbEz, _bEz, Npmlz*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcHx, _cHx, Npmlx*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcHy, _cHy, Npmly*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcHz, _cHz, Npmlz*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcEx, _cEx, Npmlx*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcEy, _cEy, Npmly*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dcEz, _cEz, Npmlz*sizeof(double), cudaMemcpyHostToDevice));

        _kpar_host->pml_Exy0 = dpml_Exy0; _kpar_host->pml_Exy1 = dpml_Exy1;
        _kpar_host->pml_Exz0 = dpml_Exz0; _kpar_host->pml_Exz1 = dpml_Exz1;
        _kpar_host->pml_Eyx0 = dpml_Eyx0; _kpar_host->pml_Eyx1 = dpml_Eyx1;
        _kpar_host->pml_Eyz0 = dpml_Eyz0; _kpar_host->pml_Eyz1 = dpml_Eyz1;
        _kpar_host->pml_Ezx0 = dpml_Ezx0; _kpar_host->pml_Ezx1 = dpml_Ezx1;
        _kpar_host->pml_Ezy0 = dpml_Ezy0; _kpar_host->pml_Ezy1 = dpml_Ezy1;

        _kpar_host->pml_Hxy0 = dpml_Hxy0; _kpar_host->pml_Hxy1 = dpml_Hxy1;
        _kpar_host->pml_Hxz0 = dpml_Hxz0; _kpar_host->pml_Hxz1 = dpml_Hxz1;
        _kpar_host->pml_Hyx0 = dpml_Hyx0; _kpar_host->pml_Hyx1 = dpml_Hyx1;
        _kpar_host->pml_Hyz0 = dpml_Hyz0; _kpar_host->pml_Hyz1 = dpml_Hyz1;
        _kpar_host->pml_Hzx0 = dpml_Hzx0; _kpar_host->pml_Hzx1 = dpml_Hzx1;
        _kpar_host->pml_Hzy0 = dpml_Hzy0; _kpar_host->pml_Hzy1 = dpml_Hzy1;

        _kpar_host->kappa_H_x = dkappa_H_x; _kpar_host->kappa_H_y = dkappa_H_y; _kpar_host->kappa_H_z = dkappa_H_z;
        _kpar_host->kappa_E_x = dkappa_E_x; _kpar_host->kappa_E_y = dkappa_E_y; _kpar_host->kappa_E_z = dkappa_E_z;
        _kpar_host->bHx = dbHx; _kpar_host->bHy = dbHy; _kpar_host->bHz = dbHz;
        _kpar_host->bEx = dbEx; _kpar_host->bEy = dbEy; _kpar_host->bEz = dbEz;
        _kpar_host->cHx = dcHx; _kpar_host->cHy = dcHy; _kpar_host->cHz = dcHz;
        _kpar_host->cEx = dcEx; _kpar_host->cEy = dcEy; _kpar_host->cEz = dcEz;

        _kpar_host->w_pml_x0 = _w_pml_x0; _kpar_host->w_pml_x1 = _w_pml_x1;
        _kpar_host->w_pml_y0 = _w_pml_y0; _kpar_host->w_pml_y1 = _w_pml_y1;
        _kpar_host->w_pml_z0 = _w_pml_z0; _kpar_host->w_pml_z1 = _w_pml_z1;
        _kpar_host->pml_xmin = pml_xmin;  _kpar_host->pml_xmax = pml_xmax;
        _kpar_host->pml_ymin = pml_ymin;  _kpar_host->pml_ymax = pml_ymax;
        _kpar_host->pml_zmin = pml_zmin;  _kpar_host->pml_zmax = pml_zmax;

        // sources
        _kpar_host->srclen = _srclen;
        gpuErrchk(cudaMalloc((void **)&dJx, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dJy, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dJz, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dMx, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dMy, sizeall * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dMz, sizeall * sizeof(complex128)));

        gpuErrchk(cudaMalloc((void **)&di0s, _srclen * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dj0s, _srclen * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dk0s, _srclen * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dIs, _srclen * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dJs, _srclen * sizeof(complex128)));
        gpuErrchk(cudaMalloc((void **)&dKs, _srclen * sizeof(complex128)));

        gpuErrchk(cudaMemcpy(di0s, _i0s, _srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dj0s, _j0s, _srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dk0s, _k0s, _srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dIs, _Is, _srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dJs, _Js, _srclen * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dKs, _Ks, _srclen * sizeof(int), cudaMemcpyHostToDevice));

        size_offset = 0;
        for(int i=0; i<_srclen;i++){
            srcsize = _Is[i]*_Js[i]*_Ks[i];
            gpuErrchk(cudaMemcpy(dJx+size_offset, _sources[i].Jx, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dJy+size_offset, _sources[i].Jy, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dJz+size_offset, _sources[i].Jz, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dMx+size_offset, _sources[i].Mx, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dMy+size_offset, _sources[i].My, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(dMz+size_offset, _sources[i].Mz, srcsize * sizeof(complex128), cudaMemcpyHostToDevice));
            size_offset += srcsize;
        }

        _kpar_host->i0s = di0s; _kpar_host->j0s = dj0s; _kpar_host->k0s = dk0s;
        _kpar_host->Is = dIs; _kpar_host->Js = dJs; _kpar_host->Ks = dKs;
        _kpar_host->Jx = dJx; _kpar_host->Jy = dJy; _kpar_host->Jz = dJz;
        _kpar_host->Mx = dMx; _kpar_host->My = dMy; _kpar_host->Mz = dMz;
        _kpar_host->src_T = _src_T; _kpar_host->src_min = _src_min; _kpar_host->src_k = _src_k;

        gpuErrchk(cudaMemcpy(_kpar_list[device_id], _kpar_host, sizeof(kernelpar), cudaMemcpyHostToDevice));
        _kpar_list_host[device_id] = _kpar_host;
    }

    // Free up memory
    delete[] _i0s; delete[] _j0s; delete[] _k0s; delete[] _Is; delete[] _Js; delete[] _Ks;
}

void fdtd::FDTD::block_CUDA_multigpu_free()
{
    int i0t{}, j0t{}, k0t{}, It{}, Jt{}, Kt{};
    for(int device_id=0; device_id<_gpus_count; device_id++){
        cudaSetDevice(device_id);

        k0t = _kpar_list_host[device_id]->k0;
        j0t = _kpar_list_host[device_id]->j0;
        i0t = _kpar_list_host[device_id]->i0;
        Kt = _kpar_list_host[device_id]->K;
        Jt = _kpar_list_host[device_id]->J;
        It = _kpar_list_host[device_id]->I;

        gpuErrchk(cudaFree(_kpar_list_host[device_id]->bc));
        // sources
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->i0s));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->j0s));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->k0s));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Is));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Js));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Ks));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Jx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Jy));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Jz));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Mx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->My));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Mz));
        // fields
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Ex));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Ey));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Ez));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Hx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Hy));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->Hz));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->E_bufferleft));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->E_bufferright));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->H_bufferleft));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->H_bufferright));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->E_ghostleft));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->E_ghostright));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->H_ghostleft));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->H_ghostright));
        // materials
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->epsx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->epsy));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->epsz));
        // PMLs
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->kappa_H_x));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->kappa_H_y));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->kappa_H_z));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->kappa_E_x));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->kappa_E_y));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->kappa_E_z));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->bHx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->bHy));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->bHz));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->bEx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->bEy));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->bEz));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->cHx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->cHy));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->cHz));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->cEx));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->cEy));
        gpuErrchk(cudaFree(_kpar_list_host[device_id]->cEz));

        if(k0t < _w_pml_x0){
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Eyx0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Ezx0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hyx0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hzx0));
        }
        if(k0t+Kt>_Nx-_w_pml_x1){
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Eyx1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Ezx1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hyx1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hzx1));
        }

        if(j0t<_w_pml_y0){
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Exy0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Ezy0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hxy0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hzy0));
        }
        if(j0t+Jt>_Ny-_w_pml_y1){
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Exy1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Ezy1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hxy1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hzy1));
        }
        if(i0t<_w_pml_z0){
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Exz0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Eyz0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hxz0));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hyz0));
        }
        if(i0t+It>_Nz-_w_pml_z1){
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Exz1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Eyz1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hxz1));
            gpuErrchk(cudaFree(_kpar_list_host[device_id]->pml_Hyz1));
        }
        // gpuErrchk(cudaFree(_kpar_list[device_id]));
    }
}

void fdtd::FDTD::block_CUDA_malloc_memcpy()
{
     // BCs
    gpuErrchk(cudaMalloc((void **)&dbc, 3 * sizeof(char)));
    gpuErrchk(cudaMemcpy(dbc, _bc, 3 * sizeof(char), cudaMemcpyHostToDevice));
    _kpar_host->bc = dbc;

    // fields
    size_t size = (_I+2)*(_J+2)*(_K+2);
    gpuErrchk(cudaMalloc((void **)&dEx, size * sizeof(double)));
    cudaMalloc((void **)&dEy, size * sizeof(double));
    cudaMalloc((void **)&dEz, size * sizeof(double));
    cudaMalloc((void **)&dHx, size * sizeof(double));
    cudaMalloc((void **)&dHy, size * sizeof(double));
    cudaMalloc((void **)&dHz, size * sizeof(double));

    gpuErrchk(cudaMemcpy(dEx, _Ex, size * sizeof(double), cudaMemcpyHostToDevice));
    cudaMemcpy(dEy, _Ey, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dEz, _Ez, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dHx, _Hx, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dHy, _Hy, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dHz, _Hz, size * sizeof(double), cudaMemcpyHostToDevice);

    _kpar_host->Ex = dEx;
    _kpar_host->Ey = dEy;
    _kpar_host->Ez = dEz;
    _kpar_host->Hx = dHx;
    _kpar_host->Hy = dHy;
    _kpar_host->Hz = dHz;

    _kpar_host->I = _I;
    _kpar_host->J = _J;
    _kpar_host->K = _K;
    _kpar_host->i0 = _i0;
    _kpar_host->j0 = _j0;
    _kpar_host->k0 = _k0;
    _kpar_host->Nx = _Nx;
    _kpar_host->Ny = _Ny;
    _kpar_host->Nz = _Nz;
    _kpar_host->size = _I*_J*_K;

    _kpar_host->dt = _dt;

    // materials
    size = _I * _J * _K;
    gpuErrchk(cudaMalloc((void **)&depsx, size * sizeof(complex128)));
    cudaMalloc((void **)&depsy, size * sizeof(complex128));
    cudaMalloc((void **)&depsz, size * sizeof(complex128));

    gpuErrchk(cudaMemcpy(depsx, _eps_x, size * sizeof(complex128), cudaMemcpyHostToDevice));
    cudaMemcpy(depsy, _eps_y, size * sizeof(complex128), cudaMemcpyHostToDevice);
    cudaMemcpy(depsz, _eps_z, size * sizeof(complex128), cudaMemcpyHostToDevice);

    _kpar_host->epsx = depsx;
    _kpar_host->epsy = depsy;
    _kpar_host->epsz = depsz;

    // PML
    int N,
        xmin = _w_pml_x0, xmax = _Nx-_w_pml_x1,
        ymin = _w_pml_y0, ymax = _Ny-_w_pml_y1,
        zmin = _w_pml_z0, zmax = _Nz-_w_pml_z1;
    int Npmlx = _w_pml_x0 + _w_pml_x1,
        Npmly = _w_pml_y0 + _w_pml_y1,
        Npmlz = _w_pml_z0 + _w_pml_z1;

    // touches xmin boudary
    if(_k0 < xmin) {
        N = _I * _J * (xmin - _k0);
        gpuErrchk(cudaMalloc((void **)&dpml_Eyx0, N * sizeof(double)));
        cudaMalloc((void **)&dpml_Ezx0, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hyx0, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hzx0, N * sizeof(double));
        gpuErrchk(cudaMemset(dpml_Eyx0, 0.0, N * sizeof(double)));
        cudaMemset(dpml_Ezx0, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hyx0, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hzx0, 0.0, N * sizeof(double));
    }
    // touches xmax boundary
    if(_k0 +_K > xmax) {
        N = _I * _J * (_k0  + _K - xmax);
        gpuErrchk(cudaMalloc((void **)&dpml_Eyx1, N * sizeof(double)));
        cudaMalloc((void **)&dpml_Ezx1, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hyx1, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hzx1, N * sizeof(double));
        gpuErrchk(cudaMemset(dpml_Eyx1, 0.0, N * sizeof(double)));
        cudaMemset(dpml_Ezx1, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hyx1, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hzx1, 0.0, N * sizeof(double));
    }
    // touches ymin boundary
    if(_j0 < ymin) {
        N = _I * _K * (ymin - _j0);
        gpuErrchk(cudaMalloc((void **)&dpml_Exy0, N * sizeof(double)));
        cudaMalloc((void **)&dpml_Ezy0, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hxy0, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hzy0, N * sizeof(double));
        gpuErrchk(cudaMemset(dpml_Exy0, 0.0, N * sizeof(double)));
        cudaMemset(dpml_Ezy0, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hxy0, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hzy0, 0.0, N * sizeof(double));
    }
    // touches ymax boundary
    if(_j0 + _J > ymax) {
        N = _I * _K * (_j0 + _J - ymax);
        gpuErrchk(cudaMalloc((void **)&dpml_Exy1, N * sizeof(double)));
        cudaMalloc((void **)&dpml_Ezy1, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hxy1, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hzy1, N * sizeof(double));
        gpuErrchk(cudaMemset(dpml_Exy1, 0.0, N * sizeof(double)));
        cudaMemset(dpml_Ezy1, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hxy1, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hzy1, 0.0, N * sizeof(double));
    }
    // touches zmin boundary
    if(_i0 < zmin) {
        N = _J * _K * (zmin - _i0);
        gpuErrchk(cudaMalloc((void **)&dpml_Exz0, N * sizeof(double)));
        cudaMalloc((void **)&dpml_Eyz0, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hxz0, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hyz0, N * sizeof(double));
        gpuErrchk(cudaMemset(dpml_Exz0, 0.0, N * sizeof(double)));
        cudaMemset(dpml_Eyz0, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hxz0, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hyz0, 0.0, N * sizeof(double));
    }
    // touches zmax boundary
    if(_i0 + _I > zmax) {
        N = _J * _K * (_i0 + _I - zmax);
        gpuErrchk(cudaMalloc((void **)&dpml_Exz1, N * sizeof(double)));
        cudaMalloc((void **)&dpml_Eyz1, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hxz1, N * sizeof(double));
        cudaMalloc((void **)&dpml_Hyz1, N * sizeof(double));
        gpuErrchk(cudaMemset(dpml_Exz1, 0.0, N * sizeof(double)));
        cudaMemset(dpml_Eyz1, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hxz1, 0.0, N * sizeof(double));
        cudaMemset(dpml_Hyz1, 0.0, N * sizeof(double));
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

    gpuErrchk(cudaMemcpy(dkappa_H_x, _kappa_H_x, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dkappa_H_y, _kappa_H_y, Npmly * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dkappa_H_z, _kappa_H_z, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dkappa_E_x, _kappa_E_x, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dkappa_E_y, _kappa_E_y, Npmly * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dkappa_E_z, _kappa_E_z, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dbHx, _bHx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dbHy, _bHy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dbHz, _bHz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dbEx, _bEx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dbEy, _bEy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dbEz, _bEz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dcHx, _cHx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dcHy, _cHy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dcHz, _cHz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dcEx, _cEx, Npmlx * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dcEy, _cEy, Npmly * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dcEz, _cEz, Npmlz * sizeof(double), cudaMemcpyHostToDevice));

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

    _kpar_host->w_pml_x0 = _w_pml_x0;
    _kpar_host->w_pml_x1 = _w_pml_x1;
    _kpar_host->w_pml_y0 = _w_pml_y0;
    _kpar_host->w_pml_y1 = _w_pml_y1;
    _kpar_host->w_pml_z0 = _w_pml_z0;
    _kpar_host->w_pml_z1 = _w_pml_z1;

}

void fdtd::FDTD::set_field_arrays(double *Ex, double *Ey, double *Ez,
                                  double *Hx, double *Hy, double *Hz)
{
    _Ex = Ex; _Ey = Ey; _Ez = Ez;
    _Hx = Hx; _Hy = Hy; _Hz = Hz;
}

void fdtd::FDTD::set_mat_arrays(complex128 *eps_x, complex128 *eps_y, complex128 *eps_z)
{
    _eps_x = eps_x; _eps_y = eps_y; _eps_z = eps_z;
//    _mu_x = mu_x; _mu_y = mu_y; _mu_z = mu_z;

}

void fdtd::FDTD::update_H(int n, double t)
{
    double odx = _R/_dx,
           ody = _R/_dy,
           odz = _R/_dz;

    int pml_xmin = _w_pml_x0, pml_xmax = _Nx-_w_pml_x1,
        pml_ymin = _w_pml_y0, pml_ymax = _Ny-_w_pml_y1,
        pml_zmin = _w_pml_z0, pml_zmax = _Nz-_w_pml_z1;

    _kpar_host->pml_xmin = pml_xmin;
    _kpar_host->pml_xmax = pml_xmax;
    _kpar_host->pml_ymin = pml_ymin;
    _kpar_host->pml_ymax = pml_ymax;
    _kpar_host->pml_zmin = pml_zmin;
    _kpar_host->pml_zmax = pml_zmax;
    _kpar_host->odx = odx;
    _kpar_host->t = t;
    _kpar_host->src_T = _src_T;
    _kpar_host->src_min = _src_min;
    _kpar_host->src_k = _src_k;

    gpuErrchk(cudaMemcpy(_kpar_device, _kpar_host, sizeof(kernelpar), cudaMemcpyHostToDevice));
    kernel_update_H <<< ceil(_I*_J*_K/128.0), 128 >>> (_kpar_device);
}

void fdtd::FDTD::set_GPUDirect()
{
    // Detect GPU Direct status
    int perfRank{0};
    bool P2P_AllEnabled{true};
    for(int src_id=0; src_id<_gpus_count; src_id++){
        for(int dst_id=src_id+1; dst_id<_gpus_count; dst_id++){
            gpuErrchk(cudaDeviceGetP2PAttribute(&perfRank, cudaDevP2PAttrAccessSupported, src_id, dst_id));
            P2P_AllEnabled = P2P_AllEnabled && perfRank;
        }
    }
    if(!P2P_AllEnabled){
        std::cout << "P2P not enabled. Preparing GPU Direct..." << std::endl;
    }

    // Prepare GPU Direct
    std::vector<std::thread> threads_nvlink;
    threads_nvlink.reserve(_gpus_count);
    for(int device_id=0; device_id<_gpus_count; device_id++)
    {
        threads_nvlink.emplace_back([=](){
            try{ prepare_nvlink(device_id, _gpus_count); }
            catch(std::runtime_error &error){ std::cerr << "Error in thread " << device_id << ":" << error.what() << std::endl; }
        });
    }
    for(auto &thread: threads_nvlink)
        thread.join();

    // Validate whether GPU Direct works 
    if(_gpus_count>1) {
        cudaP2PAssert(&_P2Pworking, _gpus_count);
        if(!_P2Pworking){
            std::cout << "P2P not working. Fallback to host-to-device copy." << std::endl;
        }
    }
}

void fdtd::FDTD::solve()
{
    // Prepare streams and events
    std::vector<cudaStream_t> compute_stream_list, sync_stream_list;
    std::vector<cudaEvent_t> e_border_event_list, h_border_event_list, e_bulk_event_list, h_bulk_event_list, e_time_event_list, h_time_event_list;
    compute_stream_list.reserve(_gpus_count);
    sync_stream_list.reserve(_gpus_count);
    e_border_event_list.reserve(_gpus_count);
    h_border_event_list.reserve(_gpus_count);
    e_bulk_event_list.reserve(_gpus_count);
    h_bulk_event_list.reserve(_gpus_count);
    e_time_event_list.reserve(_gpus_count);
    h_time_event_list.reserve(_gpus_count);
    int least_priority{}, greatest_priority{};
    gpuErrchk(cudaDeviceGetStreamPriorityRange(&least_priority, &greatest_priority));
    for(int device_id=0; device_id<_gpus_count; device_id++)
    {
        cudaSetDevice(device_id);
        cudaStream_t compute_stream, sync_stream;
        gpuErrchk(cudaStreamCreateWithPriority(&compute_stream, cudaStreamDefault, least_priority));
        gpuErrchk(cudaStreamCreateWithPriority(&sync_stream, cudaStreamDefault, greatest_priority));
        compute_stream_list.push_back(compute_stream);
        sync_stream_list.push_back(sync_stream);

        cudaEvent_t h_bulk_computed, e_bulk_computed, h_border_computed, e_border_computed, h_time_computed, e_time_computed;
        gpuErrchk(cudaEventCreateWithFlags(&h_bulk_computed, cudaEventDisableTiming));
        gpuErrchk(cudaEventCreateWithFlags(&e_bulk_computed, cudaEventDisableTiming));
        gpuErrchk(cudaEventCreateWithFlags(&h_border_computed, cudaEventDisableTiming));
        gpuErrchk(cudaEventCreateWithFlags(&e_border_computed, cudaEventDisableTiming));
        gpuErrchk(cudaEventCreateWithFlags(&h_time_computed, cudaEventDisableTiming));
        gpuErrchk(cudaEventCreateWithFlags(&e_time_computed, cudaEventDisableTiming));
        h_bulk_event_list.push_back(h_bulk_computed);
        e_bulk_event_list.push_back(e_bulk_computed);
        h_border_event_list.push_back(h_border_computed);
        e_border_event_list.push_back(e_border_computed);
        h_time_event_list.push_back(h_time_computed);
        e_time_event_list.push_back(e_time_computed);
    }

    // Prepare convergence check matrices
    double *Ex0, *Ex1, *Ex2, *Ey0, *Ey1, *Ey2, *Ez0, *Ez1, *Ez2, *phi0, *phi1, *A0, *A1;
    Ex0 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ex0, Ex0+_Nconv, 0.0);
    Ex1 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ex1, Ex1+_Nconv, 0.0);
    Ex2 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ex2, Ex2+_Nconv, 0.0);
    Ey0 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ey0, Ey0+_Nconv, 0.0);
    Ey1 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ey1, Ey1+_Nconv, 0.0);
    Ey2 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ey2, Ey2+_Nconv, 0.0);
    Ez0 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ez0, Ez0+_Nconv, 0.0);
    Ez1 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ez1, Ez1+_Nconv, 0.0);
    Ez2 = (double *)malloc(_Nconv*sizeof(double)); std::fill(Ez2, Ez2+_Nconv, 0.0);
    phi0 = (double *)malloc(_Nconv*sizeof(double)); std::fill(phi0, phi0+_Nconv, 0.0);
    phi1 = (double *)malloc(_Nconv*sizeof(double)); std::fill(phi1, phi1+_Nconv, 0.0);
    A0 = (double *)malloc(_Nconv*sizeof(double)); std::fill(A0, A0+_Nconv, 0.0);
    A1 = (double *)malloc(_Nconv*sizeof(double)); std::fill(A1, A1+_Nconv, 0.0);

    // Prepare buffer matrices for device-to-device transfer 
    std::vector<double *> Emid_list{};
    Emid_list.reserve(_gpus_count);
    for(int device_id=0; device_id<_gpus_count; device_id++){
        double *Emid = (double *)malloc(3*_ghostsize_list[device_id]*sizeof(double));
        std::fill(Emid, Emid+3*_ghostsize_list[device_id], 0.0);
        Emid_list.push_back(Emid);
    }

    double *A_change = new double{1.0f};
    double *phi_change = new double{1.0f};
    int *Tn_factor = new int{300};
    *Tn_factor = int(sqrt(pow(_X, 2) + pow(_Y, 2) +pow(_Z, 2))*12.2);
    int Tn = int(_Ncycle*3/4);
    double amp_rtol{_rtol}, phi_rtol{sqrt(_rtol)};

    // FDTD iterations
    std::vector<std::thread> threads_kernels;
    threads_kernels.reserve(_gpus_count);
    signal(SIGINT, signal_callback_handler);

    auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Start FDTD:" << ctime(&timenow);

    for(int device_id=0; device_id<_gpus_count; device_id++)
    {
        threads_kernels.emplace_back([=](){
            try{
                gpuErrchk(cudaSetDevice(device_id));
                int p{0}, conv_idx{0};
                double t0{0.0f}, t1{0.0f}, t2{0.0f};
                double phasex{0.0f}, phasey{0.0f}, phasez{0.0f}, ampx{0.0f}, ampy{0.0f}, ampz{0.0f};
                int It, Jt, Kt, i0t, j0t, k0t, I, J, K;
                int ci, cj, ck;
                It = _kpar_list_host[device_id]->I; Jt = _kpar_list_host[device_id]->J; Kt = _kpar_list_host[device_id]->K;
                i0t = _kpar_list_host[device_id]->i0; j0t = _kpar_list_host[device_id]->j0; k0t = _kpar_list_host[device_id]->k0;
                I = _kpar_list_host[device_id]->Nz; J = _kpar_list_host[device_id]->Ny; K = _kpar_list_host[device_id]->Nx;
                double *Ext, *Eyt, *Ezt, *Hxt, *Hyt, *Hzt;
                Ext = (double *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(double));
                Eyt = (double *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(double));
                Ezt = (double *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(double));
                Hxt = (double *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(double));
                Hyt = (double *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(double));
                Hzt = (double *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(double));

                int step{0};
                while(*A_change > amp_rtol || *phi_change > phi_rtol || std::isnan(*A_change) || std::isnan(*phi_change))
                {
                    _kpar_list_host[device_id]->t = step*_dt;
                    gpuErrchk(cudaMemcpyAsync(_kpar_list[device_id], _kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(h_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], e_border_event_list[device_id], 0);
                    kernel_update_H_bulk <<< _griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (_kpar_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], h_time_event_list[device_id], 0);
                    kernel_update_H_border <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                    cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);

                    if(_gpus_count>1){
                        Barrier(counter_H, mutex_H, cv_H, _gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferright, device_id, _kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->H_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else if(device_id==_gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferleft, device_id, _kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->H_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferright, device_id, _kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferleft, device_id, _kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->H_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->H_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        kernel_load_ghostlayerH <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                        cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);
                    }

                    _kpar_list_host[device_id]->t = (step+0.5)*_dt;
                    gpuErrchk(cudaMemcpyAsync(_kpar_list[device_id], _kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(e_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], h_border_event_list[device_id], 0);
                    kernel_update_E_bulk <<< _griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (_kpar_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], e_time_event_list[device_id], 0);
                    kernel_update_E_border <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                    cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    
                    if(_gpus_count>1){
                        Barrier(counter_E, mutex_E, cv_E, _gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferright, device_id, _kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->E_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else if(device_id==_gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferleft, device_id, _kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->E_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferright, device_id, _kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferleft, device_id, _kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->E_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->E_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        kernel_load_ghostlayerE <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                        cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    }
                    
                    if(p==Tn* (*Tn_factor)-1){
                        p = 0; // reset counter
                        t0=t1; t1=t2; t2= step*_dt;
                        gpuErrchk(cudaMemcpy(Ext, _kpar_list_host[device_id]->Ex, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                        gpuErrchk(cudaMemcpy(Eyt, _kpar_list_host[device_id]->Ey, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                        gpuErrchk(cudaMemcpy(Ezt, _kpar_list_host[device_id]->Ez, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                        gpuErrchk(cudaMemcpy(Hxt, _kpar_list_host[device_id]->Hx, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                        gpuErrchk(cudaMemcpy(Hyt, _kpar_list_host[device_id]->Hy, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                        gpuErrchk(cudaMemcpy(Hzt, _kpar_list_host[device_id]->Hz, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));

                        // gpuErrchk(cudaMemcpy(epsxt, _kpar_list_host[device_id]->epsx, It*Jt*Kt*sizeof(complex128), cudaMemcpyDeviceToHost));
                        gpuErrchk(cudaDeviceSynchronize());

                        for(int ii=0; ii<_Nconv; ii++){
                            conv_idx = ii * int(I*J*K/_Nconv);
                            ci = conv_idx/(J*K);
                            cj = (conv_idx%(J*K))/K;
                            ck = (conv_idx%(J*K))%K;
                            if(ci>=i0t && ci<i0t+It && cj>=j0t && cj<j0t+Jt && ck>=k0t && ck<k0t+Kt){
                                A0[ii] = A1[ii];
                                phi0[ii] = phi1[ii];
                                Ex0[ii] = Ex1[ii];
                                Ex1[ii] = Ex2[ii];
                                Ex2[ii] = Ext[(ci-i0t+1)*(Jt+2)*(Kt+2)+(cj-j0t+1)*(Kt+2)+ck-k0t+1];
                                Ey0[ii] = Ey1[ii];
                                Ey1[ii] = Ey2[ii];
                                Ey2[ii] = Eyt[(ci-i0t+1)*(Jt+2)*(Kt+2)+(cj-j0t+1)*(Kt+2)+ck-k0t+1];
                                Ez0[ii] = Ez1[ii];
                                Ez1[ii] = Ez2[ii];
                                Ez2[ii] = Ezt[(ci-i0t+1)*(Jt+2)*(Kt+2)+(cj-j0t+1)*(Kt+2)+ck-k0t+1];

                                phasex = calc_phase(t0, t1, t2, Ex0[ii], Ex1[ii], Ex2[ii]);
                                phasey = calc_phase(t0, t1, t2, Ey0[ii], Ey1[ii], Ey2[ii]);
                                phasez = calc_phase(t0, t1, t2, Ez0[ii], Ez1[ii], Ez2[ii]);
                                ampx = calc_amplitude(t0, t1, t2, Ex0[ii], Ex1[ii], Ex2[ii], phasex);
                                ampy = calc_amplitude(t0, t1, t2, Ey0[ii], Ey1[ii], Ey2[ii], phasey);
                                ampz = calc_amplitude(t0, t1, t2, Ez0[ii], Ez1[ii], Ez2[ii], phasez);
                                if(ampx<0){ ampx *= -1.0; phasex += M_PI; }
                                if(ampy<0){ ampy *= -1.0; phasey += M_PI; }
                                if(ampz<0){ ampz *= -1.0; phasez += M_PI; }
                                phi1[ii] = phasex + phasey + phasez;
                                A1[ii] = ampx + ampy + ampz;
                            }
                        }
                        Barrier(counter_conv, mutex_conv, cv_conv, _gpus_count);

                        if(device_id==0){
                            if(norm2(A0, _Nconv)<1e-12 || norm2(phi0, _Nconv)<1e-12 ){
                                *A_change = 1;
                                *phi_change = 1;
                            }else{
                                *A_change = norm2(A1, A0, _Nconv)/norm2(A0, _Nconv);
                                *phi_change = norm2(phi1, phi0, _Nconv)/norm2(phi0, _Nconv);
                            }
                            auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                            std::cout << "Step:" << step << ", A_change:" << *A_change << ", phi_change:" << *phi_change << ", Time:" << ctime(&timenow);

                            // if (*A_change<=2.0){
                            *Tn_factor = int(exp(log(*A_change)*0.30103+3.950225));
                            if(*Tn_factor<1) *Tn_factor = 1;
                            // }
                        }
                        Barrier(counter_conv2, mutex_conv2, cv_conv2, _gpus_count);
                    }else{
                        p++;
                    }
                    step++;
                }  // while end

                // FDTD_capture_t0_fields
                // for(int ii=0; ii<It; ii++){
                //     for(int jj=0; jj<Jt; jj++){
                //         for(int kk=0; kk<Kt; kk++){
                //             _Ex_t0[(ii+i0t)*J*K+(jj+j0t)*K+kk+k0t] = Ext[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+kk+1];
                //             _Ey_t0[(ii+i0t)*J*K+(jj+j0t)*K+kk+k0t] = Eyt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+kk+1];
                //             _Ez_t0[(ii+i0t)*J*K+(jj+j0t)*K+kk+k0t] = Ezt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+kk+1];
                //             _Hx_t0[(ii+i0t)*J*K+(jj+j0t)*K+kk+k0t] = Hxt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+kk+1];
                //             _Hy_t0[(ii+i0t)*J*K+(jj+j0t)*K+kk+k0t] = Hyt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+kk+1];
                //             _Hz_t0[(ii+i0t)*J*K+(jj+j0t)*K+kk+k0t] = Hzt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+kk+1];
                //         }
                //     }
                // } // FDTD_capture_t0_fields end

                // capture t0 fields on GPUs
                double *dEx_t0, *dEy_t0, *dEz_t0, *dHx_t0, *dHy_t0, *dHz_t0;
                gpuErrchk(cudaMalloc((void **)&dEx_t0, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dEy_t0, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dEz_t0, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHx_t0, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHy_t0, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHz_t0, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                // gpuErrchk(cudaMemcpy(dEx_t0, Ext, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dEy_t0, Eyt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dEz_t0, Ezt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHx_t0, Hxt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHy_t0, Hyt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHz_t0, Hzt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                gpuErrchk(cudaMemcpy(dEx_t0, _kpar_list_host[device_id]->Ex, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dEy_t0, _kpar_list_host[device_id]->Ey, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dEz_t0, _kpar_list_host[device_id]->Ez, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHx_t0, _kpar_list_host[device_id]->Hx, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHy_t0, _kpar_list_host[device_id]->Hy, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHz_t0, _kpar_list_host[device_id]->Hz, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));

                // Perform another Tn steps to get a second time point at t1
                for(int step1=0; step1<Tn; step1++){
                    _kpar_list_host[device_id]->t = (step+step1)*_dt;

                    gpuErrchk(cudaMemcpyAsync(_kpar_list[device_id], _kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(h_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], e_border_event_list[device_id], 0);
                    kernel_update_H_bulk <<< _griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (_kpar_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], h_time_event_list[device_id], 0);
                    kernel_update_H_border <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                    cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);

                    if(_gpus_count>1){
                        Barrier(counter_H, mutex_H, cv_H, _gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferright, device_id, _kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->H_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else if(device_id==_gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferleft, device_id, _kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->H_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferright, device_id, _kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferleft, device_id, _kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->H_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->H_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        kernel_load_ghostlayerH <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                        cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);
                    }

                    _kpar_list_host[device_id]->t = (step+step1+0.5)*_dt;
                    gpuErrchk(cudaMemcpyAsync(_kpar_list[device_id], _kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(e_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], h_border_event_list[device_id], 0);
                    kernel_update_E_bulk <<< _griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (_kpar_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], e_time_event_list[device_id], 0);
                    kernel_update_E_border <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                    cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    
                    if(_gpus_count>1){
                        Barrier(counter_E, mutex_E, cv_E, _gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferright, device_id, _kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->E_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else if(device_id==_gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferleft, device_id, _kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->E_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferright, device_id, _kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferleft, device_id, _kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->E_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->E_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        kernel_load_ghostlayerE <<< _griddimghost_list[device_id], blockdim,0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                        cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    }
                }
                
                // FDTD_capture_t1_fields
                Barrier(counter_t1, mutex_t1, cv_t1, _gpus_count);
                // gpuErrchk(cudaMemcpy(Ext, _kpar_list_host[device_id]->Ex, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Eyt, _kpar_list_host[device_id]->Ey, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Ezt, _kpar_list_host[device_id]->Ez, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Hxt, _kpar_list_host[device_id]->Hx, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Hyt, _kpar_list_host[device_id]->Hy, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Hzt, _kpar_list_host[device_id]->Hz, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // for(int ii=0; ii<It; ii++){
                //     for(int jj=0; jj<Jt; jj++){
                //         for(int kk=0; kk<Kt; kk++){
                //             _Ex_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Ext[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                //             _Ey_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Eyt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                //             _Ez_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Ezt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                //             _Hx_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Hxt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                //             _Hy_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Hyt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                //             _Hz_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Hzt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                //         }
                //     }
                // }

                // capture t1 fields on GPUs
                double *dEx_t1, *dEy_t1, *dEz_t1, *dHx_t1, *dHy_t1, *dHz_t1;
                gpuErrchk(cudaMalloc((void **)&dEx_t1, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dEy_t1, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dEz_t1, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHx_t1, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHy_t1, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHz_t1, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                // gpuErrchk(cudaMemcpy(dEx_t1, Ext, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dEy_t1, Eyt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dEz_t1, Ezt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHx_t1, Hxt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHy_t1, Hyt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHz_t1, Hzt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                gpuErrchk(cudaMemcpy(dEx_t1, _kpar_list_host[device_id]->Ex, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dEy_t1, _kpar_list_host[device_id]->Ey, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dEz_t1, _kpar_list_host[device_id]->Ez, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHx_t1, _kpar_list_host[device_id]->Hx, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHy_t1, _kpar_list_host[device_id]->Hy, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHz_t1, _kpar_list_host[device_id]->Hz, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));

                // Perform another Tn steps to get a third time point at t2
                for(int step2=0; step2<Tn; step2++){
                    _kpar_list_host[device_id]->t = (step+step2+Tn)*_dt;

                    gpuErrchk(cudaMemcpyAsync(_kpar_list[device_id], _kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(h_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], e_border_event_list[device_id], 0);
                    kernel_update_H_bulk <<< _griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (_kpar_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], h_time_event_list[device_id], 0);
                    kernel_update_H_border <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                    cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);

                    if(_gpus_count>1){
                        Barrier(counter_H, mutex_H, cv_H, _gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferright, device_id, _kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->H_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else if(device_id==_gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferleft, device_id, _kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->H_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], h_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferright, device_id, _kpar_list_host[device_id+1]->H_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->H_bufferleft, device_id, _kpar_list_host[device_id-1]->H_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->H_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->H_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->H_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        kernel_load_ghostlayerH <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                        cudaEventRecord(h_border_event_list[device_id], sync_stream_list[device_id]);
                    }

                    _kpar_list_host[device_id]->t = (step+step2+Tn+0.5)*_dt;
                    gpuErrchk(cudaMemcpyAsync(_kpar_list[device_id], _kpar_list_host[device_id], sizeof(kernelpar), cudaMemcpyHostToDevice, compute_stream_list[device_id]));
                    cudaEventRecord(e_time_event_list[device_id], compute_stream_list[device_id]);

                    cudaStreamWaitEvent(compute_stream_list[device_id], h_border_event_list[device_id], 0);
                    kernel_update_E_bulk <<< _griddim_list[device_id], blockdim, 0, compute_stream_list[device_id] >>> (_kpar_list[device_id]);

                    cudaStreamWaitEvent(sync_stream_list[device_id], e_time_event_list[device_id], 0);
                    kernel_update_E_border <<< _griddimghost_list[device_id], blockdim, 0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                    cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    
                    if(_gpus_count>1){
                        Barrier(counter_E, mutex_E, cv_E, _gpus_count);
                        if(device_id==0){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferright, device_id, _kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->E_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else if(device_id==_gpus_count-1){
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferleft, device_id, _kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->E_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        else{
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id+1], 0);
                            cudaStreamWaitEvent(sync_stream_list[device_id], e_border_event_list[device_id-1], 0);
                            if(_P2Pworking){
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferright, device_id, _kpar_list_host[device_id+1]->E_ghostleft, device_id+1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyPeerAsync(_kpar_list_host[device_id]->E_bufferleft, device_id, _kpar_list_host[device_id-1]->E_ghostright, device_id-1, 3*_ghostsize_list[device_id]*sizeof(double), sync_stream_list[device_id]));
                            }else{
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id+1]->E_ghostleft, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferright, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(Emid_list[device_id], _kpar_list_host[device_id-1]->E_ghostright, 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyDeviceToHost, sync_stream_list[device_id]));
                                gpuErrchk(cudaMemcpyAsync(_kpar_list_host[device_id]->E_bufferleft, Emid_list[device_id], 3*_ghostsize_list[device_id]*sizeof(double), cudaMemcpyHostToDevice, sync_stream_list[device_id]));
                            }
                        }
                        kernel_load_ghostlayerE <<< _griddimghost_list[device_id], blockdim,0, sync_stream_list[device_id] >>> (_kpar_list[device_id]);
                        cudaEventRecord(e_border_event_list[device_id], sync_stream_list[device_id]);
                    }
                }

                // Capture fields at t2
                Barrier(counter_t2, mutex_t2, cv_t2, _gpus_count);
                // gpuErrchk(cudaMemcpy(Ext, _kpar_list_host[device_id]->Ex, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Eyt, _kpar_list_host[device_id]->Ey, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Ezt, _kpar_list_host[device_id]->Ez, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Hxt, _kpar_list_host[device_id]->Hx, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Hyt, _kpar_list_host[device_id]->Hy, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));
                // gpuErrchk(cudaMemcpy(Hzt, _kpar_list_host[device_id]->Hz, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToHost));

                // capture t2 fields on GPUs
                double *dEx_t2, *dEy_t2, *dEz_t2, *dHx_t2, *dHy_t2, *dHz_t2;
                gpuErrchk(cudaMalloc((void **)&dEx_t2, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dEy_t2, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dEz_t2, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHx_t2, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHy_t2, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                gpuErrchk(cudaMalloc((void **)&dHz_t2, (It+2)*(Jt+2)*(Kt+2)*sizeof(double)));
                // gpuErrchk(cudaMemcpy(dEx_t2, Ext, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dEy_t2, Eyt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dEz_t2, Ezt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHx_t2, Hxt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHy_t2, Hyt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                // gpuErrchk(cudaMemcpy(dHz_t2, Hzt, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyHostToDevice));
                gpuErrchk(cudaMemcpy(dEx_t2, _kpar_list_host[device_id]->Ex, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dEy_t2, _kpar_list_host[device_id]->Ey, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dEz_t2, _kpar_list_host[device_id]->Ez, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHx_t2, _kpar_list_host[device_id]->Hx, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHy_t2, _kpar_list_host[device_id]->Hy, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));
                gpuErrchk(cudaMemcpy(dHz_t2, _kpar_list_host[device_id]->Hz, (It+2)*(Jt+2)*(Kt+2)*sizeof(double), cudaMemcpyDeviceToDevice));


                // Calculate complex fields
                if(device_id==0){
                    auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                    std::cout << "Calc fields:" << ctime(&timenow);
                }
                t0 = step*_dt;
                t1 = (step+Tn)*_dt;
                t2 = (step+Tn+Tn)*_dt;
                double t0H, t1H, t2H;//, phi, A, f0, f1, f2;
                t0H = t0 - _dt/2; t1H = t1 - _dt/2; t2H = t2 - _dt/2;

                complex128 *Ex_out, *Ey_out, *Ez_out, *Hx_out, *Hy_out, *Hz_out,
                            *dEx_out, *dEy_out, *dEz_out, *dHx_out, *dHy_out, *dHz_out;
                Ex_out = (complex128 *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(complex128));
                Ey_out = (complex128 *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(complex128));
                Ez_out = (complex128 *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(complex128));
                Hx_out = (complex128 *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(complex128));
                Hy_out = (complex128 *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(complex128));
                Hz_out = (complex128 *)malloc((It+2)*(Jt+2)*(Kt+2)*sizeof(complex128));
                gpuErrchk(cudaMalloc((void **)&dEx_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128)));
                gpuErrchk(cudaMalloc((void **)&dEy_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128)));
                gpuErrchk(cudaMalloc((void **)&dEz_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128)));
                gpuErrchk(cudaMalloc((void **)&dHx_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128)));
                gpuErrchk(cudaMalloc((void **)&dHy_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128)));
                gpuErrchk(cudaMalloc((void **)&dHz_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128)));

                size_t numBlock = ceil((It+2)*(Jt+2)*(Kt+2)/(double)blockdim);
                kernel_calc_complex_fields <<< numBlock, blockdim >>> (t0, t1, t2, dEx_t0, dEx_t1, dEx_t2, (It+2)*(Jt+2)*(Kt+2), dEx_out);
                kernel_calc_complex_fields <<< numBlock, blockdim >>> (t0, t1, t2, dEy_t0, dEy_t1, dEy_t2, (It+2)*(Jt+2)*(Kt+2), dEy_out);
                kernel_calc_complex_fields <<< numBlock, blockdim >>> (t0, t1, t2, dEz_t0, dEz_t1, dEz_t2, (It+2)*(Jt+2)*(Kt+2), dEz_out);
                kernel_calc_complex_fields <<< numBlock, blockdim >>> (t0H, t1H, t2H, dHx_t0, dHx_t1, dHx_t2, (It+2)*(Jt+2)*(Kt+2), dHx_out);
                kernel_calc_complex_fields <<< numBlock, blockdim >>> (t0H, t1H, t2H, dHy_t0, dHy_t1, dHy_t2, (It+2)*(Jt+2)*(Kt+2), dHy_out);
                kernel_calc_complex_fields <<< numBlock, blockdim >>> (t0H, t1H, t2H, dHz_t0, dHz_t1, dHz_t2, (It+2)*(Jt+2)*(Kt+2), dHz_out);
                gpuErrchk(cudaMemcpy(Ex_out, dEx_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128), cudaMemcpyDeviceToHost));
                gpuErrchk(cudaMemcpy(Ey_out, dEy_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128), cudaMemcpyDeviceToHost));
                gpuErrchk(cudaMemcpy(Ez_out, dEz_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128), cudaMemcpyDeviceToHost));
                gpuErrchk(cudaMemcpy(Hx_out, dHx_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128), cudaMemcpyDeviceToHost));
                gpuErrchk(cudaMemcpy(Hy_out, dHy_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128), cudaMemcpyDeviceToHost));
                gpuErrchk(cudaMemcpy(Hz_out, dHz_out, (It+2)*(Jt+2)*(Kt+2)*sizeof(complex128), cudaMemcpyDeviceToHost));

                for(int ii=0; ii<It; ii++){
                    for(int jj=0; jj<Jt; jj++){
                        for(int kk=0; kk<Kt; kk++){
                            _Ex_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Ex_out[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            _Ey_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Ey_out[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            _Ez_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Ez_out[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            _Hx_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Hx_out[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            _Hy_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Hy_out[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            _Hz_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)] = Hz_out[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                        }
                    }
                } // Calculate complex fields ends

                delete[] Ex_out; delete[] Ey_out; delete[] Ez_out; delete[] Hx_out; delete[] Hy_out; delete[] Hz_out;
                gpuErrchk(cudaFree(dEx_out)); gpuErrchk(cudaFree(dEy_out)); gpuErrchk(cudaFree(dEz_out));
                gpuErrchk(cudaFree(dHx_out)); gpuErrchk(cudaFree(dHy_out)); gpuErrchk(cudaFree(dHz_out));
                gpuErrchk(cudaFree(dEx_t0)); gpuErrchk(cudaFree(dEy_t0)); gpuErrchk(cudaFree(dEz_t0));
                gpuErrchk(cudaFree(dHx_t0)); gpuErrchk(cudaFree(dHy_t0)); gpuErrchk(cudaFree(dHz_t0));
                gpuErrchk(cudaFree(dEx_t1)); gpuErrchk(cudaFree(dEy_t1)); gpuErrchk(cudaFree(dEz_t1));
                gpuErrchk(cudaFree(dHx_t1)); gpuErrchk(cudaFree(dHy_t1)); gpuErrchk(cudaFree(dHz_t1));
                gpuErrchk(cudaFree(dEx_t2)); gpuErrchk(cudaFree(dEy_t2)); gpuErrchk(cudaFree(dEz_t2));
                gpuErrchk(cudaFree(dHx_t2)); gpuErrchk(cudaFree(dHy_t2)); gpuErrchk(cudaFree(dHz_t2));
/*
                for(int ii=0; ii<It; ii++){
                    for(int jj=0; jj<Jt; jj++){
                        for(int kk=0; kk<Kt; kk++){
                            // Ex
                            f0 = _Ex_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f1 = _Ex_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f2 = Ext[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            phi = calc_phase(t0, t1, t2, f0, f1, f2);
                            A = calc_amplitude(t0, t1, t2, f0, f1, f2, phi);
                            if(A<0){ A*=-1; phi+=M_PI; }
                            _Ex_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real = A*cos(phi);
                            _Ex_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].imag = -A*sin(phi);

                            // Ey
                            f0 = _Ey_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f1 = _Ey_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f2 = Eyt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            phi = calc_phase(t0, t1, t2, f0, f1, f2);
                            A = calc_amplitude(t0, t1, t2, f0, f1, f2, phi);
                            if(A<0){ A*=-1; phi+=M_PI; }
                            _Ey_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real = A*cos(phi);
                            _Ey_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].imag = -A*sin(phi);

                            // Ez
                            f0 = _Ez_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f1 = _Ez_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f2 = Ezt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            phi = calc_phase(t0, t1, t2, f0, f1, f2);
                            A = calc_amplitude(t0, t1, t2, f0, f1, f2, phi);
                            if(A<0){ A*=-1; phi+=M_PI; }
                            _Ez_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real = A*cos(phi);
                            _Ez_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].imag = -A*sin(phi);

                            // Hx
                            f0 = _Hx_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f1 = _Hx_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f2 = Hxt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            phi = calc_phase(t0H, t1H, t2H, f0, f1, f2);
                            A = calc_amplitude(t0H, t1H, t2H, f0, f1, f2, phi);
                            if(A<0){ A*=-1; phi+=M_PI; }
                            _Hx_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real = A*cos(phi);
                            _Hx_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].imag = -A*sin(phi);

                            // Hy
                            f0 = _Hy_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f1 = _Hy_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f2 = Hyt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            phi = calc_phase(t0H, t1H, t2H, f0, f1, f2);
                            A = calc_amplitude(t0H, t1H, t2H, f0, f1, f2, phi);
                            if(A<0){ A*=-1; phi+=M_PI; }
                            _Hy_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real = A*cos(phi);
                            _Hy_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].imag = -A*sin(phi);

                            // Hz
                            f0 = _Hz_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f1 = _Hz_t1[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real;
                            f2 = Hzt[(ii+1)*(Jt+2)*(Kt+2)+(jj+1)*(Kt+2)+(kk+1)];
                            phi = calc_phase(t0H, t1H, t2H, f0, f1, f2);
                            A = calc_amplitude(t0H, t1H, t2H, f0, f1, f2, phi);
                            if(A<0){ A*=-1; phi+=M_PI; }
                            _Hz_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].real = A*cos(phi);
                            _Hz_t0[(ii+i0t)*J*K+(jj+j0t)*K+(kk+k0t)].imag = -A*sin(phi);
                        }
                    }
                } // Calculate complex fields ends
*/
                // Free up memory
                delete[] Ext; delete[] Eyt; delete[] Ezt; delete[] Hxt; delete[] Hyt; delete[] Hzt;
                delete[] Emid_list[device_id];
            }  // try
            catch(std::runtime_error &error){
                std::cerr << "FDTD error in thread " << device_id << ":" << error.what() << std::endl;
            }
        });
    }
    for(auto &thread: threads_kernels)
        thread.join();

    timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Finish FDTD:" << ctime(&timenow);

    // Free up memory
    delete[] Ex0; delete[] Ex1; delete[] Ex2; delete[] Ey0; delete[] Ey1; delete[] Ey2; delete[] Ez0; delete[] Ez1; delete[] Ez2;
    delete[] phi0; delete[] phi1; delete[] A0; delete[] A1;
}

void fdtd::FDTD::update_E(int n, double t)
{
    double odx = _R/_dx;

    int pml_xmin = _w_pml_x0, pml_xmax = _Nx-_w_pml_x1,
        pml_ymin = _w_pml_y0, pml_ymax = _Ny-_w_pml_y1,
        pml_zmin = _w_pml_z0, pml_zmax = _Nz-_w_pml_z1;

    _kpar_host->pml_xmin = pml_xmin;
    _kpar_host->pml_xmax = pml_xmax;
    _kpar_host->pml_ymin = pml_ymin;
    _kpar_host->pml_ymax = pml_ymax;
    _kpar_host->pml_zmin = pml_zmin;
    _kpar_host->pml_zmax = pml_zmax;
    _kpar_host->odx = odx;
    _kpar_host->t = t;
    _kpar_host->src_T = _src_T;
    _kpar_host->src_min = _src_min;
    _kpar_host->src_k = _src_k;

    gpuErrchk(cudaMemcpy(_kpar_device,_kpar_host, sizeof(kernelpar), cudaMemcpyHostToDevice));
    kernel_update_E <<< ceil(_I*_J*_K/128.0), 128 >>> (_kpar_device);
}

///////////////////////////////////////////////////////////////////////////
// PML Management
///////////////////////////////////////////////////////////////////////////


void fdtd::FDTD::set_pml_widths(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    _w_pml_x0 = xmin; _w_pml_x1 = xmax;
    _w_pml_y0 = ymin; _w_pml_y1 = ymax;
    _w_pml_z0 = zmin; _w_pml_z1 = zmax;
}

void fdtd::FDTD::set_pml_properties(double sigma, double alpha, double kappa, double pow)
{
    _sigma = sigma;
    _alpha = alpha;
    _kappa = kappa;
    _pow   = pow;

    compute_pml_params();
}

void fdtd::FDTD::build_pml()
{
    int N,
        xmin = _w_pml_x0, xmax = _Nx-_w_pml_x1,
        ymin = _w_pml_y0, ymax = _Ny-_w_pml_y1,
        zmin = _w_pml_z0, zmax = _Nz-_w_pml_z1;

    // touches xmin boudary
    if(_k0 < xmin) {
        N = _I * _J * (xmin - _k0);

        // Clean up old arrays and allocate new ones
        delete [] _pml_Eyx0; _pml_Eyx0 = NULL;
        delete [] _pml_Ezx0; _pml_Ezx0 = NULL;
        _pml_Eyx0 = new double[N];
        _pml_Ezx0 = new double[N];

        delete [] _pml_Hyx0; _pml_Hyx0 = NULL;
        delete [] _pml_Hzx0; _pml_Hzx0 = NULL;
        _pml_Hyx0 = new double[N];
        _pml_Hzx0 = new double[N];
    }

    // touches xmax boundary
    if(_k0 +_K > xmax) {
        N = _I * _J * (_k0  + _K - xmax);

        // Clean up old arrays and allocate new ones
        delete [] _pml_Eyx1; _pml_Eyx1 = NULL;
        delete [] _pml_Ezx1; _pml_Ezx1 = NULL;
        _pml_Eyx1 = new double[N];
        _pml_Ezx1 = new double[N];

        delete [] _pml_Hyx1; _pml_Hyx1 = NULL;
        delete [] _pml_Hzx1; _pml_Hzx1 = NULL;
        _pml_Hyx1 = new double[N];
        _pml_Hzx1 = new double[N];
    }

    // touches ymin boundary
    if(_j0 < ymin) {
        N = _I * _K * (ymin - _j0);

        delete [] _pml_Exy0; _pml_Exy0 = NULL;
        delete [] _pml_Ezy0; _pml_Ezy0 = NULL;
        _pml_Exy0 = new double[N];
        _pml_Ezy0 = new double[N];

        delete [] _pml_Hxy0; _pml_Hxy0 = NULL;
        delete [] _pml_Hzy0; _pml_Hzy0 = NULL;
        _pml_Hxy0 = new double[N];
        _pml_Hzy0 = new double[N];
    }

    // touches ymax boundary
    if(_j0 + _J > ymax) {
        N = _I * _K * (_j0 + _J - ymax);

        delete [] _pml_Exy1; _pml_Exy1 = NULL;
        delete [] _pml_Ezy1; _pml_Ezy1 = NULL;
        _pml_Exy1 = new double[N];
        _pml_Ezy1 = new double[N];

        delete [] _pml_Hxy1; _pml_Hxy1 = NULL;
        delete [] _pml_Hzy1; _pml_Hzy1 = NULL;
        _pml_Hxy1 = new double[N];
        _pml_Hzy1 = new double[N];
    }

    // touches zmin boundary
    if(_i0 < zmin) {
        N = _J * _K * (zmin - _i0);

        delete [] _pml_Exz0; _pml_Exz0 = NULL;
        delete [] _pml_Eyz0; _pml_Eyz0 = NULL;
        _pml_Exz0 = new double[N];
        _pml_Eyz0 = new double[N];

        delete [] _pml_Hxz0; _pml_Hxz0 = NULL;
        delete [] _pml_Hyz0; _pml_Hyz0 = NULL;
        _pml_Hxz0 = new double[N];
        _pml_Hyz0 = new double[N];
    }

    // touches zmax boundary
    if(_i0 + _I > zmax) {
        N = _J * _K * (_i0 + _I - zmax);

        delete [] _pml_Hxz1; _pml_Hxz1 = NULL;
        delete [] _pml_Hyz1; _pml_Hyz1 = NULL;
        _pml_Exz1 = new double[N];
        _pml_Eyz1 = new double[N];

        delete [] _pml_Hxz1; _pml_Hxz1 = NULL;
        delete [] _pml_Hyz1; _pml_Hyz1 = NULL;
        _pml_Hxz1 = new double[N];
        _pml_Hyz1 = new double[N];
    }

    // (re)compute the spatially-dependent PML parameters
    compute_pml_params();
}

void fdtd::FDTD::reset_pml()
{
    int N,
        xmin = _w_pml_x0, xmax = _Nx-_w_pml_x1,
        ymin = _w_pml_y0, ymax = _Ny-_w_pml_y1,
        zmin = _w_pml_z0, zmax = _Nz-_w_pml_z1;

    // touches xmin boudary
    if(_k0 < xmin) {
        N = _I * _J * (xmin - _k0);
        std::fill(_pml_Eyx0, _pml_Eyx0 + N, 0);
        std::fill(_pml_Ezx0, _pml_Ezx0 + N, 0);
        std::fill(_pml_Hyx0, _pml_Hyx0 + N, 0);
        std::fill(_pml_Hzx0, _pml_Hzx0 + N, 0);
    }

    // touches xmax boundary
    if(_k0 +_K > xmax) {
        N = _I * _J * (_k0  + _K - xmax);
        std::fill(_pml_Eyx1, _pml_Eyx1 + N, 0);
        std::fill(_pml_Ezx1, _pml_Ezx1 + N, 0);
        std::fill(_pml_Hyx1, _pml_Hyx1 + N, 0);
        std::fill(_pml_Hzx1, _pml_Hzx1 + N, 0);
    }

    // touches ymin boundary
    if(_j0 < ymin) {
        N = _I * _K * (ymin - _j0);
        std::fill(_pml_Exy0, _pml_Exy0 + N, 0);
        std::fill(_pml_Ezy0, _pml_Ezy0 + N, 0);
        std::fill(_pml_Hxy0, _pml_Hxy0 + N, 0);
        std::fill(_pml_Hzy0, _pml_Hzy0 + N, 0);
    }

    // touches ymax boundary
    if(_j0 + _J > ymax) {
        N = _I * _K * (_j0 + _J - ymax);
        std::fill(_pml_Exy1, _pml_Exy1 + N, 0);
        std::fill(_pml_Ezy1, _pml_Ezy1 + N, 0);
        std::fill(_pml_Hxy1, _pml_Hxy1 + N, 0);
        std::fill(_pml_Hzy1, _pml_Hzy1 + N, 0);
    }

    // touches zmin boundary
    if(_i0 < zmin) {
        N = _J * _K * (zmin - _i0);
        std::fill(_pml_Exz0, _pml_Exz0 + N, 0);
        std::fill(_pml_Eyz0, _pml_Eyz0 + N, 0);
        std::fill(_pml_Hxz0, _pml_Hxz0 + N, 0);
        std::fill(_pml_Hyz0, _pml_Hyz0 + N, 0);
    }

    // touches zmax boundary
    if(_i0 + _I > zmax) {
        N = _J * _K * (_i0 + _I - zmax);
        std::fill(_pml_Exz1, _pml_Exz1 + N, 0);
        std::fill(_pml_Eyz1, _pml_Eyz1 + N, 0);
        std::fill(_pml_Hxz1, _pml_Hxz1 + N, 0);
        std::fill(_pml_Hyz1, _pml_Hyz1 + N, 0);
    }

}

void fdtd::FDTD::compute_pml_params()
{
    double pml_dist, pml_factor, sigma, alpha, kappa, b, c;

    // clean up the previous arrays and allocate new ones
    delete [] _kappa_H_x; _kappa_H_x = new double[_w_pml_x0 + _w_pml_x1];
    delete [] _kappa_H_y; _kappa_H_y = new double[_w_pml_y0 + _w_pml_y1];
    delete [] _kappa_H_z; _kappa_H_z = new double[_w_pml_z0 + _w_pml_z1];

    delete [] _kappa_E_x; _kappa_E_x = new double[_w_pml_x0 + _w_pml_x1];
    delete [] _kappa_E_y; _kappa_E_y = new double[_w_pml_y0 + _w_pml_y1];
    delete [] _kappa_E_z; _kappa_E_z = new double[_w_pml_z0 + _w_pml_z1];

    delete [] _bHx; _bHx = new double[_w_pml_x0 + _w_pml_x1];
    delete [] _bHy; _bHy = new double[_w_pml_y0 + _w_pml_y1];
    delete [] _bHz; _bHz = new double[_w_pml_z0 + _w_pml_z1];

    delete [] _bEx; _bEx = new double[_w_pml_x0 + _w_pml_x1];
    delete [] _bEy; _bEy = new double[_w_pml_y0 + _w_pml_y1];
    delete [] _bEz; _bEz = new double[_w_pml_z0 + _w_pml_z1];

    delete [] _cHx; _cHx = new double[_w_pml_x0 + _w_pml_x1];
    delete [] _cHy; _cHy = new double[_w_pml_y0 + _w_pml_y1];
    delete [] _cHz; _cHz = new double[_w_pml_z0 + _w_pml_z1];

    delete [] _cEx; _cEx = new double[_w_pml_x0 + _w_pml_x1];
    delete [] _cEy; _cEy = new double[_w_pml_y0 + _w_pml_y1];
    delete [] _cEz; _cEz = new double[_w_pml_z0 + _w_pml_z1];

    // calculate the PML parameters. These parameters are all functions of
    // the distance from the ONSET of the PML edge (which begins in the simulation
    // domain interior.
    // Note: PML parameters are ordered such that distance from PML onset
    // always increases with index.
    
    // setup xmin PML parameters
    for(int k = 0; k < _w_pml_x0; k++) {
        pml_dist = double(k - 0.5)/_w_pml_x0; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);
        if(pml_factor < 0) pml_factor = 0;

        // compute H coefficients
        sigma = _sigma * pml_factor;
        alpha = _alpha * (1-pml_factor);
        kappa = (_kappa-1.0) * pml_factor+1.0;
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_H_x[k] = kappa;
        _bHx[k] = b;
        _cHx[k] = c;

        pml_dist = double(k)/_w_pml_x0; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        // compute E coefficients
        sigma = _sigma * pml_factor;
        alpha = _alpha * (1-pml_factor);
        kappa = (_kappa-1) * pml_factor+1;
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_E_x[k] = kappa;
        _bEx[k] = b;
        _cEx[k] = c;

    }
    for(int k = 0; k < _w_pml_x1; k++) {
        // compute H coefficients
        pml_dist = double(k + 0.5)/_w_pml_x1; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_H_x[_w_pml_x0 + k] = kappa;
        _bHx[_w_pml_x0 + k] = b;
        _cHx[_w_pml_x0 + k] = c;

        //compute E coefficients
        pml_dist = double(k)/_w_pml_x1; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_E_x[_w_pml_x0 + k] = kappa;
        _bEx[_w_pml_x0 + k] = b;
        _cEx[_w_pml_x0 + k] = c;
    }
    for(int j = 0; j < _w_pml_y0; j++) {
        // calc H coefficients
        pml_dist = double(j - 0.5)/_w_pml_y0; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);
        if(pml_factor < 0) pml_factor = 0;

        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_H_y[j] = kappa;
        _bHy[j] = b;
        _cHy[j] = c;

        // calc E coefficients
        pml_dist = double(j)/_w_pml_y0; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_E_y[j] = kappa;
        _bEy[j] = b;
        _cEy[j] = c;
    
    }
    for(int j = 0; j < _w_pml_y1; j++) {
         // calc H coeffs
         pml_dist = double(j + 0.5)/_w_pml_y1; // distance from pml edge
         pml_factor = pml_ramp(pml_dist);

         sigma = _sigma * pml_factor;
         kappa = (_kappa-1) * pml_factor+1;
         alpha = _alpha * (1-pml_factor);
         b = exp(-_dt*(sigma/kappa + alpha));
         if(b == 1) c = 0;
         else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_H_y[_w_pml_y0 + j] = kappa;
        _bHy[_w_pml_y0 + j] = b;
        _cHy[_w_pml_y0 + j] = c;

        // compute E coefficients
        pml_dist = double(j)/_w_pml_y1; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha); 

        _kappa_E_y[_w_pml_y0 + j] = kappa;
        _bEy[_w_pml_y0 + j] = b;
        _cEy[_w_pml_y0 + j] = c;
    }

    for(int i = 0; i < _w_pml_z0; i++) {
        // calc H coeffs
        pml_dist = double(i)/_w_pml_z0; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c= 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_H_z[i] = kappa;
        _bHz[i] = b;
        _cHz[i] = c;

        // calc E coeffs
        pml_dist = double(i+0.5)/_w_pml_z0; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        // compute coefficients
        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_E_z[i] = kappa;
        _bEz[i] = b;
        _cEz[i] = c;
    }

    for(int i = 0; i < _w_pml_z1; i++) {
        // calc H coeffs
        pml_dist = double(i)/_w_pml_z1; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);

        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_H_z[_w_pml_z0 + i] = kappa;
        _bHz[_w_pml_z0 + i] = b;
        _cHz[_w_pml_z0 + i] = c;

        // calc E coeffs
        pml_dist = double(i - 0.5)/_w_pml_z1; // distance from pml edge
        pml_factor = pml_ramp(pml_dist);
        if(pml_factor < 0) pml_factor = 0;

        // compute coefficients
        sigma = _sigma * pml_factor;
        kappa = (_kappa-1) * pml_factor+1;
        alpha = _alpha * (1-pml_factor);
        b = exp(-_dt*(sigma/kappa + alpha));
        if(b == 1) c = 0;
        else c = (b - 1)*sigma / (sigma*kappa + kappa*kappa*alpha);

        _kappa_E_z[_w_pml_z0 + i] = kappa;
        _bEz[_w_pml_z0 + i] = b;
        _cEz[_w_pml_z0 + i] = c;
    }
}

double fdtd::FDTD::pml_ramp(double pml_dist)
{
    return std::pow(pml_dist, _pow);
}

///////////////////////////////////////////////////////////////////////////
// Amp/Phase Calculation management Management
///////////////////////////////////////////////////////////////////////////
void fdtd::FDTD::set_t0_arrays(complex128 *Ex_t0, complex128 *Ey_t0, complex128 *Ez_t0,
                                complex128 *Hx_t0, complex128 *Hy_t0, complex128 *Hz_t0)
{
    _Ex_t0 = Ex_t0; _Ey_t0 = Ey_t0; _Ez_t0 = Ez_t0;
    _Hx_t0 = Hx_t0; _Hy_t0 = Hy_t0; _Hz_t0 = Hz_t0;
}

void fdtd::FDTD::set_t1_arrays(complex128 *Ex_t1, complex128 *Ey_t1, complex128 *Ez_t1,
complex128 *Hx_t1, complex128 *Hy_t1, complex128 *Hz_t1)
{
_Ex_t1 = Ex_t1; _Ey_t1 = Ey_t1; _Ez_t1 = Ez_t1;
_Hx_t1 = Hx_t1; _Hy_t1 = Hy_t1; _Hz_t1 = Hz_t1;
}

void fdtd::FDTD::capture_pbox_fields(std::complex<double> *Ex_full, std::complex<double> *Ey_full,
                                    std::complex<double> *Ez_full, std::complex<double> *Ex_pbox,
                                    std::complex<double> *Ey_pbox, std::complex<double> *Ez_pbox)
{
    int ind_full, ind_pbox;

    for(int i = _i1; i < _i2; i++) {
        for(int j = 0; j < _J; j++) {
            for(int k = 0; k < _K; k++) {
                ind_full = i*_J*_K + j*_K + k;
                ind_pbox = (i-_i1)*_J*_K + j*_K + k;

                Ex_pbox[ind_pbox] = Ex_full[ind_full];
                Ey_pbox[ind_pbox] = Ey_full[ind_full];
                Ez_pbox[ind_pbox] = Ez_full[ind_full];
            }
        }
    }
}

void fdtd::FDTD::capture_t0_fields()
{
    int ind_local, ind_global;

    for(int i = 0; i < _I; i++) {
        for(int j = 0; j < _J; j++) {
            for(int k = 0; k < _K; k++) {
                ind_local = (i+1)*(_J+2)*(_K+2) + (j+1)*(_K+2) + k + 1;
                ind_global = i*_J*_K + j*_K + k;

                // Copy the fields at the current time to the auxillary arrays
                _Ex_t0[ind_global] = _Ex[ind_local];
                _Ey_t0[ind_global] = _Ey[ind_local];
                _Ez_t0[ind_global] = _Ez[ind_local];

                _Hx_t0[ind_global] = _Hx[ind_local];
                _Hy_t0[ind_global] = _Hy[ind_local];
                _Hz_t0[ind_global] = _Hz[ind_local];
            }
        }
    }

}

void fdtd::FDTD::capture_t1_fields()
{
    int ind_local, ind_global;

    for(int i = 0; i < _I; i++) {
        for(int j = 0; j < _J; j++) {
            for(int k = 0; k < _K; k++) {
                ind_local = (i+1)*(_J+2)*(_K+2) + (j+1)*(_K+2) + k + 1;
                ind_global = i*_J*_K + j*_K + k;

                // Copy the fields at the current time to the auxillary arrays
                _Ex_t1[ind_global] = _Ex[ind_local];
                _Ey_t1[ind_global] = _Ey[ind_local];
                _Ez_t1[ind_global] = _Ez[ind_local];

                _Hx_t1[ind_global] = _Hx[ind_local];
                _Hy_t1[ind_global] = _Hy[ind_local];
                _Hz_t1[ind_global] = _Hz[ind_local];
            }
        }
    }

}

void fdtd::FDTD::calc_complex_fields(double t0, double t1)
{
    double f0, f1, phi, A, t0H, t1H;
    int ind_local, ind_global;

    t0H = t0 - 0.5*_dt;
    t1H = t1 - 0.5*_dt;

    for(int i = 0; i < _I; i++) {
        for(int j = 0; j < _J; j++) {
            for(int k = 0; k < _K; k++) {
                ind_local = (i+1)*(_J+2)*(_K+2) + (j+1)*(_K+2) + k + 1;
                ind_global = i*_J*_K + j*_K + k;
                
                // Compute amplitude and phase for Ex
                // Note: we are careful to assume exp(-i*w*t) time dependence
                f0 = _Ex_t0[ind_global].real;
                f1 = _Ex[ind_local];
                phi = calc_phase(t0, t1, f0, f1);
                A = calc_amplitude(t0, t1, f0, f1, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Ex_t0[ind_global].real = A*cos(phi);
                _Ex_t0[ind_global].imag = -A*sin(phi); 

                // Ey
                f0 = _Ey_t0[ind_global].real;
                f1 = _Ey[ind_local];
                phi = calc_phase(t0, t1, f0, f1);
                A = calc_amplitude(t0, t1, f0, f1, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Ey_t0[ind_global].real = A*cos(phi);
                _Ey_t0[ind_global].imag = -A*sin(phi); 

                // Ez
                f0 = _Ez_t0[ind_global].real;
                f1 = _Ez[ind_local];
                phi = calc_phase(t0, t1, f0, f1);
                A = calc_amplitude(t0, t1, f0, f1, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Ez_t0[ind_global].real = A*cos(phi);
                _Ez_t0[ind_global].imag = -A*sin(phi); 

                // Hx
                f0 = _Hx_t0[ind_global].real;
                f1 = _Hx[ind_local];
                phi = calc_phase(t0H, t1H, f0, f1);
                A = calc_amplitude(t0H, t1H, f0, f1, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Hx_t0[ind_global].real = A*cos(phi);
                _Hx_t0[ind_global].imag = -A*sin(phi); 

                // Hy
                f0 = _Hy_t0[ind_global].real;
                f1 = _Hy[ind_local];
                phi = calc_phase(t0H, t1H, f0, f1);
                A = calc_amplitude(t0H, t1H, f0, f1, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Hy_t0[ind_global].real = A*cos(phi);
                _Hy_t0[ind_global].imag = -A*sin(phi); 

                // Hz
                f0 = _Hz_t0[ind_global].real;
                f1 = _Hz[ind_local];
                phi = calc_phase(t0H, t1H, f0, f1);
                A = calc_amplitude(t0H, t1H, f0, f1, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Hz_t0[ind_global].real = A*cos(phi);
                _Hz_t0[ind_global].imag = -A*sin(phi); 
            }
        }
    }

}


void fdtd::FDTD::calc_complex_fields(double t0, double t1, double t2)
{
    double f0, f1, f2, phi, A, t0H, t1H, t2H;
    int ind_local, ind_global;

    t0H = t0 - 0.5*_dt;
    t1H = t1 - 0.5*_dt;
    t2H = t2 - 0.5*_dt;

    for(int i = 0; i < _I; i++) {
        for(int j = 0; j < _J; j++) {
            for(int k = 0; k < _K; k++) {
                ind_local = (i+1)*(_J+2)*(_K+2) + (j+1)*(_K+2) + k + 1;
                ind_global = i*_J*_K + j*_K + k;

                // Compute amplitude and phase for Ex
                // Note: we are careful to assume exp(-i*w*t) time dependence
                f0 = _Ex_t0[ind_global].real;
                f1 = _Ex_t1[ind_global].real;
                f2 = _Ex[ind_local];
                phi = calc_phase(t0, t1, t2, f0, f1, f2);
                A = calc_amplitude(t0, t1, t2, f0, f1, f2, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Ex_t0[ind_global].real = A*cos(phi);
                _Ex_t0[ind_global].imag = -A*sin(phi); 

                // Ey
                f0 = _Ey_t0[ind_global].real;
                f1 = _Ey_t1[ind_global].real;
                f2 = _Ey[ind_local];
                phi = calc_phase(t0, t1, t2, f0, f1, f2);
                A = calc_amplitude(t0, t1, t2, f0, f1, f2, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Ey_t0[ind_global].real = A*cos(phi);
                _Ey_t0[ind_global].imag = -A*sin(phi); 

                // Ez
                f0 = _Ez_t0[ind_global].real;
                f1 = _Ez_t1[ind_global].real;
                f2 = _Ez[ind_local];
                phi = calc_phase(t0, t1, t2, f0, f1, f2);
                A = calc_amplitude(t0, t1, t2, f0, f1, f2, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Ez_t0[ind_global].real = A*cos(phi);
                _Ez_t0[ind_global].imag = -A*sin(phi); 

                // Hx
                f0 = _Hx_t0[ind_global].real;
                f1 = _Hx_t1[ind_global].real;
                f2 = _Hx[ind_local];
                phi = calc_phase(t0H, t1H, t2H, f0, f1, f2);
                A = calc_amplitude(t0H, t1H, t2H, f0, f1, f2, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Hx_t0[ind_global].real = A*cos(phi);
                _Hx_t0[ind_global].imag = -A*sin(phi); 

                // Hy
                f0 = _Hy_t0[ind_global].real;
                f1 = _Hy_t1[ind_global].real;
                f2 = _Hy[ind_local];
                phi = calc_phase(t0H, t1H, t2H, f0, f1, f2);
                A = calc_amplitude(t0H, t1H, t2H, f0, f1, f2, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Hy_t0[ind_global].real = A*cos(phi);
                _Hy_t0[ind_global].imag = -A*sin(phi); 

                // Hz
                f0 = _Hz_t0[ind_global].real;
                f1 = _Hz_t1[ind_global].real;
                f2 = _Hz[ind_local];
                phi = calc_phase(t0H, t1H, t2H, f0, f1, f2);
                A = calc_amplitude(t0H, t1H, t2H, f0, f1, f2, phi);
                if(A < 0) {
                    A *= -1;
                    phi += M_PI;
                }
                _Hz_t0[ind_global].real = A*cos(phi);
                _Hz_t0[ind_global].imag = -A*sin(phi);

            }
        }
    }
}

inline double fdtd::calc_phase(double t0, double t1, double f0, double f1)
{
    if(f0 == 0.0 and f1 == 0) {
        return 0.0;
    }
    else {
        return atan((f1*sin(t0)-f0*sin(t1))/(f0*cos(t1)-f1*cos(t0)));
    }
}

inline double fdtd::calc_amplitude(double t0, double t1, double f0, double f1, double phase)
{
    if(f0*f0 > f1*f1) {
        return f1 / (sin(t1)*cos(phase) + cos(t1)*sin(phase));
    }
    else {
        return f0 / (sin(t0)*cos(phase) + cos(t0)*sin(phase));
    }
}

inline double fdtd::calc_phase(double t0, double t1, double t2, double f0, double f1, double f2)
{
    double f10 = f1 - f0,
           f21 = f2 - f1;

    if(f10 == 0 && f21 == 0) {
        return 0.0;
    }
    else {
        return atan2(f10*(sin(t2)-sin(t1)) - f21*(sin(t1)-sin(t0)), 
                     f21*(cos(t1)-cos(t0)) - f10*(cos(t2)-cos(t1)));
    }
}

inline double fdtd::calc_amplitude(double t0, double t1, double t2, double f0, double f1, double f2, double phase)
{
    double f21 = f2 - f1,
           f10 = f1 - f0;

    if(f21 == 0 && f10 == 0) {
        return 0.0;
    }
    else if(f21*f21 >= f10*f10) {
        return f21 / (cos(phase)*(sin(t2)-sin(t1)) + sin(phase)*(cos(t2)-cos(t1)));
    }
    else {
        return f10 / (cos(phase)*(sin(t1)-sin(t0)) + sin(phase)*(cos(t1)-cos(t0)));
    }
}

///////////////////////////////////////////////////////////////////////////
// Source management
///////////////////////////////////////////////////////////////////////////
void fdtd::FDTD::add_source(complex128 *Jx, complex128 *Jy, complex128 *Jz,
                            complex128 *Mx, complex128 *My, complex128 *Mz,
                            int i0, int j0, int k0, int I, int J, int K,
                            bool calc_phase)
{
    int ind=0;
    double real, imag;
    SourceArray src = {Jx, Jy, Jz, Mx, My, Mz, i0, j0, k0, I, J, K};

    // these source arrays may *actually* be complex-valued. In the time
    // domain, complex values correspond to temporal phase shifts. We need
    // to convert the complex value to an amplitude and phase. Fortunately,
    // we can use the memory that is already allocated for these values.
    // Specifically, we use src_array.real = amplitude and
    // src_array.imag = phase
    //
    // Important note: EMopt assumes the time dependence is exp(-i*omega*t).
    // In order to account for this minus sign, we need to invert the sign
    // of the calculated phase.
    if(calc_phase) {

    for(int i = 0; i < I; i++) {
        for(int j = 0; j < J; j++) {
            for(int k = 0; k < K; k++) {
                ind = i*J*K + j*K + k;

                
                // Jx
                real = Jx[ind].real;
                imag = Jx[ind].imag;

                Jx[ind].real = sqrt(real*real + imag*imag);
                if(imag == 0 && real == 0) Jx[ind].imag = 0.0;
                else Jx[ind].imag = -1*atan2(imag, real);

                // Jy
                real = Jy[ind].real;
                imag = Jy[ind].imag;

                Jy[ind].real = sqrt(real*real + imag*imag);
                if(imag == 0 && real == 0) Jy[ind].imag = 0.0;
                else Jy[ind].imag = -1*atan2(imag, real);

                // Jz
                real = Jz[ind].real;
                imag = Jz[ind].imag;

                Jz[ind].real = sqrt(real*real + imag*imag);
                if(imag == 0 && real == 0) Jz[ind].imag = 0.0;
                else Jz[ind].imag = -1*atan2(imag, real);

                // Mx
                real = Mx[ind].real;
                imag = Mx[ind].imag;

                Mx[ind].real = sqrt(real*real + imag*imag);
                if(imag == 0 && real == 0) Mx[ind].imag = 0.0;
                else Mx[ind].imag = -1*atan2(imag, real);

                // My
                real = My[ind].real;
                imag = My[ind].imag;

                My[ind].real = sqrt(real*real + imag*imag);
                if(imag == 0 && real == 0) My[ind].imag = 0.0;
                else My[ind].imag = -1*atan2(imag, real);

                // Mz
                real = Mz[ind].real;
                imag = Mz[ind].imag;

                Mz[ind].real = sqrt(real*real + imag*imag);
                if(imag == 0 && real == 0) Mz[ind].imag = 0.0;
                else Mz[ind].imag = -1*atan2(imag, real);
                
            }
        }
    }
    }

    _sources.push_back(src);

    _srclen += 1;
    _kpar_host->srclen += 1;
}

void fdtd::FDTD::clear_sources()
{
    _sources.clear();

    _srclen = 0;
    _kpar_host->srclen = 0;
}

void fdtd::FDTD::set_source_properties(double src_T, double src_min)
{
    _src_T = src_T;
    _src_min = src_min;
    _src_k = src_T*src_T / log((1+src_min)/src_min);
    //_src_k = 6.0 / src_T; // rate of src turn on
    //_src_n0 = 1.0 / _src_k * log((1.0-src_min)/src_min); // src delay
}

inline double fdtd::FDTD::src_func_t(int n, double t, double phase)
{
    //return sin(t + phase) / (1.0 + exp(-_src_k*(n-_src_n0)));
    if(t <= _src_T)
        return sin(t + phase)*((1+_src_min) * exp(-(t-_src_T)*(t-_src_T) / _src_k) - _src_min);
    else
        return sin(t + phase);
}

///////////////////////////////////////////////////////////////////////////
// Boundary Conditions
///////////////////////////////////////////////////////////////////////////
void fdtd::FDTD::set_bc(char* newbc)
{
    for(int i = 0; i < 3; i++){
        _bc[i] = newbc[i];
    }
}

///////////////////////////////////////////////////////////////////////////
// ctypes interface
///////////////////////////////////////////////////////////////////////////

fdtd::FDTD* FDTD_new()
{
    return new fdtd::FDTD();
}

void FDTD_set_wavelength(fdtd::FDTD* fdtd, double wavelength)
{
    fdtd->set_wavelength(wavelength);
}

void FDTD_set_physical_dims(fdtd::FDTD* fdtd, 
                            double X, double Y, double Z,
                            double dx, double dy, double dz)
{
    fdtd->set_physical_dims(X, Y, Z, dx, dy, dz);
}

void FDTD_set_grid_dims(fdtd::FDTD* fdtd, int Nx, int Ny, int Nz)
{
    fdtd->set_grid_dims(Nx, Ny, Nz);
}

void FDTD_set_local_grid(fdtd::FDTD* fdtd, 
                         int k0, int j0, int i0,
                         int K, int J, int I)
{
    fdtd->set_local_grid(k0, j0, i0, K, J, I);
}

void FDTD_set_local_grid_perturb(fdtd::FDTD* fdtd,
                         int i1, int i2)
{
    fdtd->set_local_grid_perturb(i1, i2);
}

void FDTD_set_dt(fdtd::FDTD* fdtd, double dt)
{
    fdtd->set_dt(dt);
}

void FDTD_set_rtol(fdtd::FDTD* fdtd, double rtol)
{
    fdtd->set_rtol(rtol);
}

void FDTD_set_Ncycle(fdtd::FDTD* fdtd, double Ncycle)
{
    fdtd->set_Ncycle(Ncycle);
}

void FDTD_set_complex_eps(fdtd::FDTD* fdtd, bool complex_eps)
{
    fdtd->set_complex_eps(complex_eps);
}

void FDTD_set_gpus_count(fdtd::FDTD* fdtd, int gpus_count)
{
    fdtd->set_gpus_count(gpus_count);
}

void FDTD_set_GPUDirect(fdtd::FDTD* fdtd)
{
    fdtd->set_GPUDirect();
}

void FDTD_set_domain_decomposition(fdtd::FDTD* fdtd, char domain_decomp)
{
    fdtd->set_domain_decomposition(domain_decomp);
}

void FDTD_copyCUDA_field_arrays(fdtd::FDTD* fdtd)
{
    fdtd->copyCUDA_field_arrays();
}

void FDTD_block_CUDA_free(fdtd::FDTD* fdtd)
{
    fdtd->block_CUDA_free();
}

void FDTD_block_CUDA_src_free(fdtd::FDTD* fdtd)
{
    fdtd->block_CUDA_src_free();
}

void FDTD_block_CUDA_malloc_memcpy(fdtd::FDTD* fdtd)
{
    fdtd->block_CUDA_malloc_memcpy();
}

void FDTD_block_CUDA_src_malloc_memcpy(fdtd::FDTD* fdtd)
{
    fdtd->block_CUDA_src_malloc_memcpy();
}

void FDTD_block_CUDA_multigpu_init(fdtd::FDTD* fdtd)
{
    fdtd->block_CUDA_multigpu_init();
}

void FDTD_block_CUDA_multigpu_free(fdtd::FDTD* fdtd)
{
    fdtd->block_CUDA_multigpu_free();
}

void FDTD_solve(fdtd::FDTD* fdtd)
{
    fdtd->solve();
}

void FDTD_calc_ydAx(fdtd::FDTD* fdtd, size_t size, size_t Nx, size_t Ny, size_t Nz, size_t i0, size_t i1, size_t i2,
                std::complex<double> *ydAx,
                std::complex<double> *Ex_adj, std::complex<double> *Ey_adj, std::complex<double> *Ez_adj,
                std::complex<double> *Ex_fwd, std::complex<double> *Ey_fwd, std::complex<double> *Ez_fwd,
                std::complex<double> *epsx0, std::complex<double> *epsy0, std::complex<double> *epsz0,
                std::complex<double> *epsxp, std::complex<double> *epsyp, std::complex<double> *epszp)
{
    fdtd->calc_ydAx(size, Nx, Ny, Nz, i0, i1, i2, ydAx, Ex_adj, Ey_adj, Ez_adj, Ex_fwd, Ey_fwd, Ez_fwd,
                epsx0, epsy0, epsz0, epsxp, epsyp, epszp);
}

void FDTD_set_field_arrays(fdtd::FDTD* fdtd,
                           double *Ex, double *Ey, double *Ez,
                           double *Hx, double *Hy, double *Hz)
{
    fdtd->set_field_arrays(Ex, Ey, Ez, Hx, Hy, Hz);
}

void FDTD_set_mat_arrays(fdtd::FDTD* fdtd,
                         complex128 *eps_x, complex128 *eps_y, complex128 *eps_z
                         )
{
    fdtd->set_mat_arrays(eps_x, eps_y, eps_z);
}

void FDTD_update_H(fdtd::FDTD* fdtd, int n, double t)
{
    fdtd->update_H(n, t);
}

void FDTD_update_E(fdtd::FDTD* fdtd, int n, double t)
{
    fdtd->update_E(n, t);
}

void FDTD_set_pml_widths(fdtd::FDTD* fdtd, int xmin, int xmax,
                                           int ymin, int ymax,
                                           int zmin, int zmax)
{
    fdtd->set_pml_widths(xmin, xmax, ymin, ymax, zmin, zmax);
}

void FDTD_set_pml_properties(fdtd::FDTD* fdtd, double sigma, double alpha,
                                               double kappa, double pow)
{
    fdtd->set_pml_properties(sigma, alpha, kappa, pow);
}

void FDTD_build_pml(fdtd::FDTD* fdtd)
{
    fdtd->build_pml();
}

void FDTD_reset_pml(fdtd::FDTD* fdtd)
{
    fdtd->reset_pml();
}

void FDTD_set_t0_arrays(fdtd::FDTD* fdtd,
                         complex128 *Ex_t0, complex128 *Ey_t0, complex128 *Ez_t0,
                         complex128 *Hx_t0, complex128 *Hy_t0, complex128 *Hz_t0)
{
    fdtd->set_t0_arrays(Ex_t0, Ey_t0, Ez_t0, Hx_t0, Hy_t0, Hz_t0);
}

void FDTD_set_t1_arrays(fdtd::FDTD* fdtd,
                         complex128 *Ex_t1, complex128 *Ey_t1, complex128 *Ez_t1,
                         complex128 *Hx_t1, complex128 *Hy_t1, complex128 *Hz_t1)
{
    fdtd->set_t1_arrays(Ex_t1, Ey_t1, Ez_t1, Hx_t1, Hy_t1, Hz_t1);
}

double FDTD_calc_phase_2T(double t0, double t1, double f0, double f1)
{
    return fdtd::calc_phase(t0, t1, f0, f1);
}

double FDTD_calc_amplitude_2T(double t0, double t1, double f0, double f1, double phase)
{
    return fdtd::calc_amplitude(t0, t1, f0, f1, phase);
}

double FDTD_calc_phase_3T(double t0, double t1, double t2, double f0, double f1, double f2)
{
    return fdtd::calc_phase(t0, t1, t2, f0, f1, f2);
}

double FDTD_calc_amplitude_3T(double t0, double t1, double t2, double f0, double f1, double f2, double phase)
{
    return fdtd::calc_amplitude(t0, t1, t2, f0, f1, f2, phase);
}

void FDTD_capture_pbox_fields(fdtd::FDTD* fdtd, std::complex<double> *Ex_full, std::complex<double> *Ey_full,
                            std::complex<double> *Ez_full, std::complex<double> *Ex_pbox,
                            std::complex<double> *Ey_pbox, std::complex<double> *Ez_pbox)
{
    fdtd->capture_pbox_fields(Ex_full, Ey_full, Ez_full, Ex_pbox, Ey_pbox, Ez_pbox);
}

void FDTD_capture_t0_fields(fdtd::FDTD* fdtd)
{
    fdtd->capture_t0_fields();
}

void FDTD_capture_t1_fields(fdtd::FDTD* fdtd)
{
    fdtd->capture_t1_fields();
}


void FDTD_calc_complex_fields_2T(fdtd::FDTD* fdtd, double t0, double t1)
{
    fdtd->calc_complex_fields(t0, t1);
}

void FDTD_calc_complex_fields_3T(fdtd::FDTD* fdtd, double t0, double t1, double t2)
{
    fdtd->calc_complex_fields(t0, t1, t2);
}

void FDTD_add_source(fdtd::FDTD* fdtd,
                     complex128 *Jx, complex128 *Jy, complex128 *Jz,
                     complex128 *Mx, complex128 *My, complex128 *Mz,
                     int i0, int j0, int k0, int I, int J, int K, bool calc_phase)
{
    fdtd->add_source(Jx, Jy, Jz, Mx, My, Mz, i0, j0, k0, I, J, K, calc_phase);
}

void FDTD_clear_sources(fdtd::FDTD* fdtd)
{
    fdtd->clear_sources();
}

void FDTD_set_source_properties(fdtd::FDTD* fdtd, double src_T, double src_min)
{
    fdtd->set_source_properties(src_T, src_min);
}

double FDTD_src_func_t(fdtd::FDTD* fdtd, int n, double t, double phase)
{
    return fdtd->src_func_t(n, t, phase);
}

void FDTD_set_bc(fdtd::FDTD* fdtd, char* newbc)
{
    fdtd->set_bc(newbc);
}

// Ghost communication helper functions
void FDTD_copy_to_ghost_comm(double* src, complex128* ghost, int I, int J, int K)
{
    unsigned int nstart = 0,
                 ind_ijk, ind_ghost;

    // copy xmin
    for(int i = 0; i < I; i++) {
        for(int j = 0; j < J; j++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + (j+1)*(K+2) + 1;
            ind_ghost = nstart + j + i*J;

            ghost[ind_ghost] = src[ind_ijk];
        }
    }

    // copy xmax
    nstart = I*J;
    for(int i = 0; i < I; i++) {
        for(int j = 0; j < J; j++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + (j+1)*(K+2) + K;
            ind_ghost = nstart + j + i*J;

            ghost[ind_ghost] = src[ind_ijk];
        }
    }

    // copy ymin
    nstart = 2*I*J;
    for(int i = 0; i < I; i++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + 1*(K+2) + k + 1;
            ind_ghost = nstart + k + i*K;

            ghost[ind_ghost] = src[ind_ijk];
        }
    }

    // copy ymax
    nstart = 2*I*J + I*K;
    for(int i = 0; i < I; i++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + J*(K+2) + k + 1;
            ind_ghost = nstart + k + i*K;

            ghost[ind_ghost] = src[ind_ijk];
        }
    }

    // copy zmin
    nstart = 2*I*J + 2*I*K;
    for(int j = 0; j < J; j++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = 1*(J+2)*(K+2) + (j+1)*(K+2) + k + 1;
            ind_ghost = nstart + k + j*K;

            ghost[ind_ghost] = src[ind_ijk];
        }
    }

    // copy zmax
    nstart = 2*I*J + 2*I*K + J*K;
    for(int j = 0; j < J; j++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = I*(J+2)*(K+2) + (j+1)*(K+2) + k + 1;
            ind_ghost = nstart + k + j*K;

            ghost[ind_ghost] = src[ind_ijk];
        }
    }
}

void FDTD_copy_from_ghost_comm(double* dest, complex128* ghost, int I, int J, int K)
{
    unsigned int nstart = 2*I*J + 2*I*K + 2*J*K,
                 ind_ijk, ind_ghost;

    // copy xmin
    for(int i = 0; i < I; i++) {
        for(int j = 0; j < J; j++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + (j+1)*(K+2) + 0;
            ind_ghost = nstart + j + i*J;

            dest[ind_ijk] = ghost[ind_ghost].real;
        }
    }

    // copy xmax
    nstart = 2*I*J + 2*I*K + 2*J*K + I*J;
    for(int i = 0; i < I; i++) {
        for(int j = 0; j < J; j++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + (j+1)*(K+2) + K+1;
            ind_ghost = nstart + j + i*J;

            dest[ind_ijk] = ghost[ind_ghost].real;
        }
    }

    // copy ymin
    nstart = 2*I*J + 2*I*K + 2*J*K + 2*I*J;
    for(int i = 0; i < I; i++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + 0*(K+2) + k + 1;
            ind_ghost = nstart + k + i*K;

            dest[ind_ijk] = ghost[ind_ghost].real;
        }
    }

    // copy ymax
    nstart = 2*I*J + 2*I*K + 2*J*K + 2*I*J + I*K;
    for(int i = 0; i < I; i++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = (i+1)*(J+2)*(K+2) + (J+1)*(K+2) + k + 1;
            ind_ghost = nstart + k + i*K;

            dest[ind_ijk] = ghost[ind_ghost].real;
        }
    }

    // copy zmin
    nstart = 2*I*J + 2*I*K + 2*J*K + 2*I*J + 2*I*K;
    for(int j = 0; j < J; j++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = 0*(J+2)*(K+2) + (j+1)*(K+2) + k + 1;
            ind_ghost = nstart + k + j*K;

            dest[ind_ijk] = ghost[ind_ghost].real;
        }
    }

    // copy zmax
    nstart = 2*I*J + 2*I*K + 2*J*K + 2*I*J + 2*I*K + J*K;
    for(int j = 0; j < J; j++) {
        for(int k = 0; k < K; k++) {
            ind_ijk = (I+1)*(J+2)*(K+2) + (j+1)*(K+2) + k + 1;
            ind_ghost = nstart + k + j*K;

            dest[ind_ijk] = ghost[ind_ghost].real;
        }
    }
}
