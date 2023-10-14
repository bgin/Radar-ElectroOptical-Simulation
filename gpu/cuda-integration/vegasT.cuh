#ifndef VEGAS_VEGAS_T_CUH
#define VEGAS_VEGAS_T_CUH

/*

code works for gaussian and sin using switch statement. device pointerr/template
slow down the code by 2x

chunksize needs to be tuned based on the ncalls. For now hardwired using a
switch statement


nvcc -O2 -DCUSTOM -o vegas vegasT.cu -arch=sm_70
OR
nvcc -O2 -DCURAND -o vegas vegasT.cu -arch=sm_70

example run command

nvprof ./vegas 0 6 0.0  10.0  1.0E+09  10 0 0 0

nvprof  ./vegas 1 9 -1.0  1.0  1.0E+07 15 10 10

nvprof ./vegas 2 2 -1.0 1.0  1.0E+09 1 0 0

Last three arguments are: total iterations, iteration

*/

//__device__ long idum = -1;

#define OUTFILEVAR 0

#include "common/cuda/cudaMemoryUtil.h"
#include "cuda/pagani/quad/quad.h"
#include "common/cuda/Volume.cuh"
#include "common/cuda/cudaApply.cuh"
#include "common/cuda/cudaArray.cuh"
#include "cuda/mcubes/seqCodesDefs.hh"
#include "cuda/mcubes/util/vegas_utils.cuh"
#include "cuda/mcubes/util/verbose_utils.cuh"
#include <assert.h>
#include <chrono>
#include <ctime>
#include <curand_kernel.h>
#include <fstream>
#include <inttypes.h>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cuda_profiler_api.h>

#include "common/integration_result.hh"

#define WARP_SIZE 32
#define BLOCK_DIM_X 128

// Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError()                                                       \
  {                                                                            \
    cudaError_t e = cudaGetLastError();                                        \
    if (e != cudaSuccess) {                                                    \
      printf("Cuda failure %s:%d: '%s'\n",                                     \
             __FILE__,                                                         \
             __LINE__,                                                         \
             cudaGetErrorString(e));                                           \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

namespace cuda_mcubes {

#define NR_END 1
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

  __inline__ __device__ double
  warpReduceSum(double val)
  {

    // could there be an issue if block has fewer than 32 threads running?
    // at least with 1 thread and warpReduceSm commneted out, we still ahve
    // chi-sq issues and worse absolute error

    val += __shfl_down_sync(0xffffffff, val, 16, WARP_SIZE);
    val += __shfl_down_sync(0xffffffff, val, 8, WARP_SIZE);
    val += __shfl_down_sync(0xffffffff, val, 4, WARP_SIZE);
    val += __shfl_down_sync(0xffffffff, val, 2, WARP_SIZE);
    val += __shfl_down_sync(0xffffffff, val, 1, WARP_SIZE);
    return val;
  }

  __inline__ __device__ double
  blockReduceSum(double val)
  {

    /*static*/ __shared__ double shared[32]; // Shared mem for 32 partial sums
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;

    val = warpReduceSum(val); // Each warp performs partial reduction

    if (lane == 0)
      shared[wid] = val; // Write reduced value to shared memory

    __syncthreads(); // Wait for all partial reductions

    // read from shared memory only if that warp existed
    val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

    if (wid == 0) {
      val = warpReduceSum(val); // Final reduce within first warp
    }
    __syncthreads(); // added by Ioannis due to cuda-memcheck racecheck
                     // reporting race between read/write
    return val;
  }

  __inline__ __device__ __host__ void
  get_indx(uint32_t ms, uint32_t* da, int ND, int NINTV)
  {
    // called like :    get_indx(m * chunkSize, &kg[1], ndim, ng);
    uint32_t dp[Internal_Vegas_Params::get_MXDIM()];
    uint32_t j, t0, t1;
    uint32_t m = ms;
    dp[0] = 1;
    dp[1] = NINTV;

    for (j = 0; j < ND - 2; j++) {
      dp[j + 2] = dp[j + 1] * NINTV;
    }

    for (j = 0; j < ND; j++) {
      t0 = dp[ND - j - 1];
      t1 = m / t0;
      da[j] = 1 + t1;
      m = m - t1 * t0;
    }
  }

  void
  Test_get_indx(int ndim,
                int ng,
                uint32_t totalNumThreads,
                int chunkSize,
                int it,
                std::ofstream& interval_myfile)
  {
    constexpr int mxdim_p1 = 21; // Internal_Vegas_Params::get_MXDIM_p1();

    if (it == 1)
      interval_myfile << "m, kg[1], kg[2], kg[3], it\n";

    for (uint32_t m = 0; m < totalNumThreads; m++) {
      uint32_t kg[mxdim_p1];
      get_indx(m, &kg[1], ndim, ng);

      interval_myfile << m << ",";
      for (int ii = 1; ii <= ndim; ii++)
        interval_myfile << kg[ii] << ",";
      interval_myfile << it << "\n";
    }
  }

  template <int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = Curand_generator>
  __inline__ __device__ void
  Setup_Integrand_Eval(Random_num_generator<GeneratorType>* rand_num_generator,
                       double xnd,
                       double dxg,
                       const double* const xi,
                       const double* const regn,
                       const double* const dx,
                       const uint32_t* const kg,
                       int* const ia,
                       double* const x,
                       double& wgt,
                       int npg,
                       int sampleID,
                       uint32_t cube_id,
                       double* randoms = nullptr)
  {
    constexpr int ndmx = Internal_Vegas_Params::get_NDMX();
    constexpr int ndmx1 = Internal_Vegas_Params::get_NDMX_p1();

#pragma unroll ndim // unroll needed
    for (int j = 1; j <= ndim; j++) {

      const double ran00 = (*rand_num_generator)();

      // if constexpr(DEBUG_MCUBES) {
      //  if(randoms != nullptr){
      //      size_t nums_per_cube = npg*ndim;
      //      size_t nums_per_sample = ndim;
      //      size_t index = cube_id*nums_per_cube +
      //      nums_per_sample*(sampleID-1) + j-1;
      //      //randoms[index] = ran00;
      //  }
      // }

      const double xn = (kg[j] - ran00) * dxg + 1.0;
      double rc = 0., xo = 0.;
      ia[j] = IMAX(IMIN((int)(xn), ndmx), 1);

      // TODO redundant memory access below xi[j * ndmx1 + ia[j] - 1]

      if (ia[j] > 1) {
        const double binA = (xi[j * ndmx1 + ia[j]]);
        const double binB = (xi[j * ndmx1 + ia[j] - 1]);
        xo = binA - binB;              // bin
                                       // length
        rc = binB + (xn - ia[j]) * xo; // scaling ran00 to bin bounds
      } else {
        xo = (xi[j * ndmx1 + ia[j]]);
        rc = (xn - ia[j]) * xo;
      }

      x[j] = regn[j] + rc * (dx[j]);
      wgt *= xo * xnd; // xnd is number of bins, xo is the length of the bin,
                       // xjac is 1/num_calls
    }
  }

  template <typename IntegT,
            int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = Curand_generator>
  __device__ void
  Process_npg_samples(IntegT* d_integrand,
                      int npg,
                      double xnd,
                      double xjac,
                      Random_num_generator<GeneratorType>* rand_num_generator,
                      double dxg,
                      const double* const regn,
                      const double* const dx,
                      const double* const xi,
                      const uint32_t* const kg,
                      int* const ia,
                      double* const x,
                      double& wgt,
                      double* d,
                      double& fb,
                      double& f2b,
                      uint32_t cube_id,
                      double* randoms = nullptr,
                      double* funcevals = nullptr)
  {
    constexpr int mxdim_p1 = Internal_Vegas_Params::get_MXDIM_p1();

    for (int k = 1; k <= npg; k++) {

      double wgt = xjac;
      Setup_Integrand_Eval<ndim, DEBUG_MCUBES, GeneratorType>(
        rand_num_generator,
        xnd,
        dxg,
        xi,
        regn,
        dx,
        kg,
        ia,
        x,
        wgt,
        npg,
        k,
        cube_id,
        randoms);

      gpu::cudaArray<double, ndim> xx;
#pragma unroll ndim // unroll needed
      for (int i = 0; i < ndim; i++) {
        xx[i] = x[i + 1];
      }

      double tmp = 0.;
      /*const double*/ tmp += gpu::apply(*d_integrand, xx);

      const double f = wgt * tmp;

      // if constexpr(DEBUG_MCUBES){
      // if constexpr(DEBUG_MCUBES){
      //     if(funcevals != nullptr){
      //         size_t nums_evals_per_cube = npg;
      //         size_t index = cube_id*nums_evals_per_cube + (k-1);
      //         //funcevals[index] = f;
      //     }
      // }

      double f2 = f * f;
      fb += f;
      f2b += f2;

      for (int j = 1; j <= ndim; j++) {
        const int index = ia[j] * mxdim_p1 + j;
        // d[ia[j] * mxdim_p1 + j] += f2; //costed less to do_eval
        atomicAdd(&d[index], f2);
      }
    }
  }

  template <typename IntegT,
            int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = Curand_generator>
  __inline__ __device__ void
  Process_chunks(IntegT* d_integrand,
                 int chunkSize,
                 int lastChunk,
                 int ng,
                 int npg,
                 Random_num_generator<GeneratorType>* rand_num_generator,
                 double dxg,
                 double xnd,
                 double xjac,
                 const double* const regn,
                 const double* const dx,
                 const double* const xi,
                 uint32_t* const kg,
                 int* const ia,
                 double* const x,
                 double& wgt,
                 double* d,
                 double& fbg,
                 double& f2bg,
                 size_t cube_id_offset,
                 double* randoms = nullptr,
                 double* funcevals = nullptr)
  {

    for (int t = 0; t < chunkSize; t++) {
      double fb = 0.,
             f2b = 0.0; // init to zero for each interval processed by thread
      uint32_t cube_id = cube_id_offset + t;

      // can't use if(is_same<GeneratorType, Custom_generator>) in device code
      // can't do if statement checking whether typename GeneratorType ==
      // Curand_generator if constexpr (mcubes::is_same<GeneratorType,
      // Custom_generator>())
      if constexpr (mcubes::TypeChecker<GeneratorType, Custom_generator>::
                      is_custom_generator()) {
        rand_num_generator->SetSeed(cube_id);
      }

      Process_npg_samples<IntegT, ndim, DEBUG_MCUBES, GeneratorType>(
        d_integrand,
        npg,
        xnd,
        xjac,
        rand_num_generator,
        dxg,
        regn,
        dx,
        xi,
        kg,
        ia,
        x,
        wgt,
        d,
        fb,
        f2b,
        cube_id,
        randoms,
        funcevals);

      f2b = sqrt(f2b * npg);
      f2b = (f2b - fb) * (f2b + fb);

      if (f2b <= 0.0) {
        f2b = TINY;
      }

      fbg += fb;
      f2bg += f2b;

      for (int k = ndim; k >= 1; k--) {
        kg[k] %= ng;

        if (++kg[k] != 1)
          break;
      }
    }
  }

  template <typename IntegT,
            int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = Curand_generator>
  __global__ void
  vegas_kernel(IntegT* d_integrand,
               int ng,
               int npg,
               double xjac,
               double dxg,
               double* result_dev,
               double xnd,
               double* xi,
               double* d,
               double* dx,
               double* regn,
               int ncubes,
               int iter,
               double sc,
               double sci,
               double ing,
               int chunkSize,
               uint32_t totalNumThreads,
               int LastChunk,
               unsigned int seed_init,
               double* randoms = nullptr,
               double* funcevals = nullptr)
  {
    constexpr int mxdim_p1 = Internal_Vegas_Params::get_MXDIM_p1();
    uint32_t m = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t tx = threadIdx.x;
    double wgt;
    uint32_t kg[mxdim_p1];
    int ia[mxdim_p1];
    double x[mxdim_p1];
    double fbg = 0., f2bg = 0.;

    if (m < totalNumThreads) {

      size_t cube_id_offset =
        (blockIdx.x * blockDim.x + threadIdx.x) * chunkSize;

      if (m == totalNumThreads - 1)
        chunkSize = LastChunk;

      Random_num_generator<GeneratorType> rand_num_generator(seed_init);
      get_indx(cube_id_offset, &kg[1], ndim, ng);

      Process_chunks<IntegT, ndim, DEBUG_MCUBES, GeneratorType>(
        d_integrand,
        chunkSize,
        LastChunk,
        ng,
        npg,
        &rand_num_generator,
        dxg,
        xnd,
        xjac,
        regn,
        dx,
        xi,
        kg,
        ia,
        x,
        wgt,
        d,
        fbg,
        f2bg,
        cube_id_offset,
        randoms,
        funcevals);
    }

    // testing if synch is needed
    __syncthreads();
    fbg = blockReduceSum(fbg);
    f2bg = blockReduceSum(f2bg);

    if (tx == 0) {
      // result_dev[0] += fbg;
      // result_dev[1] += f2bg;

      atomicAdd(&result_dev[0], fbg);
      atomicAdd(&result_dev[1], f2bg);
    }
    // end of subcube if
  }

  template <typename IntegT,
            int ndim,
            typename GeneratorType = Curand_generator>
  __global__ void
  vegas_kernelF(IntegT* d_integrand,
                int ng,
                int npg,
                double xjac,
                double dxg,
                double* result_dev,
                double xnd,
                double* xi,
                double* d,
                double* dx,
                double* regn,
                int ncubes,
                int iter,
                double sc,
                double sci,
                double ing, // not needed?
                int chunkSize,
                uint32_t totalNumThreads,
                int LastChunk,
                unsigned int seed_init)
  {

    constexpr int ndmx = 500;    // Internal_Vegas_Params::get_NDMX();
    constexpr int ndmx_p1 = 501; // Internal_Vegas_Params::get_NDMX_p1();
    constexpr int mxdim_p1 = 21; // Internal_Vegas_Params::get_MXDIM_p1();

    uint32_t m = blockIdx.x * blockDim.x + threadIdx.x;
    int tx = threadIdx.x;
    size_t cube_id_offset = (blockIdx.x * blockDim.x + threadIdx.x) * chunkSize;

    double fb, f2b, wgt, xn, xo, rc, f, f2, ran00;
    uint32_t kg[mxdim_p1];
    int iaj;
    double x[mxdim_p1];
    int k, j;
    double fbg = 0., f2bg = 0.;

    if (m < totalNumThreads) {

      if (m == totalNumThreads - 1)
        chunkSize = LastChunk;

      Random_num_generator<GeneratorType> rand_num_generator(seed_init);

      fbg = f2bg = 0.0;
      get_indx(cube_id_offset, &kg[1], ndim, ng);

      for (int t = 0; t < chunkSize; t++) {
        fb = f2b = 0.0;

        if constexpr (mcubes::TypeChecker<GeneratorType, Custom_generator>::
                        is_custom_generator()) {
          rand_num_generator.SetSeed(cube_id_offset + t);
        }

        for (k = 1; k <= npg; k++) {
          wgt = xjac;
          for (j = 1; j <= ndim; j++) {

            ran00 = rand_num_generator();
            xn = (kg[j] - ran00) * dxg + 1.0;
            iaj = IMAX(IMIN((int)(xn), ndmx), 1);

            if (iaj > 1) {
              xo = xi[j * ndmx_p1 + iaj] - xi[j * ndmx_p1 + iaj - 1];
              rc = xi[j * ndmx_p1 + iaj - 1] + (xn - iaj) * xo;
            } else {
              xo = xi[j * ndmx_p1 + iaj];
              rc = (xn - iaj) * xo;
            }

            x[j] = regn[j] + rc * dx[j];
            wgt *= xo * xnd;
          }

          double tmp;
          gpu::cudaArray<double, ndim> xx;
          for (int i = 0; i < ndim; i++) {
            xx[i] = x[i + 1];
          }

          tmp = gpu::apply(*d_integrand, xx);

          f = wgt * tmp;
          f2 = f * f;

          fb += f;
          f2b += f2;

        } // end of npg loop

        f2b = sqrt(f2b * npg);
        // double example = f2b;
        f2b = (f2b - fb) * (f2b + fb);

        if (f2b <= 0.0)
          f2b = TINY;

        fbg += fb;
        f2bg += f2b;

        for (int k = ndim; k >= 1; k--) {
          kg[k] %= ng;
          if (++kg[k] != 1)
            break;
        }

      } // end of chunk for loop
    }

    fbg = blockReduceSum(fbg);
    f2bg = blockReduceSum(f2bg);

    if (tx == 0) {
      // printf("Block %i done\n", blockIdx.x);
      // result_dev[0] += fbg;
      // result_dev[1] += f2bg;
      atomicAdd(&result_dev[0], fbg);
      atomicAdd(&result_dev[1], f2bg);
    }

    // end of subcube if
  }

  __inline__ void
  rebin(double rc, int nd, double r[], double xin[], double xi[])
  {
    int i, k = 0;
    double dr = 0.0, xn = 0.0, xo = 0.0;

    for (i = 1; i < nd; i++) {
      while (rc > dr) {
        dr += r[++k];
      }
      if (k > 1)
        xo = xi[k - 1];
      xn = xi[k];
      dr -= rc;

      xin[i] = xn - (xn - xo) * dr / r[k];
    }

    for (i = 1; i < nd; i++)
      xi[i] = xin[i];
    xi[nd] = 1.0;
  }

  template <typename IntegT,
            int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = typename ::Curand_generator>
  void
  vegas(IntegT const& integrand,
        double epsrel,
        double epsabs,
        double ncall,
        double* tgral,
        double* sd,
        double* chi2a,
        int* status,
        size_t* iters,
        int titer,
        int itmax,
        int skip,
        quad::Volume<double, ndim> const* vol)
  {
    // Mcubes_state mcubes_state(ncall, ndim);
    // all of the ofstreams below will be removed, replaced by DataLogger
    auto t0 = std::chrono::high_resolution_clock::now();

    constexpr int ndmx = 500;    // Internal_Vegas_Params::get_NDMX();
    constexpr int ndmx_p1 = 501; // Internal_Vegas_Params::get_NDMX_p1();
    constexpr int mxdim_p1 = 21; // Internal_Vegas_Params::get_MXDIM_p1();

    // IntegT* d_integrand = quad::cuda_copy_to_device(integrand);
    IntegT* d_integrand = quad::cuda_copy_to_managed(integrand);
    double regn[2 * mxdim_p1];

    for (int j = 1; j <= ndim; j++) {
      regn[j] = vol->lows[j - 1];
      regn[j + ndim] = vol->highs[j - 1];
    }

    int i, it, j, nd, ndo, ng, npg;
    double calls, dv2g, dxg, rc, ti, tsi, wgt, xjac, xn, xnd, xo;
    double k, ncubes;
    double schi, si, swgt;
    double result[2];
    double *d, *dt, *dx, *r, *x, *xi, *xin;
    int* ia;

    d =
      (double*)malloc(sizeof(double) * (ndmx_p1) * (mxdim_p1)); // contributions
    dt = (double*)malloc(sizeof(double) * (mxdim_p1));          // for cpu-only
    dx = (double*)malloc(sizeof(double) *
                         (mxdim_p1)); // length of integ-space at each dim
    r = (double*)malloc(sizeof(double) * (ndmx_p1));
    x = (double*)malloc(sizeof(double) * (mxdim_p1));
    xi = (double*)malloc(sizeof(double) * (mxdim_p1) *
                         (ndmx_p1)); // right bin coord
    xin = (double*)malloc(sizeof(double) * (ndmx_p1));
    ia = (int*)malloc(sizeof(int) * (mxdim_p1));

    // code works only  for (2 * ng - NDMX) >= 0)

    ndo = 1;
    for (j = 1; j <= ndim; j++) {
      xi[j * ndmx_p1 + 1] =
        1.0; // this index is the first for each bin for each dimension
    }

    si = swgt = schi = 0.0;
    nd = ndmx;
    ng = 1;
    ng = (int)pow(ncall / 2.0 + 0.25, 1.0 / ndim); // why do we add .25?
    for (k = 1, i = 1; i < ndim; i++) {
      k *= ng;
    }

    double sci = 1.0 / k; // I dont' think that's used anywhere
    double sc = k;        // I dont' think that's used either
    k *= ng;
    ncubes = k;

    npg = IMAX(ncall / k, 2);
    // assert(npg == Compute_samples_per_cube(ncall, ncubes)); //to replace line
    // directly above assert(ncubes == ComputeNcubes(ncall, ndim)); //to replace
    // line directly above

    calls = (double)npg * (double)k;
    dxg = 1.0 / ng;

    double ing = dxg;
    for (dv2g = 1, i = 1; i <= ndim; i++)
      dv2g *= dxg;
    dv2g = (calls * dv2g * calls * dv2g) / npg / npg / (npg - 1.0);

    xnd = nd;
    dxg *= xnd;
    xjac = 1.0 / calls;
    for (j = 1; j <= ndim; j++) {
      dx[j] = regn[j + ndim] - regn[j];
      xjac *= dx[j];
    }

    for (i = 1; i <= IMAX(nd, ndo); i++)
      r[i] = 1.0;
    for (j = 1; j <= ndim; j++) {
      rebin(ndo / xnd, nd, r, xin, &xi[j * ndmx_p1]);
    }

    ndo = nd;

    double *d_dev, *dx_dev, *x_dev, *xi_dev, *regn_dev, *result_dev;
    int* ia_dev;

    cudaMalloc((void**)&result_dev, sizeof(double) * 2);
    cudaCheckError();
    cudaMalloc((void**)&d_dev, sizeof(double) * (ndmx_p1) * (mxdim_p1));
    cudaCheckError();
    cudaMalloc((void**)&dx_dev, sizeof(double) * (mxdim_p1));
    cudaCheckError();
    cudaMalloc((void**)&x_dev, sizeof(double) * (mxdim_p1));
    cudaCheckError();
    cudaMalloc((void**)&xi_dev, sizeof(double) * (mxdim_p1) * (ndmx_p1));
    cudaCheckError();
    cudaMalloc((void**)&regn_dev, sizeof(double) * ((ndim * 2) + 1));
    cudaCheckError();
    cudaMalloc((void**)&ia_dev, sizeof(int) * (mxdim_p1));
    cudaCheckError();

    cudaMemcpy(dx_dev, dx, sizeof(double) * (mxdim_p1), cudaMemcpyHostToDevice);
    cudaCheckError();
    cudaMemcpy(x_dev, x, sizeof(double) * (mxdim_p1), cudaMemcpyHostToDevice);
    cudaCheckError();
    cudaMemcpy(regn_dev,
               regn,
               sizeof(double) * ((ndim * 2) + 1),
               cudaMemcpyHostToDevice);
    cudaCheckError();

    cudaMemset(ia_dev, 0, sizeof(int) * (mxdim_p1));

    int chunkSize = GetChunkSize(ncall);
    // chunkSize = 4;
    uint32_t totalNumThreads =
      (uint32_t)((ncubes /*+ chunkSize - 1*/) / chunkSize);

    uint32_t totalCubes = totalNumThreads * chunkSize; // even-split cubes
    int extra = ncubes - totalCubes;                   // left-over cubes
    int LastChunk = extra + chunkSize; // last chunk of last thread

    Kernel_Params params(ncall, chunkSize, ndim);
    // uint32_t nBlocks =
    //   ((uint32_t)(((ncubes + BLOCK_DIM_X - 1) / BLOCK_DIM_X)) / chunkSize) +
    //   1; // compute blocks based on chunk_size, ncubes, and block_dim_x

    // uint32_t nBlocks = ncubes % BLOCK_DIM_X == 0 ? (ncubes/BLOCK_DIM_X) :
    // (ncubes/BLOCK_DIM_X) + 1; uint32_t nThreads = BLOCK_DIM_X;

    IterDataLogger<DEBUG_MCUBES> data_collector(
      totalNumThreads, chunkSize, extra, npg, ndim);

    for (it = 1; it <= itmax && (*status) == 1; (*iters)++, it++) {
      ti = tsi = 0.0;
      for (j = 1; j <= ndim; j++) {
        for (i = 1; i <= nd; i++)
          d[i * mxdim_p1 + j] = 0.0;
      }

      cudaMemcpy(xi_dev,
                 xi,
                 sizeof(double) * (mxdim_p1) * (ndmx_p1),
                 cudaMemcpyHostToDevice);
      cudaCheckError(); // bin bounds
      cudaMemset(
        d_dev, 0, sizeof(double) * (ndmx_p1) * (mxdim_p1)); // bin contributions
      cudaMemset(result_dev, 0, 2 * sizeof(double));

      using MilliSeconds =
        std::chrono::duration<double, std::chrono::milliseconds::period>;

      MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
      unsigned int seed = /*static_cast<unsigned int>(time_diff.count()) +*/
        static_cast<unsigned int>(it);

      cudaEvent_t start, stop;
      cudaEventCreate(&start);
      cudaEventCreate(&stop);

      cudaProfilerStart();
      cudaEventRecord(start);
      vegas_kernel<IntegT, ndim, DEBUG_MCUBES, GeneratorType>
        <<<params.nBlocks, params.nThreads>>>(d_integrand,
                                              ng,
                                              npg,
                                              xjac,
                                              dxg,
                                              result_dev,
                                              xnd,
                                              xi_dev,
                                              d_dev,
                                              dx_dev,
                                              regn_dev,
                                              ncubes,
                                              it,
                                              sc,
                                              sci,
                                              ing,
                                              chunkSize,
                                              totalNumThreads,
                                              LastChunk,
                                              seed + it,
                                              data_collector.randoms,
                                              data_collector.funcevals);
      cudaDeviceSynchronize();
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      float kernel_time = 0;
      cudaEventElapsedTime(&kernel_time, start, stop);
      // std::cout<< "vegas_kernel:" << params.nBlocks << "," << kernel_time <<
      // std::endl; std::cout<<"\textra cubes:"<< extra << "," << "time:" <<
      // kernel_time << std::endl;
      params.set_alternate(chunkSize);
      cudaProfilerStop();
      cudaMemcpy(xi,
                 xi_dev,
                 sizeof(double) * (mxdim_p1) * (ndmx_p1),
                 cudaMemcpyDeviceToHost);
      cudaCheckError();

      cudaMemcpy(d,
                 d_dev,
                 sizeof(double) * (ndmx_p1) * (mxdim_p1),
                 cudaMemcpyDeviceToHost);

      cudaCheckError(); // we do need to the contributions for the rebinning
      cudaMemcpy(
        result, result_dev, sizeof(double) * 2, cudaMemcpyDeviceToHost);

      ti = result[0];
      tsi = result[1];
      tsi *= dv2g;

      if (it > skip) {
        wgt = 1.0 / tsi;
        si += wgt * ti;
        schi += wgt * ti * ti;
        swgt += wgt;
        *tgral = si / swgt;
        *chi2a = (schi - si * (*tgral)) / (static_cast<double>(it) - 0.9999);
        if (*chi2a < 0.0)
          *chi2a = 0.0;
        *sd = sqrt(1.0 / swgt);
        tsi = sqrt(tsi);
        *status = GetStatus(*tgral, *sd, it, epsrel, epsabs);
        // printf("%i %.15f +- %.15f iteration: %.15f +- %.15f chi:%.15f\n", it,
        // *tgral, *sd, ti, sqrt(tsi), *chi2a);
      }

      if constexpr (DEBUG_MCUBES == true) {
        data_collector.PrintBins(it, xi, d, ndim);
        // data_collector.PrintRandomNums(it, ncubes, npg, ndim);
        // data_collector.PrintFuncEvals(it, ncubes, npg, ndim);
        data_collector.PrintIterResults(it, *tgral, *sd, *chi2a, ti, tsi);
      }

      // replace above with datalogger.print();
      for (j = 1; j <= ndim; j++) {
        xo = d[1 * mxdim_p1 + j]; // bin 1 of dim j
        xn = d[2 * mxdim_p1 + j]; // bin 2 of dim j
        d[1 * mxdim_p1 + j] = (xo + xn) / 2.0;
        dt[j] = d[1 * mxdim_p1 + j]; // set dt sum to contribution of bin 1

        for (i = 2; i < nd; i++) {
          rc = xo + xn;
          xo = xn;
          xn = d[(i + 1) * mxdim_p1 + j];
          d[i * mxdim_p1 + j] = (rc + xn) / 3.0;
          dt[j] += d[i * mxdim_p1 + j];
        }

        d[nd * mxdim_p1 + j] = (xo + xn) / 2.0; // do bin nd last
        dt[j] += d[nd * mxdim_p1 + j];
      }

      for (j = 1; j <= ndim; j++) {
        if (dt[j] > 0.0) { // enter if there is any contribution only
          rc = 0.0;
          for (i = 1; i <= nd; i++) {
            // if(d[i*mxdim_p1+j]<TINY) d[i*mxdim_p1+j]=TINY;
            // if(d[i*mxdim_p1+j]<TINY) printf("d[%i]:%.15e\n", i*mxdim_p1+j,
            // d[i*mxdim_p1+j]); printf("d[%i]:%.15e\n", i*mxdim_p1+j,
            // d[i*mxdim_p1+j]);
            r[i] = pow((1.0 - d[i * mxdim_p1 + j] / dt[j]) /
                         (log(dt[j]) - log(d[i * mxdim_p1 + j])),
                       1.5 /*Internal_Vegas_Params::get_ALPH()*/);
            rc += r[i]; // rc is it the total number of sub-increments
          }
          rebin(
            rc / xnd,
            nd,
            r,
            xin,
            &xi[j * ndmx_p1]); // first bin of each dimension is at a diff index
        }
      }

    } // end of iterations

    //  Start of iterations without adjustment

    cudaMemcpy(xi_dev,
               xi,
               sizeof(double) * (mxdim_p1) * (ndmx_p1),
               cudaMemcpyHostToDevice);
    cudaCheckError();

    for (it = itmax + 1; it <= titer && (*status); (*iters)++, it++) {
      ti = tsi = 0.0;
      cudaMemset(result_dev, 0, 2 * sizeof(double));

      using MilliSeconds =
        std::chrono::duration<double, std::chrono::milliseconds::period>;
      MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
      unsigned int seed = static_cast<unsigned int>(time_diff.count()) +
                          static_cast<unsigned int>(it);
      vegas_kernelF<IntegT, ndim, GeneratorType>
        <<<params.nBlocks, params.nThreads>>>(d_integrand,
                                              ng,
                                              npg,
                                              xjac,
                                              dxg,
                                              result_dev,
                                              xnd,
                                              xi_dev,
                                              d_dev,
                                              dx_dev,
                                              regn_dev,
                                              ncubes,
                                              it,
                                              sc,
                                              sci,
                                              ing,
                                              chunkSize,
                                              totalNumThreads,
                                              LastChunk,
                                              seed + it);
      cudaDeviceSynchronize();
      cudaMemcpy(
        result, result_dev, sizeof(double) * 2, cudaMemcpyDeviceToHost);

      ti = result[0];
      tsi = result[1];
      tsi *= dv2g;
      // printf("iter %d  integ = %.15e   std = %.15e var:%.15e dv2g:%f\n", it,
      // ti, sqrt(tsi), tsi, dv2g);

      wgt = 1.0 / tsi;
      si += wgt * ti;
      schi += wgt * ti * ti;
      swgt += wgt;
      *tgral = si / swgt;

      *chi2a = (schi - si * (*tgral)) / (it - 0.9999);
      if (*chi2a < 0.0)
        *chi2a = 0.0;
      *sd = sqrt(1.0 / swgt);
      tsi = sqrt(tsi);
      *status = GetStatus(*tgral, *sd, it, epsrel, epsabs);
      // printf("%i, %.15f,  %.15f, %.15f, %.15f, %.15f\n", it, *tgral, *sd, ti,
      // sqrt(tsi), *chi2a);
      if constexpr (DEBUG_MCUBES) {
        data_collector.PrintBins(it, xi, d, ndim);
        // data_collector.PrintRandomNums(it, ncubes, npg, ndim);
        // data_collector.PrintFuncEvals(it, ncubes, npg, ndim);
        data_collector.PrintIterResults(it, *tgral, *sd, *chi2a, ti, tsi);
      }
      // printf("cummulative ti:%5d, integral: %.15e, sd:%.4e,chi_sq:%9.2g\n",
      // it, *tgral, *sd, *chi2a);
    } // end of iterations

    free(d);
    free(dt);
    free(dx);
    free(ia);
    free(x);
    free(xi);
    free(xin);
    free(r);

    d_integrand->~IntegT();
    cudaFree(d_dev);
    cudaFree(dx_dev);
    cudaFree(ia_dev);
    cudaFree(x_dev);
    cudaFree(xi_dev);
    cudaFree(regn_dev);
    cudaFree(result_dev);
    cudaFree(d_integrand);
  }

  template <typename IntegT,
            int NDIM,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = typename ::Curand_generator>
  numint::integration_result
  integrate(IntegT& ig,
            double epsrel,
            double epsabs,
            double ncall,
            quad::Volume<double, NDIM> const* volume,
            int totalIters = 15,
            int adjustIters = 15,
            int skipIters = 5)
  {

    numint::integration_result result;
    result.status = 1;
    vegas<IntegT, NDIM, DEBUG_MCUBES, GeneratorType>(ig,
                                                     epsrel,
                                                     epsabs,
                                                     ncall,
                                                     &result.estimate,
                                                     &result.errorest,
                                                     &result.chi_sq,
                                                     &result.status,
                                                     &result.iters,
                                                     totalIters,
                                                     adjustIters,
                                                     skipIters,
                                                     volume);
    return result;
  }

  class mcubes {
  public:
    mcubes(double epsrel,
           double epsabs,
           double ncall,
           int totalIters = 15,
           int adjustIters = 15,
           int skipIters = 5)
      : epsrel_(epsrel)
      , epsabs_(epsabs)
      , ncall_(ncall)
      , totalIters_(totalIters)
      , adjustIters_(adjustIters)
      , skipIters_(skipIters)
    {}

    template <typename F, int NDIM>
    numint::integration_result integrate(
      F& integrand,
      quad::Volume<double, NDIM> const* volume);

  private:
    double epsrel_;
    double epsabs_;
    double ncall_;
    int totalIters_;
    int adjustIters_;
    int skipIters_;
  };

  template <typename F, int NDIM>
  numint::integration_result
  mcubes::integrate(F& integrand, quad::Volume<double, NDIM> const* volume)
  {
    return mcubes::integrate<F, NDIM>(integrand,
                                      epsrel_,
                                      epsabs_,
                                      ncall_,
                                      volume,
                                      totalIters_,
                                      adjustIters_,
                                      skipIters_);
  }

  template <typename IntegT,
            int NDIM,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = typename ::Curand_generator>
  numint::integration_result
  simple_integrate(IntegT const& integrand,
                   double epsrel,
                   double epsabs,
                   double ncall,
                   quad::Volume<double, NDIM> const* volume,
                   int totalIters = 15,
                   int adjustIters = 15,
                   int skipIters = 5)
  {

    numint::integration_result result;
    result.status = 1;

    do {
      vegas<IntegT, NDIM, DEBUG_MCUBES, GeneratorType>(integrand,
                                                       epsrel,
                                                       epsabs,
                                                       ncall,
                                                       &result.estimate,
                                                       &result.errorest,
                                                       &result.chi_sq,
                                                       &result.status,
                                                       &result.iters,
                                                       totalIters,
                                                       adjustIters,
                                                       skipIters,
                                                       volume);
    } while (result.status == 1 && AdjustParams(ncall, totalIters) == true);

    return result;
  }

}
#endif
