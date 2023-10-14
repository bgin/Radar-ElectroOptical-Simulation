#ifndef VEGAS_T_1D_CUH
#define VEGAS_T_1D_CUH

/*

nvcc -O2 -DCUSTOM -o vegas1D vegas_mcubes_1D.cu -arch=sm_70
OR
nvcc -O2 -DCURAND -o vegas1D vegas_mcubes_1D.cu -arch=sm_70

example run command

nvprof ./vegas1D 0 6 0.0  10.0  2.0E+09  58, 0, 0

nvprof  ./vegas1D 1 9 -1.0  1.0  1.0E+07 15 10 10

nvprof ./vegas1D 2 2 -1.0 1.0  1.0E+09 1 0 0

Last three arguments are: total iterations, iteration

 */
#include <chrono>
#include <stdio.h>
// #include <malloc.h>
#include "Volume.cuh"
#include "cudaApply.cuh"
#include "cudaArray.cuh"
#include "vegas_utils.cuh"
#include "verbose_utils.cuh"
#include <ctime>
#include <curand_kernel.h>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "integration_result.hh"

#define WARP_SIZE 32
#define BLOCK_DIM_X 128

namespace mcubes1D {

  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;

  __inline__ __device__ double
  warpReduceSum(double val)
  {
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

    static __shared__ double shared[32]; // Shared mem for 32 partial sums
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;

    val = warpReduceSum(val); // Each warp performs partial reduction

    if (lane == 0)
      shared[wid] = val; // Write reduced value to shared memory

    __syncthreads(); // Wait for all partial reductions

    // read from shared memory only if that warp existed
    val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

    if (wid == 0)
      val = warpReduceSum(val); // Final reduce within first warp

    return val;
  }

  __inline__ __device__ void
  get_indx(int ms, uint32_t* da, int ND, int NINTV)
  {
    int dp[/*Internal_Vegas_Params::get_MXDIM()*/ 21];
    int j, t0, t1;
    int m = ms;
    dp[0] = 1;
    dp[1] = NINTV;

    for (j = 0; j < ND - 2; j++) {
      dp[j + 2] = dp[j + 1] * NINTV;
    }
    //
    for (j = 0; j < ND; j++) {
      t0 = dp[ND - j - 1];
      t1 = m / t0;
      da[j] = 1 + t1;
      m = m - t1 * t0;
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
               double ncubes,
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

    constexpr int ndmx_p1 = Internal_Vegas_Params::get_NDMX_p1();
    constexpr int mxdim_p1 = Internal_Vegas_Params::get_MXDIM_p1();
    constexpr int ndmx = Internal_Vegas_Params::get_NDMX();

    uint32_t m = blockIdx.x * blockDim.x + threadIdx.x;
    int tx = threadIdx.x;

    double fb, f2b /*, wgt, xn, xo, rc, f, f2*/;
    uint32_t kg[mxdim_p1];
    int ia[mxdim_p1];
    double x[mxdim_p1];
    int k, j;
    double fbg = 0., f2bg = 0.;
    gpu::cudaArray<double, ndim> xx;

    if (m < totalNumThreads) {
      Random_num_generator<GeneratorType> rand_num_generator(seed_init);

      get_indx(m * chunkSize, &kg[1], ndim, ng);
      size_t cube_id_offset = m * chunkSize;
      if (m == totalNumThreads - 1)
        chunkSize = LastChunk;

      fbg = f2bg = 0.0;
      for (int t = 0; t < chunkSize; t++) {
        fb = f2b = 0.0;

        if constexpr (mcubes::TypeChecker<GeneratorType, Custom_generator>::
                        is_custom_generator()) {
          // size_t cube_id_offset = m*chunkSize;
          uint32_t cube_id = cube_id_offset + t;
          rand_num_generator.SetSeed(cube_id);
        }

        for (k = 1; k <= npg; k++) {
          double wgt = xjac;
          for (j = 1; j <= ndim; j++) {

            const double ran00 = rand_num_generator();

            if constexpr (DEBUG_MCUBES) {
              if (randoms != nullptr) {
                // size_t cube_id_offset = m*chunkSize;
                uint32_t cube_id = cube_id_offset + t;
                size_t nums_per_cube = npg * ndim;
                size_t nums_per_sample = ndim;
                size_t index =
                  cube_id * nums_per_cube + nums_per_sample * (k - 1) + j - 1;
                randoms[index] = ran00;
              }
            }

            const double xn = (kg[j] - ran00) * dxg + 1.0;
            double rc = 0., xo = 0.;
            ia[j] = IMAX(IMIN((int)(xn), ndmx), 1);

            if (ia[j] > 1) {
              xo = xi[j * ndmx_p1 + ia[j]] - xi[j * ndmx_p1 + ia[j] - 1];
              rc = xi[j * ndmx_p1 + ia[j] - 1] + (xn - ia[j]) * xo;
            } else {
              xo = xi[j * ndmx_p1 + ia[j]];
              rc = (xn - ia[j]) * xo;
            }

            x[j] = regn[j] + rc * dx[j];
            wgt *= xo * xnd;
          }

          for (int dim = 0; dim <= ndim; dim++)
            xx[dim] = x[dim + 1];

          const double tmp = gpu::apply(*d_integrand, xx);
          const double f = wgt * tmp;

          if constexpr (DEBUG_MCUBES) {
            if (funcevals != nullptr) {
              // size_t cube_id_offset = m*chunkSize;
              uint32_t cube_id = cube_id_offset + t;
              size_t nums_evals_per_cube = npg;
              size_t index = cube_id * nums_evals_per_cube + (k - 1);
              funcevals[index] = f;
            }
          }

          double f2 = f * f;

          fb += f;
          f2b += f2;

          atomicAdd(&d[ia[1] * mxdim_p1 + 1], f2);
        } // end of npg loop

        f2b = sqrt(f2b * npg);
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
    }   // end of subcube if

    fbg = blockReduceSum(fbg);
    f2bg = blockReduceSum(f2bg);

    if (tx == 0) {
      atomicAdd(&result_dev[0], fbg);
      atomicAdd(&result_dev[1], f2bg);
    }
  }

  template <typename IntegT,
            int ndim,
            bool DEBUG_MCUBES = false,
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
                double ncubes,
                int iter,
                double sc,
                double sci,
                double ing,
                int chunkSize,
                uint32_t totalNumThreads,
                int LastChunk,
                unsigned int seed_init)
  {

    constexpr int ndmx_p1 = Internal_Vegas_Params::get_NDMX_p1();
    constexpr int ndmx = Internal_Vegas_Params::get_NDMX();
    constexpr int mxdim_p1 = Internal_Vegas_Params::get_MXDIM_p1();

    uint32_t m = blockIdx.x * blockDim.x + threadIdx.x;
    int tx = threadIdx.x;

    double fb, f2b, wgt, xn, xo, rc, f, f2, ran00;
    uint32_t kg[mxdim_p1];
    int iaj;
    double x[mxdim_p1];
    int k, j;
    double fbg = 0., f2bg = 0.;
    gpu::cudaArray<double, ndim> xx;

    if (m < totalNumThreads) {
      get_indx(m * chunkSize, &kg[1], ndim, ng);

      if (m == totalNumThreads - 1)
        chunkSize = LastChunk;
      Random_num_generator<GeneratorType> rand_num_generator(seed_init);

      fbg = f2bg = 0.0;
      for (int t = 0; t < chunkSize; t++) {
        fb = f2b = 0.0;

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

          for (int i = 0; i < ndim; i++) {
            xx[i] = x[i + 1];
          }
          double tmp = gpu::apply(*d_integrand, xx);
          f = wgt * tmp;
          f2 = f * f;

          fb += f;
          f2b += f2;

        } // end of npg loop

        f2b = sqrt(f2b * npg);
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
    }   // end of subcube if

    fbg = blockReduceSum(fbg);
    f2bg = blockReduceSum(f2bg);

    if (tx == 0) {
      atomicAdd(&result_dev[0], fbg);
      atomicAdd(&result_dev[1], f2bg);
    }
  }

  void
  rebin(double rc, int nd, double r[], double xin[], double xi[])
  {
    int i, k = 0;
    double dr = 0.0, xn = 0.0, xo = 0.0;
    for (i = 1; i < nd; i++) {
      while (rc > dr)
        dr += r[++k];
      if (k > 1)
        xo = xi[k - 1];
      xn = xi[k];
      dr -= rc;
      xin[i] = xn - (xn - xo) * dr / r[k];
    }

    for (i = 1; i < nd; i++)
      xi[i] = xin[i];
    xi[nd] = 1.0;
    // for (i=1;i<=nd;i++) printf("bins edges: %.10f\n", xi[i]);
    // printf("---------------------\n");
  }

  template <typename IntegT,
            int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = Curand_generator>
  void
  vegas1D(IntegT integrand,
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
    auto t0 = std::chrono::high_resolution_clock::now();

    constexpr int mxdim_p1 = Internal_Vegas_Params::get_MXDIM_p1();
    constexpr int ndmx = Internal_Vegas_Params::get_NDMX();
    constexpr int ndmx_p1 = Internal_Vegas_Params::get_NDMX_p1();

    double regn[2 * mxdim_p1];
    int i, it, j /*, k*/, nd, ndo, ng, npg /*, ncubes*/;
    double k, ncubes;
    double calls, dv2g, dxg, rc, ti, tsi, wgt, xjac, xn, xnd, xo;

    for (int j = 1; j <= ndim; j++) {
      regn[j] = vol->lows[j - 1];
      regn[j + ndim] = vol->highs[j - 1];
    }

    double schi, si, swgt;
    double result[2];
    double *d, *dt, *dx, *r, *x, *xi, *xin;
    int* ia;

    d = (double*)malloc(sizeof(double) * (ndmx_p1) * (mxdim_p1));
    dt = (double*)malloc(sizeof(double) * (mxdim_p1));
    dx = (double*)malloc(sizeof(double) * (mxdim_p1));
    r = (double*)malloc(sizeof(double) * (ndmx_p1));
    x = (double*)malloc(sizeof(double) * (mxdim_p1));
    xi = (double*)malloc(sizeof(double) * (mxdim_p1) * (ndmx_p1));
    xin = (double*)malloc(sizeof(double) * (ndmx_p1));
    ia = (int*)malloc(sizeof(int) * (mxdim_p1));

    // code works only  for (2 * ng - NDMX) >= 0)

    ndo = 1;
    for (j = 1; j <= ndim; j++)
      xi[j * ndmx_p1 + 1] = 1.0;

    si = swgt = schi = 0.0;
    nd = ndmx;
    ng = 1;
    ng = (int)pow(ncall / 2.0 + 0.25, 1.0 / ndim);
    for (k = 1, i = 1; i < ndim; i++)
      k *= ng;

    double sci = 1.0 / k;
    double sc = k;
    k *= ng;

    ncubes = k;

    npg = IMAX(ncall / k, 2);
    calls = (double)npg * k;
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
    for (j = 1; j <= ndim; j++)
      rebin(ndo / xnd, nd, r, xin, &xi[j * ndmx_p1]);

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

    uint32_t totalNumThreads = (uint32_t)((ncubes) / chunkSize);
    uint32_t totalCubes = totalNumThreads * chunkSize;
    int extra = ncubes - totalCubes;
    int LastChunk = extra + chunkSize;
    uint32_t nBlocks =
      ((uint32_t)(((ncubes + BLOCK_DIM_X - 1) / BLOCK_DIM_X)) / chunkSize) + 1;
    uint32_t nThreads = BLOCK_DIM_X;

    IterDataLogger<DEBUG_MCUBES> data_collector(
      totalNumThreads, chunkSize, extra, npg, ndim);
    IntegT* d_integrand = cuda_copy_to_managed(integrand);

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
      cudaCheckError();

      cudaMemset(d_dev, 0, sizeof(double) * (ndmx_p1) * (mxdim_p1));
      cudaMemset(result_dev, 0, 2 * sizeof(double));

      MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
      unsigned int seed = static_cast<unsigned int>(time_diff.count()) +
                          static_cast<unsigned int>(it);

      vegas_kernel<IntegT, ndim, DEBUG_MCUBES, GeneratorType>
        <<<nBlocks, nThreads>>>(d_integrand,
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

      cudaMemcpy(xi,
                 xi_dev,
                 sizeof(double) * (mxdim_p1) * (ndmx_p1),
                 cudaMemcpyDeviceToHost);
      cudaCheckError();

      cudaMemcpy(d,
                 d_dev,
                 sizeof(double) * (ndmx_p1) * (mxdim_p1),
                 cudaMemcpyDeviceToHost);
      cudaCheckError();
      cudaMemcpy(
        result, result_dev, sizeof(double) * 2, cudaMemcpyDeviceToHost);

      ti = result[0];
      tsi = result[1];
      tsi *= dv2g;
      // printf("iter = %d  integ = %e   std = %e\n", it, ti, sqrt(tsi));

      if (it > skip) {
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
        // printf("it %d\n", it);
        // printf("%5d,%14.7g,%9.2g,%9.2g\n", it, *tgral, *sd, *chi2a);
        *status = GetStatus(*tgral, *sd, it, epsrel, epsabs);
      }
      // printf("%3d   %e  %e\n", it, ti, tsi);

      if constexpr (DEBUG_MCUBES == true) {
        data_collector.PrintBins(it, xi, d, ndim);
        // data_collector.PrintRandomNums(it, ncubes, npg, ndim);
        // data_collector.PrintFuncEvals(it, ncubes, npg, ndim);
        data_collector.PrintIterResults(it, *tgral, *sd, *chi2a, ti, tsi);
      }

      for (j = 1; j <= 1; j++) {
        xo = d[1 * mxdim_p1 + j];
        xn = d[2 * mxdim_p1 + j];
        d[1 * mxdim_p1 + j] = (xo + xn) / 2.0;
        dt[j] = d[1 * mxdim_p1 + j];

        for (i = 2; i < nd; i++) {
          rc = xo + xn;
          xo = xn;
          xn = d[(i + 1) * mxdim_p1 + j];
          d[i * mxdim_p1 + j] = (rc + xn) / 3.0;
          dt[j] += d[i * mxdim_p1 + j];
        }

        d[nd * mxdim_p1 + j] = (xo + xn) / 2.0;
        dt[j] += d[nd * mxdim_p1 + j];
      }

      for (j = 1; j <= 1; j++) {
        if (dt[j] > 0.0) {
          rc = 0.0;
          for (i = 1; i <= nd; i++) {
            // if (d[i * mxdim_p1 + j] < TINY) d[i * mxdim_p1 + j] = TINY;
            r[i] = pow((1.0 - d[i * mxdim_p1 + j] / dt[j]) /
                         (log(dt[j]) - log(d[i * mxdim_p1 + j])),
                       Internal_Vegas_Params::get_ALPH());
            rc += r[i];
          }
          rebin(rc / xnd, nd, r, xin, &xi[j * ndmx_p1]);
        }
      }

      for (j = 2; j <= ndim; j++) {
        for (i = 1; i <= nd; i++) {
          xi[j * ndmx_p1 + i] = xi[ndmx_p1 + i];
        }
      }

    } // end of iteration loop

    cudaMemcpy(xi_dev,
               xi,
               sizeof(double) * (mxdim_p1) * (ndmx_p1),
               cudaMemcpyHostToDevice);
    cudaCheckError();

    for (it = itmax + 1; it <= titer && (*status) == 1; (*iters)++, it++) {
      ti = tsi = 0.0;
      cudaMemset(result_dev, 0, 2 * sizeof(double));

      using MilliSeconds =
        std::chrono::duration<double, std::chrono::milliseconds::period>;
      MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
      unsigned int seed = static_cast<unsigned int>(time_diff.count()) +
                          static_cast<unsigned int>(it);
      vegas_kernelF<IntegT, ndim, DEBUG_MCUBES, GeneratorType>
        <<<nBlocks, nThreads>>>(d_integrand,
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

      if constexpr (DEBUG_MCUBES == true) {
        data_collector.PrintBins(it, xi, d, ndim);
        // data_collector.PrintRandomNums(it, ncubes, npg, ndim);
        // data_collector.PrintFuncEvals(it, ncubes, npg, ndim);
        data_collector.PrintIterResults(it, *tgral, *sd, *chi2a, ti, tsi);
      }

      //  printf("%5d,%14.7g,%9.4g,%9.2g\n", it, *tgral, *sd, *chi2a);

    } // end of iteration

    free(d);
    free(dt);
    free(dx);
    free(ia);
    free(x);
    free(xi);
    free(xin);
    free(r);

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
            int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = Curand_generator>
  numint::integration_result
  integrate1D(IntegT integrand,
              double epsrel,
              double epsabs,
              double ncall,
              quad::Volume<double, ndim> const* volume,
              int titer,
              int itmax,
              int skip)
  {
    numint::integration_result res;
    res.status = 1;

    vegas1D<IntegT, ndim, DEBUG_MCUBES, GeneratorType>(integrand,
                                                       epsrel,
                                                       epsabs,
                                                       ncall,
                                                       &res.estimate,
                                                       &res.errorest,
                                                       &res.chi_sq,
                                                       &res.status,
                                                       &res.iters,
                                                       titer,
                                                       itmax,
                                                       skip,
                                                       volume);
    return res;
  }

  template <typename IntegT,
            int ndim,
            bool DEBUG_MCUBES = false,
            typename GeneratorType = Curand_generator>
  numint::integration_result
  simple_integrate1D(IntegT integrand,
                     double epsrel,
                     double epsabs,
                     double ncall,
                     quad::Volume<double, ndim> const* volume,
                     int titer,
                     int itmax,
                     int skip)
  {
    numint::integration_result res;
    res.status = 1;
    do {
      vegas1D<IntegT, ndim, DEBUG_MCUBES, GeneratorType>(integrand,
                                                         epsrel,
                                                         epsabs,
                                                         ncall,
                                                         &res.estimate,
                                                         &res.errorest,
                                                         &res.chi_sq,
                                                         &res.status,
                                                         &res.iters,
                                                         titer,
                                                         itmax,
                                                         skip,
                                                         volume);
    } while (res.status == 1 && AdjustParams(ncall, titer) == true);

    return res;
  }

}

#endif
