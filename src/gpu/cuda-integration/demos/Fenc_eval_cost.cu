#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"
#include <chrono>
#include <iostream>
#include <string>
#include "math.h"

class SinSum6D {
public:
  __device__ __host__ double
  operator()(double x, double y, double z, double k, double l, double m)
  {
    return sin(x + y + z + k + l + m);
  }
};

class Gauss9D {
public:
  __device__ __host__ double
  operator()(double x,
             double y,
             double z,
             double k,
             double l,
             double m,
             double n,
             double o,
             double p)
  {
    double sum = pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(k, 2) + pow(l, 2) +
                 pow(m, 2) + pow(n, 2) + pow(o, 2) + pow(p, 2);
    return exp(-1 * sum / (2 * pow(0.01, 2))) *
           (1 / pow(sqrt(2 * M_PI) * 0.01, 9));
  }
};

class GENZ_3_8D {
public:
  __device__ __host__ double
  operator()(double x,
             double y,
             double z,
             double w,
             double v,
             double u,
             double t,
             double s)
  {
    return pow(1 + 8 * s + 7 * t + 6 * u + 5 * v + 4 * w + 3 * x + 2 * y + z,
               -9);
  }
};

class GENZ_4_8D {
public:
  __device__ __host__ double
  operator()(double x,
             double y,
             double z,
             double w,
             double v,
             double k,
             double m,
             double n)
  {
    double beta = .5;
    return exp(-1.0 *
               (pow(25, 2) * pow(x - beta, 2) + pow(25, 2) * pow(y - beta, 2) +
                pow(25, 2) * pow(z - beta, 2) + pow(25, 2) * pow(w - beta, 2) +
                pow(25, 2) * pow(v - beta, 2) + pow(25, 2) * pow(k - beta, 2) +
                pow(25, 2) * pow(m - beta, 2) + pow(25, 2) * pow(n - beta, 2)));
  }
};

class GENZ_6_6D {
public:
  __device__ __host__ double
  operator()(double u, double v, double w, double x, double y, double z)
  {
    if (z > .9 || y > .8 || x > .7 || w > .6 || v > .5 || u > .4)
      return 0.;
    else
      return exp(10 * z + 9 * y + 8 * x + 7 * w + 6 * v +
                 5 * u) /*/1.5477367885091207413e8*/;
  }
};

double
num_sub_cubes(const double ncall, const int ndim)
{

  const int ng =
    static_cast<int>(pow(ncall / 2.0 + 0.25, 1.0 / static_cast<double>(ndim)));
  int k = 1;
  int i = 1;
  for (k = 1, i = 1; i < ndim; i++) {
    k *= ng;
  }

  k *= ng;
  const double ncubes = static_cast<double>(k);
  return ncubes;
};

int
samples_per_cube(double ncall, double ncubes)
{
  return static_cast<int>(IMAX(ncall / ncubes, 2));
}

double
num_func_evals(size_t ncall, int ndim)
{
  const double ncubes = num_sub_cubes(ncall, ndim);
  const int spc = samples_per_cube(ncall, ncubes);
  return ncubes * spc;
}

double
num_random_nums(double ncall, int ndim)
{
  const double fevals = num_func_evals(ncall, ndim);
  return fevals * static_cast<double>(ndim);
}

__global__ void
generate_random_nums(double* rands,
                     double min,
                     double max,
                     unsigned int seed_init)
{
  size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

  // printf("tid:%lu\n", tid);
  Random_num_generator<Curand_generator> rand_num_generator(seed_init);

  rands[tid] = rand_num_generator();
}

double
avg_gaussian_cost_per_iter(double ncall)
{
  // Gauss9D gauss_integrand;
  // constexpr int ndim = 9;
  return 0.;
}

template <typename IntegT, int ndim>
__global__ void
just_random_gen(IntegT* d_integrand,
                int npg,
                int chunkSize,
                uint32_t totalNumThreads,
                int LastChunk,
                unsigned int seed_init)
{
  uint32_t m = blockIdx.x * blockDim.x + threadIdx.x;

  if (m < totalNumThreads) {
    Random_num_generator<Curand_generator> rand_num_generator(seed_init);

    for (int t = 0; t < chunkSize; t++) {
      for (int k = 1; k <= npg; k++) {

        gpu::cudaArray<double, ndim> x;
        for (int i = 0; i < ndim; i++) {
          x[i] = rand_num_generator();
        }
        // gpu::apply(*d_integrand, x);
      } // end npg
    }   // end chunkSize
  }     // close if statement
}

template <typename IntegT, int ndim>
__global__ void
random_eval(IntegT* d_integrand,
            int npg,
            int chunkSize,
            uint32_t totalNumThreads,
            int LastChunk,
            unsigned int seed_init)
{

  uint32_t m = blockIdx.x * blockDim.x + threadIdx.x;

  if (m < totalNumThreads) {
    Random_num_generator<Curand_generator> rand_num_generator(seed_init);

    for (int t = 0; t < chunkSize; t++) {
      for (int k = 1; k <= npg; k++) {

        gpu::cudaArray<double, ndim> x;
        for (int i = 0; i < ndim; i++) {
          x[i] = rand_num_generator();
        }
        gpu::apply(*d_integrand, x);
      } // end npg
    }   // end chunkSize
  }     // close if statement
}

double
avg_sinsum_cost_per_iter()
{
  auto t0 = std::chrono::high_resolution_clock::now();
  SinSum6D integrand;
  SinSum6D* d_integrand = quad::cuda_copy_to_managed(integrand);
  constexpr int ndim = 6;
  double ncall = 2.0e9;
  // size_t BLOCK_DIM_X = 128;

  int chunkSize = GetChunkSize(ncall);
  double ncubes = num_sub_cubes(ncall, ndim);

  uint32_t totalNumThreads = (uint32_t)((ncubes) / chunkSize);
  uint32_t totalCubes = totalNumThreads * chunkSize;
  int extra = ncubes - totalCubes;
  int LastChunk = extra + chunkSize;
  uint32_t nBlocks =
    ((uint32_t)(((ncubes + BLOCK_DIM_X - 1) / BLOCK_DIM_X)) / chunkSize) + 1;
  uint32_t nThreads = BLOCK_DIM_X;
  int npg = samples_per_cube(ncall, ncubes);

  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
  unsigned int seed =
    static_cast<unsigned int>(time_diff.count()) + static_cast<unsigned int>(0);

  for (int i = 0; i < 10000; i++) {
    random_eval<SinSum6D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  for (int i = 0; i < 10000; i++) {
    just_random_gen<SinSum6D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  return 0.;
}

double
avg_gauss9D_cost_per_iter()
{
  auto t0 = std::chrono::high_resolution_clock::now();
  Gauss9D integrand;
  Gauss9D* d_integrand = quad::cuda_copy_to_managed(integrand);
  constexpr int ndim = 9;
  double ncall = 1.0e8;
  // size_t BLOCK_DIM_X = 128;

  int chunkSize = GetChunkSize(ncall);
  double ncubes = num_sub_cubes(ncall, ndim);

  uint32_t totalNumThreads = (uint32_t)((ncubes) / chunkSize);
  uint32_t totalCubes = totalNumThreads * chunkSize;
  int extra = ncubes - totalCubes;
  int LastChunk = extra + chunkSize;
  uint32_t nBlocks =
    ((uint32_t)(((ncubes + BLOCK_DIM_X - 1) / BLOCK_DIM_X)) / chunkSize) + 1;
  uint32_t nThreads = BLOCK_DIM_X;
  int npg = samples_per_cube(ncall, ncubes);

  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
  unsigned int seed =
    static_cast<unsigned int>(time_diff.count()) + static_cast<unsigned int>(0);

  for (int i = 0; i < 1000; i++) {
    random_eval<Gauss9D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  for (int i = 0; i < 1000; i++) {
    just_random_gen<Gauss9D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  return 0.;
}

double
avg_genz6_6D_cost_per_iter()
{
  auto t0 = std::chrono::high_resolution_clock::now();
  GENZ_6_6D integrand;
  GENZ_6_6D* d_integrand = quad::cuda_copy_to_managed(integrand);
  constexpr int ndim = 6;
  double ncall = 2.0e9;
  // size_t BLOCK_DIM_X = 128;

  int chunkSize = GetChunkSize(ncall);
  double ncubes = num_sub_cubes(ncall, ndim);

  uint32_t totalNumThreads = (uint32_t)((ncubes) / chunkSize);
  uint32_t totalCubes = totalNumThreads * chunkSize;
  int extra = ncubes - totalCubes;
  int LastChunk = extra + chunkSize;
  uint32_t nBlocks =
    ((uint32_t)(((ncubes + BLOCK_DIM_X - 1) / BLOCK_DIM_X)) / chunkSize) + 1;
  uint32_t nThreads = BLOCK_DIM_X;
  int npg = samples_per_cube(ncall, ncubes);

  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
  unsigned int seed =
    static_cast<unsigned int>(time_diff.count()) + static_cast<unsigned int>(0);

  for (int i = 0; i < 1000; i++) {
    random_eval<GENZ_6_6D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  for (int i = 0; i < 1000; i++) {
    just_random_gen<GENZ_6_6D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  return 0.;
}

double
avg_genz4_8D_cost_per_iter()
{
  auto t0 = std::chrono::high_resolution_clock::now();
  GENZ_4_8D integrand;
  GENZ_4_8D* d_integrand = quad::cuda_copy_to_managed(integrand);
  constexpr int ndim = 8;
  double ncall = 1.0e8;
  // size_t BLOCK_DIM_X = 128;

  int chunkSize = GetChunkSize(ncall);
  double ncubes = num_sub_cubes(ncall, ndim);

  uint32_t totalNumThreads = (uint32_t)((ncubes) / chunkSize);
  uint32_t totalCubes = totalNumThreads * chunkSize;
  int extra = ncubes - totalCubes;
  int LastChunk = extra + chunkSize;
  uint32_t nBlocks =
    ((uint32_t)(((ncubes + BLOCK_DIM_X - 1) / BLOCK_DIM_X)) / chunkSize) + 1;
  uint32_t nThreads = BLOCK_DIM_X;
  int npg = samples_per_cube(ncall, ncubes);

  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
  unsigned int seed =
    static_cast<unsigned int>(time_diff.count()) + static_cast<unsigned int>(0);

  for (int i = 0; i < 1000; i++) {
    random_eval<GENZ_4_8D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  for (int i = 0; i < 1000; i++) {
    just_random_gen<GENZ_4_8D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  return 0.;
}

double
avg_genz3_8D_cost_per_iter()
{
  auto t0 = std::chrono::high_resolution_clock::now();
  GENZ_3_8D integrand;
  GENZ_3_8D* d_integrand = quad::cuda_copy_to_managed(integrand);
  constexpr int ndim = 8;
  double ncall = 1.0e8;
  // size_t BLOCK_DIM_X = 128;

  int chunkSize = GetChunkSize(ncall);
  double ncubes = num_sub_cubes(ncall, ndim);

  uint32_t totalNumThreads = (uint32_t)((ncubes) / chunkSize);
  uint32_t totalCubes = totalNumThreads * chunkSize;
  int extra = ncubes - totalCubes;
  int LastChunk = extra + chunkSize;
  uint32_t nBlocks =
    ((uint32_t)(((ncubes + BLOCK_DIM_X - 1) / BLOCK_DIM_X)) / chunkSize) + 1;
  uint32_t nThreads = BLOCK_DIM_X;
  int npg = samples_per_cube(ncall, ncubes);

  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  MilliSeconds time_diff = std::chrono::high_resolution_clock::now() - t0;
  unsigned int seed =
    static_cast<unsigned int>(time_diff.count()) + static_cast<unsigned int>(0);

  for (int i = 0; i < 1000; i++) {
    random_eval<GENZ_3_8D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  for (int i = 0; i < 1000; i++) {
    just_random_gen<GENZ_3_8D, ndim><<<nBlocks, nThreads>>>(
      d_integrand, npg, chunkSize, totalNumThreads, LastChunk, seed);
    cudaDeviceSynchronize();
  }

  return 0.;
}

int
main()
{
  // avg_sinsum_cost_per_iter();
  // avg_gauss9D_cost_per_iter();
  // avg_genz6_6D_cost_per_iter();
  //   avg_genz3_8D_cost_per_iter();
  avg_genz4_8D_cost_per_iter();
  return 0;
}