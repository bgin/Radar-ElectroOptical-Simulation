#ifndef CUDACUHRE_QUAD_QUAD_h
#define CUDACUHRE_QUAD_QUAD_h

#define TIMING_DEBUG 1
#define BLOCK_SIZE 64
#define SM_REGION_POOL_SIZE 128

#define GLOBAL_ERROR 1
#define MAX_GLOBALPOOL_SIZE 2048

#include <fstream>
#include <string>
#include <vector>
#include "integration_result.hh"

using TYPE = double;

static int FIRST_PHASE_MAXREGIONS = (1 << 14);

// Utilities
#include "cudaArchUtil.h"
#include "cudaDebugUtil.h"

class VerboseResults {
public:
  std::vector<std::vector<double>> funcEvaluationPoints;
  std::vector<double> results;
  size_t numFuncEvals = 0;
  size_t NDIM = 0;
};

template <typename T>
struct Structures {
  T* gpuG = nullptr;
  T* cRuleWt = nullptr;
  T* GPUScale = nullptr;
  T* GPUNorm = nullptr;
  int* gpuGenPos = nullptr;
  int* gpuGenPermGIndex = nullptr;
  int* gpuGenPermVarCount = nullptr;
  int* gpuGenPermVarStart = nullptr;
  size_t* cGeneratorCount = nullptr;
};

struct Result {
  double avg = 0., err = 0.;
  int bisectdim = 0;
};

struct Bounds {
  double lower, upper;
};

struct GlobalBounds {
  double unScaledLower, unScaledUpper;
};

template <int dim>
struct Region {
  int div;
  Result result;
  Bounds bounds[dim];
};

#define NRULES 5

namespace pagani {
  template <size_t ndim>
  __host__ __device__ constexpr size_t
  CuhreFuncEvalsPerRegion()
  {
    return (1 + 2 * ndim + 2 * ndim + 2 * ndim + 2 * ndim +
            2 * ndim * (ndim - 1) + 4 * ndim * (ndim - 1) +
            4 * ndim * (ndim - 1) * (ndim - 2) / 3 + (1 << ndim));
  }
}

template <int NDIM>
struct Snapshot {
  __host__
  Snapshot(int* iterations, int size)
  {
    numSnapshots = size;
    total_regions = 0;
    currArrayHead = 0;

    for (int i = 0; i < size; i++)
      total_regions += iterations[i];

    cudaMalloc((void**)&arr, sizeof(Region<NDIM>) * total_regions);
    cudaMalloc((void**)&sizes, sizeof(int) * numSnapshots);
    cudaMemcpy(
      sizes, iterations, sizeof(int) * numSnapshots, cudaMemcpyHostToDevice);
  }

  __host__ void
  Save(std::string baseFileName)
  {

    Region<NDIM>* h_arr = 0;
    int* h_sizes = 0;

    h_arr = (Region<NDIM>*)malloc(sizeof(Region<NDIM>) * total_regions);
    h_sizes = (int*)malloc(sizeof(int) * numSnapshots);

    cudaMemcpy(
      h_arr, arr, sizeof(Region<NDIM>) * total_regions, cudaMemcpyDeviceToHost);
    cudaMemcpy(
      h_sizes, sizes, sizeof(int) * numSnapshots, cudaMemcpyDeviceToHost);
    int index = 0;

    for (int i = 0; i < numSnapshots; i++) {
      std::string filename = baseFileName + std::to_string(h_sizes[i]) + ".csv";
      std::ofstream outfile(filename.c_str());
      int current_size = h_sizes[i];
      int snapShotStartIndex = 0;

      for (int j = 0; j < i; j++)
        snapShotStartIndex += h_sizes[j];

      for (; index < current_size + snapShotStartIndex; index++) {
        outfile << h_arr[index].result.avg << "," << h_arr[index].result.err
                << ",";
        for (int dim = 0; dim < NDIM; dim++) {
          outfile << h_arr[index].bounds[dim].upper << ","
                  << h_arr[index].bounds[dim].lower << ",";
        }
        outfile << h_arr[index].div << "," << -1 << std::endl;
      }
      outfile.close();
    }

    free(h_arr);
    free(h_sizes);
    cudaFree(sizes);
    cudaFree(arr);
  }

  __host__ __device__
  Snapshot()
  {
    numSnapshots = 0;
    arr = nullptr;
    sizes = nullptr;
  }

  int currArrayHead;
  int numSnapshots;
  Region<NDIM>* arr;
  int* sizes;
  int total_regions;
};

#endif
