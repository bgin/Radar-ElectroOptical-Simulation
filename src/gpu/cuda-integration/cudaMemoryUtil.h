#ifndef CUDACUHRE_QUAD_UTIL_CUDAMEMORY_UTIL_H
#define CUDACUHRE_QUAD_UTIL_CUDAMEMORY_UTIL_H

#include "cudaDebugUtil.h"
#include <cuda.h>

namespace quad {

  class Managed {
  public:
    void*
    operator new(size_t len)
    {
      void* ptr;
      cudaMallocManaged(&ptr, len);
      cudaDeviceSynchronize();
      return ptr;
    }

    void
    operator delete(void* ptr)
    {
      cudaDeviceSynchronize();
      cudaFree(ptr);
    }
  };

  template <typename T>
  class MemoryUtil {};

  template <typename T>
  class HostMemory : public MemoryUtil<T> {
  public:
    void*
    AllocateMemory(void* ptr, size_t n)
    {
      ptr = malloc(n);
      return ptr;
    }

    void
    ReleaseMemory(void* ptr)
    {
      free(ptr);
    }
  };

  template <typename T>
  class DeviceMemory : public MemoryUtil<T> {
  public:
    double
    GetFreeMemPercentage()
    {
      size_t free_physmem, total_physmem;
      QuadDebugExit(cudaMemGetInfo(&free_physmem, &total_physmem));
      return (double)free_physmem / total_physmem;
    }

    size_t
    GetAmountFreeMem()
    {
      size_t free_physmem, total_physmem;
      QuadDebugExit(cudaMemGetInfo(&free_physmem, &total_physmem));
      return free_physmem;
    }

    cudaError_t
    AllocateMemory(void** d_ptr, size_t n)
    {
      return cudaMalloc(d_ptr, n);
    }

    cudaError_t
    AllocateUnifiedMemory(void** d_ptr, size_t n)
    {
      return cudaMallocManaged(d_ptr, n);
    }

    cudaError_t
    ReleaseMemory(void* d_ptr)
    {
      return cudaFree(d_ptr);
    }

    cudaError_t
    SetHeapSize(size_t hSize = (size_t)2 * 1024 * 1024 * 1024)
    {
      return cudaDeviceSetLimit(cudaLimitMallocHeapSize, hSize);
    }

    cudaError_t
    CopyHostToDeviceConstantMemory(const char* d_ptr, void* h_ptr, size_t n)
    {
      return cudaMemcpyToSymbol(d_ptr, h_ptr, n);
    }

    cudaError_t
    CopyHostToDeviceConstantMemory(const void* d_ptr, void* h_ptr, size_t n)
    {
      return cudaMemcpyToSymbol(d_ptr, h_ptr, n);
    }

    //@brief Initialize Device
    cudaError_t
    DeviceInit(int dev = -1, int verbose = 0)
    {
      cudaError_t error = cudaSuccess;

      do {
        int deviceCount;
        error = QuadDebug(cudaGetDeviceCount(&deviceCount));
        if (error)
          break;
        if (deviceCount == 0) {
          fprintf(stderr, "No devices supporting CUDA.\n");
          exit(1);
        }

        if ((dev > deviceCount - 1) || (dev < 0)) {
          dev = 0;
        }

        // error = QuadDebug(cudaSetDevice(dev));
        if (error)
          break;

        size_t free_physmem, total_physmem;
        QuadDebugExit(cudaMemGetInfo(&free_physmem, &total_physmem));

        cudaDeviceProp deviceProp;
        error = QuadDebug(cudaGetDeviceProperties(&deviceProp, dev));
        if (error)
          break;

        if (deviceProp.major < 1) {
          fprintf(stderr, "Device does not support CUDA.\n");
          exit(1);
        }
        if (false && verbose) {
          printf("Using device %d: %s (SM%d, %d SMs, %lld free / %lld total MB "
                 "physmem, ECC %s)\n",
                 dev,
                 deviceProp.name,
                 deviceProp.major * 100 + deviceProp.minor * 10,
                 deviceProp.multiProcessorCount,
                 (unsigned long long)free_physmem / 1024 / 1024,
                 (unsigned long long)total_physmem / 1024 / 1024,
                 (deviceProp.ECCEnabled) ? "on" : "off");
          fflush(stdout);
        }

      } while (0);
      return error;
    }
  };

  template <class T>
  T*
  cuda_malloc_managed(size_t size)
  {
    CudaCheckError();
    T* temp = nullptr;
    auto rc = cudaMallocManaged(&temp, sizeof(T) * size);
    if (rc != cudaSuccess) {
      temp = nullptr;

      size_t free_physmem, total_physmem;
      cudaMemGetInfo(&free_physmem, &total_physmem);
      printf("cuda_malloc_managed(size) allocating size %lu free mem:%lu\n",
             size,
             free_physmem);

      CudaCheckError();
      throw std::bad_alloc();
    }
    return temp;
  }

  template <class T>
  T*
  cuda_malloc_managed()
  {
    T* temp = nullptr;
    auto rc = cudaMallocManaged(&temp, sizeof(T));
    if (rc != cudaSuccess) {
      size_t free_physmem, total_physmem;
      cudaMemGetInfo(&free_physmem, &total_physmem);
      printf("cuda_malloc_managed() allocating size %lu free mem:%lu\n",
             sizeof(T),
             free_physmem);
      CudaCheckError();
      throw std::bad_alloc();
    }
    CudaCheckError();
    return temp;
  }

  template <class T>
  T*
  cuda_copy_to_managed(T const& on_host)
  {
    T* buffer = cuda_malloc_managed<T>();
    CudaCheckError();
    try {
      new (buffer) T(on_host);
      CudaCheckError();
    }
    catch (...) {
      cudaFree(buffer);
      throw;
    }
    return buffer;
  }

  template <class T>
  T*
  cuda_malloc(size_t size)
  {
    T* temp;
    auto rc = cudaMalloc((void**)&temp, sizeof(T) * size);
    if (rc != cudaSuccess) {
      throw std::bad_alloc();
    }
    return temp;
  }

  template <typename T>
  void
  cuda_memcpy_to_device(T* dest, T* src, size_t size)
  {
    auto rc = cudaMemcpy(dest, src, sizeof(T) * size, cudaMemcpyHostToDevice);
    if (rc != cudaSuccess) {
      printf("error in cuda_mempcy_to_device with host src\n");
      throw std::bad_alloc();
      abort();
    }
  }

  template <typename T>
  void
  cuda_memcpy_to_device(T* dest, const T* src, size_t size)
  {
    auto rc = cudaMemcpy(dest, src, sizeof(T) * size, cudaMemcpyHostToDevice);
    if (rc != cudaSuccess) {
      printf("error in cuda_mempcy_to_device with host src\n");
      throw std::bad_alloc();
      abort();
    }
  }

  template <typename T>
  void
  cuda_memcpy_device_to_device(T* dest, T* src, size_t size)
  {
    auto rc = cudaMemcpy(dest, src, sizeof(T) * size, cudaMemcpyDeviceToDevice);
    if (rc != cudaSuccess) {
      printf("error in cuda_memcpy_device_to_device\n");
      throw std::bad_alloc();
      abort();
    }
  }

  template <class T>
  T*
  cuda_copy_to_device(T const& on_host)
  {
    T* buffer = cuda_malloc<T>(1);
    const T* hp = &on_host;
    cuda_memcpy_to_device<T>(buffer, hp, 1);
    CudaCheckError();
    return buffer;
  }

  inline size_t
  get_free_mem()
  {
    size_t free_physmem, total_physmem;
    QuadDebugExit(cudaMemGetInfo(&free_physmem, &total_physmem));
    return free_physmem;
  }
  
  
template <typename T>
T*
copy_to_host(T* src, size_t size)
{
  T* dest = new T[size];
  auto rc = cudaMemcpy(dest, src, sizeof(T) * size, cudaMemcpyDeviceToHost);
  if (rc != cudaSuccess)
    throw std::bad_alloc();
  return dest;
}

template <typename T>
void
cuda_memcpy_to_host(T* dest, T const* src, size_t n_elements)
{
  auto rc =
    cudaMemcpy(dest, src, sizeof(T) * n_elements, cudaMemcpyDeviceToHost);
  if (rc != cudaSuccess)
    throw std::bad_alloc();
}


template <typename T>
void
cuda_memcpy_device_to_device(T* dest, T const* src, size_t size)
{
  auto rc = cudaMemcpy(dest, src, sizeof(T) * size, cudaMemcpyDeviceToDevice);
  if (rc != cudaSuccess)
    throw std::bad_alloc();
}

template <typename T>
struct Range {
  Range() = default;
  Range(T l, T h) : low(l), high(h) {}
  T low = 0., high = 0.;
};

template <typename T>
__global__ void
device_print_array(T* arr, size_t size)
{
  for (size_t i = 0; i < size; ++i)
    printf(
      "arr[%lu]:%i\n", i, arr[i]); // can't print arbitrary types from device,
                                   // must fix to do std::cout from host
}

template <typename T>
void
print_device_array(T* arr, size_t size)
{
  device_print_array<T><<<1, 1>>>(arr, size);
  cudaDeviceSynchronize();
}

template <class T>
T*
host_alloc(size_t size)
{
  T* temp = new T[size];
  if (temp == nullptr) {
	printf("cannot allocate on host\n");
    throw std::bad_alloc();
  }
  return temp;
}

// rename to  free_and_reallocate, delete copy_size, unecessary
template <typename T>
void
ExpandcuArray(T*& array, int currentSize, int newSize)
{
  int copy_size = std::min(currentSize, newSize);
  T* temp = cuda_malloc<T>(newSize);
  cudaFree(array);
  array = temp;
}

template <typename IntegT>
IntegT*
make_gpu_integrand(const IntegT& integrand)
{
  IntegT* d_integrand;
  cudaMallocManaged((void**)&d_integrand, sizeof(IntegT));
  memcpy(d_integrand, &integrand, sizeof(IntegT));
  return d_integrand;
}

template <typename T>
__global__ void
set_array_to_value(T* array, size_t size, T val)
{
  size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < size) {
    array[tid] = val;
  }
}

template <typename T>
__global__ void
set_array_range_to_value(T* array,
                         size_t first_to_change,
                         size_t last_to_change,
                         size_t total_size,
                         T val)
{
  size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= first_to_change && tid <= last_to_change && tid < total_size) {
    array[tid] = val;
  }
}

template <typename T>
void
set_device_array_range(T* arr,
                       size_t first_to_change,
                       size_t last_to_change,
                       size_t size,
                       T val)
{
  size_t num_threads = 64;
  size_t num_blocks = size / num_threads + ((size % num_threads) ? 1 : 0);
  set_array_range_to_value<T><<<num_blocks, num_threads>>>(
    arr, first_to_change, last_to_change, size, val);
  cudaDeviceSynchronize();
}

template <typename T>
void
set_device_array(T* arr, size_t size, T val)
{
  size_t num_threads = 64;
  size_t num_blocks = size / num_threads + ((size % num_threads) ? 1 : 0);
  set_array_to_value<T><<<num_blocks, num_threads>>>(arr, size, val);
  cudaDeviceSynchronize();
}

template <typename T, typename C = T>
bool
array_values_smaller_than_val(T* dev_arr, size_t dev_arr_size, C val)
{
  double* host_arr = host_alloc<double>(dev_arr_size);
  cudaMemcpy(
    host_arr, dev_arr, sizeof(double) * dev_arr_size, cudaMemcpyDeviceToHost);

  for (int i = 0; i < dev_arr_size; i++) {
    if (host_arr[i] >= static_cast<T>(val))
      return false;
  }
  return true;
}

template <typename T, typename C = T>
bool
array_values_larger_than_val(T* dev_arr, size_t dev_arr_size, C val)
{
  double* host_arr = host_alloc<double>(dev_arr_size);
  cudaMemcpy(
    host_arr, dev_arr, sizeof(double) * dev_arr_size, cudaMemcpyDeviceToHost);

  for (int i = 0; i < dev_arr_size; i++) {
    if (host_arr[i] < static_cast<T>(val)) {
      std::cout << "host_arr[" << i << "]:" << host_arr[i] << " val:" << val
                << "\n";
      return false;
    }
  }
  return true;
}


}

#endif
