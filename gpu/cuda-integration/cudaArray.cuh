#ifndef CUDACUHRE_QUAD_UTIL_CUDAARRAY_CUH
#define CUDACUHRE_QUAD_UTIL_CUDAARRAY_CUH

#include <cstring>
#include <array>
#include "quad.h"
#include "cudaMemoryUtil.h"

// cudaArray is meant to allow in-kernel use of functions that expect std::array interface, e.g. std::forward

namespace gpu {
  template <typename T, std::size_t s>
  class cudaArray {
  public:
    void
    Initialize(T const* initData)
    {
      std::memcpy(data, initData, sizeof(T) * s);
    }

    __host__ __device__ const T*
    begin() const
    {
      return &data[0];
    }

    __host__ __device__ const T*
    end() const
    {
      return (&data[0] + s);
    }

    __host__ __device__ constexpr std::size_t
    size() const
    {
      return s;
    }

    __host__ __device__ T&
    operator[](std::size_t i)
    {
      return data[i];
    }

    __host__ __device__ T const&
    operator[](std::size_t i) const
    {
      return data[i];
    }

    T data[s];
  };

  template <typename T>
  class cudaDynamicArray {
  public:
    __host__ __device__
    cudaDynamicArray(const cudaDynamicArray& a)
    {
#ifndef __CUDA_ARCH__
      N = a.N;
	  _data = quad::cuda_malloc_managed<T>(N);
      memcpy(_data, a._data, sizeof(T) * a.N);
#else
      // can't instantiate on device and then access on host
      N = a.N;
      _data = new T[a.N];
      memcpy(_data, a._data, sizeof(T) * a.N);
#endif
    }

    __host__ __device__
    cudaDynamicArray()
    {
      _data = nullptr;
      N = 0;
    }

    // make everything host device

    cudaDynamicArray(T const* initData, size_t s)
    {
      N = s;
	  _data = quad::cuda_malloc_managed<T>(s);
	  quad::cuda_memcpy_to_device<T>(_data, initData, s);
    }

	explicit cudaDynamicArray(size_t s)
    {
      N = s;
	  _data = quad::cuda_malloc_managed<T>(s);
    }

    void
    Reserve(size_t s)
    {
      N = s;
	  _data = quad::cuda_malloc_managed<T>(s);
    }

    __host__ __device__ ~cudaDynamicArray()
    {
#ifndef __CUDACC__
      cudaFree(_data);
#endif
    }

    __host__ __device__ const T*
    begin() const
    {
      return &_data[0];
    }

    __host__ __device__ const T*
    end() const
    {
      return (&_data[0] + N);
    }

    __host__ __device__ constexpr std::size_t
    size() const
    {
      return N;
    }

	__host__ __device__ 
	T*
	data(){
		return _data;
	}
	
    __host__ __device__ T&
    operator[](std::size_t i)
    {
      return _data[i];
    }

    __host__ __device__ T
    operator[](std::size_t i) const
    {
      return _data[i];
    }
	
	
	private: 
	
		T* _data;
		size_t N;
  }; // cudaDynamicArray

};

#endif
