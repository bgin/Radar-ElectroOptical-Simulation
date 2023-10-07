#ifndef CUDACUHRE_QUAD_UTIL_CUDAAPPLY_CUH
#define CUDACUHRE_QUAD_UTIL_CUDAAPPLY_CUH

#include "cudaArray.cuh"
namespace gpu {

  namespace detail {
    template <class F, typename T, size_t N, std::size_t... I>
    __device__ T
    apply_impl(F&& f,
               gpu::cudaArray<T, N> const& data,
               std::index_sequence<I...>)
    {
      return f(data[I]...);
    };
  }

  template <class F, typename T, size_t N>
  __device__ T
  // Unsure if we need to pass 'f' by value, for GPU execution
  apply(F&& f, gpu::cudaArray<T, N> const& data)
  {
    return detail::apply_impl(
      std::forward<F>(f), data, std::make_index_sequence<N>());
  }
}

#endif
