#ifndef CUDACUHRE_QUAD_UTIL_CUDAARCH_UTIL_H
#define CUDACUHRE_QUAD_UTIL_CUDAARCH_UTIL_H
namespace quad {

  /// QUAD_PTX_ARCH reflects the PTX version targeted by the active compiler
  /// pass (or zero during the host pass).
#ifndef __CUDA_ARCH__
#define QUAD_PTX_ARCH 0
#else
#define QUAD_PTX_ARCH __CUDA_ARCH__
#endif
}

#endif
