#include "cuda/mcubes/mcubesSeq.hh"
#include <iostream>
#include <tuple>
#include <utility>

class GENZ_3_3D {
public:
  __device__ __host__ double
  operator()(double x, double y, double z)
  {
    return pow(1 + 3 * x + 2 * y + z, -4);
  }
};

template <typename T, int NDIM>
struct Volume {

  T lows[NDIM] = {0.0};
  T highs[NDIM];

  __host__
  Volume()
  {
    for (T& x : highs)
      x = 1.0;
  }

  __host__
  Volume(std::array<T, NDIM> l, std::array<T, NDIM> h)
  {
    std::memcpy(lows, l.data(), NDIM * sizeof(T));
    std::memcpy(highs, h.data(), NDIM * sizeof(T));
  }

  __host__ __device__
  Volume(T const* l, T const* h)
  {
    std::memcpy(lows, l, NDIM * sizeof(T));
    std::memcpy(highs, h, NDIM * sizeof(T));
  }
};

int
main()
{
  double epsrel = 1e-3;

  double ncall = 1.0e7;
  int titer = 20;
  constexpr int ndim = 3;

  double lows[] = {0., 0., 0.};
  double highs[] = {1., 1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_3_3D integrand;
  auto res = seq_mcubes_integrate<GENZ_3_3D, ndim>(
    integrand, ndim, epsrel, 1.e-12, ncall, &volume, titer);

  std::cout.precision(17);
  // std::cout<<res.estimate<< "," res.errorest
  // std::cout<<res<<"\n";
  return 0;
}
