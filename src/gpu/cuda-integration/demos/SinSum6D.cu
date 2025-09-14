#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"
#include <chrono>
#include <iostream>
#include <string>

class SinSum6D {
public:
  __device__ __host__ double
  operator()(double x, double y, double z, double k, double l, double m)
  {
    return sin(x + y + z + k + l + m);
  }
};

int
main(int argc, char** argv)
{
  double epsrel = 1.e-3;
  double epsabs = 1.e-20;

  constexpr int ndim = 6;
  double ncall = 2.0e9;
  int titer = 10;
  int itmax = 0;
  int skip = 0;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = -49.165073;
  std::cout << "id, estimate, std, chi, iters, adj_iters, skip_iters, ncall, "
               "time, abserr, relerr\n";

  double lows[] = {0., 0., 0., 0., 0., 0.};
  double highs[] = {10., 10., 10., 10., 10., 10.};
  quad::Volume<double, ndim> volume(lows, highs);
  SinSum6D integrand;
  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;

  constexpr bool MCUBES_DEBUG = false;
  auto t0 = std::chrono::high_resolution_clock::now();
  auto res = cuda_mcubes::integrate<SinSum6D, ndim, MCUBES_DEBUG>(
    integrand,
    epsrel,
    epsabs,
    params.ncall,
    &volume,
    params.t_iter,
    params.num_adjust_iters,
    params.num_skip_iters);
  MilliSeconds dt = std::chrono::high_resolution_clock::now() - t0;

  std::cout.precision(15);
  std::cout << "SinSum6D"
            << "," << epsrel << "," << std::scientific << true_value << ","
            << std::scientific << res.estimate << "," << std::scientific
            << res.errorest << "," << res.chi_sq << "," << params.t_iter << ","
            << params.num_adjust_iters << "," << params.num_skip_iters << ","
            << res.iters << "," << params.ncall << "," << res.neval << ","
            << dt.count() << "," << res.status << "\n";

  return 0;
}