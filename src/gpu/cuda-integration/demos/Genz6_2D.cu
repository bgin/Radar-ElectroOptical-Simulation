#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"

class GENZ_6_2D {
public:
  __device__ __host__ double
  operator()(double y, double z)
  {
    if (z > .9 || y > .8)
      return 0.;
    else
      return exp(10 * z + 9 * y);
  }
};

int
main(int argc, char** argv)
{
  double epsrel = 1e-3;
  double epsrel_min = 1.e-6;
  constexpr int ndim = 2;

  double ncall = 1.0e7;
  int titer = 100;
  int itmax = 20;
  int skip = 5;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = 120489.75982636053604;

  double lows[] = {0., 0.};
  double highs[] = {1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_6_2D integrand;

  print_mcubes_header();
  while (mcubes_time_and_call<GENZ_6_2D, ndim>(
           integrand, epsrel, true_value, "GENZ_6_2D", params, &volume) ==
           true &&
         epsrel >= epsrel_min) {
    epsrel /= 5.;
  }

  return 0;
}