#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"

class GENZ_3_6D {
public:
  __device__ __host__ double
  operator()(double x, double y, double z, double w, double v, double u)
  {
    return pow(1 + 6 * u + 5 * v + 4 * w + 3 * x + 2 * y + z, -7);
  }
};

int
main(int argc, char** argv)
{
  double epsrel = 1e-3;
  double epsrel_min = 1.e-6;
  constexpr int ndim = 6;

  double ncall = 1.0e7;
  int titer = 100;
  int itmax = 20;
  int skip = 5;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = 7.1790160638199853886e-7;

  double lows[] = {0., 0., 0., 0., 0., 0.};
  double highs[] = {1., 1., 1., 1., 1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_3_6D integrand;

  print_mcubes_header();
  while (mcubes_time_and_call<GENZ_3_6D, ndim>(
           integrand, epsrel, true_value, "GENZ_3_6D", params, &volume) ==
           true &&
         epsrel >= epsrel_min) {
    epsrel /= 5.;
  }

  return 0;
}