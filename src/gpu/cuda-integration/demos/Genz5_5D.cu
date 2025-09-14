#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"

class GENZ_5_5D {
public:
  __device__ __host__ double
  operator()(double x, double y, double z, double k, double m)
  {
    double beta = .5;
    double t1 = -10. * fabs(x - beta) - 10. * fabs(y - beta) -
                10. * fabs(z - beta) - 10. * fabs(k - beta) -
                10. * fabs(m - beta);
    return exp(t1);
  }
};

int
main(int argc, char** argv)
{
  double epsrel = 1.e-3;
  double epsrel_min = 1e-9;
  constexpr int ndim = 5;

  double ncall = 1.0e6;
  int titer = 100;
  int itmax = 20;
  int skip = 5;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = 0.0003093636;

  double lows[] = {0., 0., 0., 0., 0.};
  double highs[] = {1., 1., 1., 1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_5_5D integrand;

  print_mcubes_header();
  // std::array<double, 10> required_ncall =
  // {1.e7, 1.e7, 1.e7, 1.e7, 1.e7, 1.e9, 3.e9, 8.e9, 8.e9, 8.e9};

  bool success = false;
  do {
    params.ncall = ncall; // required_ncall[expID];
    success = mcubes_time_and_call<GENZ_5_5D, ndim>(
      integrand, epsrel, true_value, "GENZ_5_5D", params, &volume);
    epsrel /= 5.;
  } while (epsrel >= epsrel_min && success == true);

  return 0;
}