#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"

class GENZ_5_8D {
public:
  __device__ __host__ double
  operator()(double x,
             double y,
             double z,
             double k,
             double m,
             double n,
             double p,
             double q)
  {
    double beta = .5;
    double t1 = -10. * fabs(x - beta) - 10. * fabs(y - beta) -
                10. * fabs(z - beta) - 10. * fabs(k - beta) -
                10. * fabs(m - beta) - 10. * fabs(n - beta) -
                10. * fabs(p - beta) - 10. * fabs(q - beta);
    return exp(t1);
  }
};

int
main(int argc, char** argv)
{
  double epsrel = 1.e-3;
  double epsrel_min = 1e-9;
  constexpr int ndim = 8;

  double ncall = 1.0e6;
  int titer = 100;
  int itmax = 20;
  int skip = 5;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = 2.425217625641885e-06;

  double lows[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double highs[] = {1., 1., 1., 1., 1., 1., 1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_5_8D integrand;

  print_mcubes_header();
  std::array<double, 10> required_ncall = {
    1.e6, 1.e6, 1.e6, 1.e7, 1.e9, 1.e9, 5.e9, 8.e9, 8.e9, 8.e9};

  bool success = false;
  size_t curr_epsrel = 0;
  do {
    params.ncall = required_ncall[curr_epsrel];
    for (int run = 0; run < 100; run++) {
      success = mcubes_time_and_call<GENZ_5_8D, ndim, false, Custom_generator>(
        integrand, epsrel, true_value, "f5, 8", params, &volume);
      if (!success)
        break;
    }
    epsrel /= 5.;
    curr_epsrel++;
    if (curr_epsrel > required_ncall.size())
      break;
  } while (epsrel >= epsrel_min && success == true);
  return 0;
}
