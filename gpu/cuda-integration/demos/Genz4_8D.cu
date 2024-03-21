#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"

class GENZ_4_8D {
public:
  __device__ __host__ double
  operator()(double x,
             double y,
             double z,
             double w,
             double v,
             double k,
             double m,
             double n)
  {
    // double alpha = 25.;
    double beta = .5;
    return exp(-1.0 *
               (pow(25, 2) * pow(x - beta, 2) + pow(25, 2) * pow(y - beta, 2) +
                pow(25, 2) * pow(z - beta, 2) + pow(25, 2) * pow(w - beta, 2) +
                pow(25, 2) * pow(v - beta, 2) + pow(25, 2) * pow(k - beta, 2) +
                pow(25, 2) * pow(m - beta, 2) + pow(25, 2) * pow(n - beta, 2)));
  }
};

int
main(int argc, char** argv)
{
  double epsrel = 1e-3;
  double epsrel_min = 1e-9;
  constexpr int ndim = 8;

  double ncall = 1.0e6;
  int titer = 100;
  int itmax = 20;
  int skip = 5;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = (6.383802190004379e-10);

  double lows[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double highs[] = {1., 1., 1., 1., 1., 1., 1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_4_8D integrand;

  print_mcubes_header();
  //  std::array<double, 6> required_ncall =
  //  {1.e7, 1.e7, 1.e7, 1.e9, 2.e9, 6.e9};

  bool success = false;
  // size_t expID = 0;
  do {
    // params.ncall = ncall;//required_ncall[expID];
    for (int run = 0; run < 100; run++) {
      success = mcubes_time_and_call<GENZ_4_8D, ndim>(
        integrand, epsrel, true_value, "f4 8D", params, &volume);
      if (!success)
        break;
    }
    epsrel /= 5.;
    // expID++;
  } while (epsrel >= epsrel_min && success == true);
  return 0;
}
