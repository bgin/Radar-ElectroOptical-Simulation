#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"

class GENZ_3_3D {
public:
  __device__ __host__ double
  operator()(double x, double y, double z)
  {
    return pow(1 + 3 * x + 2 * y + z, -4);
  }
};

int
main(int argc, char** argv)
{
  double epsrel = 1.e-3;
  double epsrel_min = 1e-9;
  constexpr int ndim = 3;

  double ncall = 1.e7;
  int titer = 100;
  int itmax = 20;
  int skip = 0;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = 0.010846560846560846561;
  double lows[] = {0., 0., 0.};
  double highs[] = {1., 1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_3_3D integrand;

  print_mcubes_header();
  std::array<double, 7> required_ncall =
    //{1.e7, 1.e7, 1.e7, 1.e7, 1.e7, 1.e8, 1.e9};
    {1.e7, 1.e7, 1.e7, 1.e7, 1.e7, 1.e8, 1.e9};
  bool success = false;
  size_t curr_epsrel = 0;

  do {
    params.ncall = required_ncall[curr_epsrel];
    for (int run = 0; run < 1; run++) {
      success = mcubes_time_and_call<GENZ_3_3D, ndim, false, Custom_generator>(
        integrand, epsrel, true_value, "f3, 3", params, &volume);
      if (!success)
        break;
    }
    break;
    epsrel /= 5.;
    curr_epsrel++;
    if (curr_epsrel > required_ncall.size())
      break;
  } while (epsrel >= epsrel_min && success == true);
  return 0;
}
