#include "cuda/mcubes/demos/demo_utils.cuh"
#include "cuda/mcubes/vegasT.cuh"

class GENZ_1_8D {
public:
  __device__ __host__ double
  operator()(double s,
             double t,
             double u,
             double v,
             double w,
             double x,
             double y,
             double z)
  {
    return cos(s + 2. * t + 3. * u + 4. * v + 5. * w + 6. * x + 7. * y +
               8. * z);
  }
};
int
main(int argc, char** argv)
{
  double epsrel = 1e-3;
  double epsrel_min = 1e-9;
  constexpr int ndim = 8;

  double ncall = 1.0e9;
  int titer = 100;
  int itmax = 50;
  int skip = 5;
  VegasParams params(ncall, titer, itmax, skip);

  double true_value = (1. / 315.) * sin(1.) * sin(3. / 2.) * sin(2.) *
                      sin(5. / 2.) * sin(3.) * sin(7. / 2.) * sin(4.) *
                      (sin(37. / 2.) - sin(35. / 2.));

  double lows[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double highs[] = {1., 1., 1., 1., 1., 1., 1., 1.};
  quad::Volume<double, ndim> volume(lows, highs);
  GENZ_1_8D integrand;

  /*std::map<double, double> required_ncall {
    {1.e-3, 6.e9}, {2.e-4, 6.e9}, {4.e-5, 6.e9}, {8.e-6, 6.e9}, {1.6e-6, 6.e9}
  };*/

  print_mcubes_header();
  bool success = false;
  do {
    params.ncall = ncall; // required_ncall[epsrel];
    success = mcubes_time_and_call<GENZ_1_8D, ndim>(
      integrand, epsrel, true_value, "Genz1_8D", params, &volume);
    epsrel /= 5.;
  } while (epsrel >= epsrel_min && success == true);

  return 0;
}
