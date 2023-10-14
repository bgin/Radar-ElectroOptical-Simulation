#ifndef VEGAS_DEMO_UTILS_CUH
#define VEGAS_DEMO_UTILS_CUH

#include "vegasT.cuh"
#include "vegasT1D.cuh"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

using std::cout;
using std::chrono::duration;
using std::chrono::high_resolution_clock;

bool
ApproxEqual(double a, double b, double epsilon = 1.e-5)
{
  if (std::abs(a - b) <= epsilon)
    return true;
  return false;
}

struct VegasParams {

  VegasParams(double callsPerIter,
              double total_iters,
              double adjust_iters,
              int skipIters)
    : ncall(callsPerIter)
    , t_iter(total_iters)
    , num_adjust_iters(adjust_iters)
    , num_skip_iters(skipIters){};

  double ncall = 1.e7;
  int t_iter = 70;
  int num_adjust_iters = 40;
  int num_skip_iters = 5;
};

void
print_mcubes_header()
{
  std::cout << "id, epsrel, integral, estimate, errorest, chi, iters, "
               "adj_iters, skip_iters, completed_iters, ncall, neval,"
               "time, status\n";
}

template <typename F,
          int ndim,
          bool MCUBES_DEBUG = false,
          typename GeneratorType = Curand_generator>
bool
mcubes_time_and_call(F integrand,
                     double epsrel,
                     double correct_answer,
                     char const* integralName,
                     VegasParams& params,
                     quad::Volume<double, ndim>* volume)
{
  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  // We make epsabs so small that epsrel is always the stopping condition.
  double constexpr epsabs = 1.0e-20;
  bool success = false;
  int run = 0;

  do {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto res = cuda_mcubes::integrate<F, ndim, MCUBES_DEBUG, GeneratorType>(
      integrand,
      epsrel,
      epsabs,
      params.ncall,
      volume,
      params.t_iter,
      params.num_adjust_iters,
      params.num_skip_iters);
    MilliSeconds dt = std::chrono::high_resolution_clock::now() - t0;
    success = (res.status == 0);
    std::cout.precision(15);

    if (success)
      std::cout << integralName << "," << epsrel << "," << std::scientific
                << correct_answer << "," << std::scientific << res.estimate
                << "," << std::scientific << res.errorest << "," << res.chi_sq
                << "," << params.t_iter << "," << params.num_adjust_iters << ","
                << params.num_skip_iters << "," << res.iters << ","
                << params.ncall << "," << res.neval << "," << dt.count() << ","
                << res.status << "\n";

    if (run == 0 && !success)
      AdjustParams(params.ncall, params.t_iter);
    if (success)
      run++;
  } while (success == false &&
           CanAdjustNcallOrIters(params.ncall, params.t_iter) == true);

  return success;
}

template <typename F,
          int ndim,
          bool MCUBES_DEBUG = false,
          typename GeneratorType = Curand_generator>
bool
mcubes1D_time_and_call(F integrand,
                       double epsrel,
                       double correct_answer,
                       char const* integralName,
                       VegasParams& params,
                       quad::Volume<double, ndim>* volume)
{
  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  // We make epsabs so small that epsrel is always the stopping condition.
  double constexpr epsabs = 1.0e-20;
  bool success = false;
  int run = 0;
  do {

    double exp_epsrel = epsrel /**.5*/;
    // std::cout<<"Trying with ncall:"<<ncall<<"\n";
    auto t0 = std::chrono::high_resolution_clock::now();
    auto res = cuda_mcubes::integrate<F, ndim, MCUBES_DEBUG, GeneratorType>(
      integrand,
      exp_epsrel,
      epsabs,
      params.ncall,
      volume,
      params.t_iter,
      params.num_adjust_iters,
      params.num_skip_iters);
    MilliSeconds dt = std::chrono::high_resolution_clock::now() - t0;
    success = (res.status == 0);
    std::cout.precision(15);

    if (success)
      std::cout << integralName << "," << epsrel << "," << std::scientific
                << correct_answer << "," << std::scientific << res.estimate
                << "," << std::scientific << res.errorest << "," << res.chi_sq
                << "," << params.t_iter << "," << params.num_adjust_iters << ","
                << params.num_skip_iters << "," << res.iters << ","
                << params.ncall << "," << res.neval << "," << dt.count() << ","
                << res.status << "\n";

    if (run == 0 && !success)
      AdjustParams(params.ncall, params.t_iter);
    if (success)
      run++;
  } while (success == false &&
           CanAdjustNcallOrIters(params.ncall, params.t_iter) == true);

  return success;
}

template <typename F,
          int ndim,
          bool MCUBES_DEBUG = false,
          typename GeneratorType = Curand_generator>
bool
common_header_mcubes_time_and_call(F integrand,
                                   double epsrel,
                                   double correct_answer,
                                   double difficulty,
                                   std::string alg_id,
                                   std::string integ_id,
                                   VegasParams& params,
                                   quad::Volume<double, ndim>* volume)
{
  using MilliSeconds =
    std::chrono::duration<double, std::chrono::milliseconds::period>;
  // We make epsabs so small that epsrel is always the stopping condition.
  double constexpr epsabs = 1.0e-20;
  bool success = false;
  int run = 0;

  do {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto res = cuda_mcubes::integrate<F, ndim, MCUBES_DEBUG, GeneratorType>(
      integrand,
      epsrel,
      epsabs,
      params.ncall,
      volume,
      params.t_iter,
      params.num_adjust_iters,
      params.num_skip_iters);
    MilliSeconds dt = std::chrono::high_resolution_clock::now() - t0;
    success = (res.status == 0);
    std::cout.precision(17);

    if (success)
      std::cout << integ_id << "," << std::scientific << alg_id << ","
                << difficulty << "," << epsrel << "," << epsabs << ","
                << std::scientific << correct_answer << "," << std::scientific
                << res.estimate << "," << std::scientific << res.errorest << ","
                << dt.count() << "," << res.status << "\n";

    if (run == 0 && !success)
      AdjustParams(params.ncall, params.t_iter);
    if (success)
      run++;
  } while (success == false &&
           CanAdjustNcallOrIters(params.ncall, params.t_iter) == true);

  return success;
}

template <typename F, int ndim>
bool
signle_invocation_time_and_call(F integrand,
                                double epsrel,
                                double correct_answer,
                                char const* integralName,
                                VegasParams& params,
                                quad::Volume<double, ndim>* volume,
                                int num_repeats = 100)
{

  bool success = false;
  for (int i = 0; i < num_repeats; ++i) {

    constexpr bool MCUBES_DEBUG = false;
    using MilliSeconds =
      std::chrono::duration<double, std::chrono::milliseconds::period>;
    // We make epsabs so small that epsrel is always the stopping condition.
    double constexpr epsabs = 1.0e-20;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto res = cuda_mcubes::integrate<F, ndim, MCUBES_DEBUG, Custom_generator>(
      integrand,
      epsrel,
      epsabs,
      params.ncall,
      volume,
      params.t_iter,
      params.num_adjust_iters,
      params.num_skip_iters);
    MilliSeconds dt = std::chrono::high_resolution_clock::now() - t0;
    success = (res.status == 0);
    std::cout.precision(15);

    std::cout << "estimates:" << std::scientific << res.estimate << ","
              << params.ncall << "\t" << dt.count() << std::endl;
    /*if (success)
      std::cout << integralName << "," << epsrel << "," << std::scientific
                << correct_answer << "," << std::scientific << res.estimate
                << "," << std::scientific << res.errorest << "," << res.chi_sq
                << "," << params.t_iter << "," << params.num_adjust_iters << ","
                << params.num_skip_iters << "," << res.iters << ","
                << params.ncall << "," << res.neval << "," << dt.count() << ","
                << res.status << "\n";*/
  }

  return success;
}

template <typename F, int ndim>
void
call_mcubes_kernel(int num_repeats)
{
  std::array<double, 4> required_ncall = {1.e8, 1.e9, 2.e9, 3.e9};
  double ncall = 1.0e8;
  int titer = 1;
  int itmax = 1;
  int skip = 0;
  VegasParams params(ncall, titer, itmax, skip);
  F integrand;
  quad::Volume<double, ndim> volume;
  size_t run = 0;
  double epsrel = 1.e-3;
  double true_value = 0.;
  for (auto num_samples : required_ncall) {
    params.ncall = num_samples;

    signle_invocation_time_and_call<F, ndim>(
      integrand, epsrel, true_value, "f", params, &volume, num_repeats);
    run++;
    if (run > required_ncall.size())
      break;
  }
}

#endif
