#ifndef VERBOSE_UTILS_CUH
#define VERBOSE_UTILS_CUH

#include "cudaMemoryUtil.h"
#include "vegas_utils.cuh"

// this isn't needed anymore
std::ofstream
GetOutFileVar(std::string filename)
{
  std::ofstream myfile;
  myfile.open(filename.c_str());
  return myfile;
}

template <bool DEBUG_MCUBES>
class IterDataLogger {
  std::ofstream myfile_bin_bounds;
  std::ofstream myfile_randoms;
  std::ofstream myfile_funcevals;
  std::ofstream interval_myfile;
  std::ofstream iterations_myfile;

public:
  double* randoms = nullptr;
  double* funcevals = nullptr;

  IterDataLogger(uint32_t totalNumThreads,
                 int chunkSize,
                 int extra,
                 int npg,
                 int ndim)
  {
    if constexpr (DEBUG_MCUBES) {
      randoms = quad::cuda_malloc_managed<double>(
        (totalNumThreads * chunkSize + extra) * npg * ndim);
      funcevals = quad::cuda_malloc_managed<double>(
        (totalNumThreads * chunkSize + extra) * npg);

      myfile_bin_bounds.open("pmcubes_bin_bounds.csv");
      myfile_bin_bounds << "it, cube, chunk, sample, dim, ran00\n";
      myfile_bin_bounds.precision(15);

      myfile_randoms.open("pmcubes_random_nums.csv");
      myfile_randoms << "it, cube, chunk, sample, dim, ran00\n";
      myfile_randoms.precision(15);

      myfile_funcevals.open("pmcubes_funcevals.csv");
      myfile_funcevals << "it, cube, chunk, sample, funceval\n";

      interval_myfile.open("pmcubes_intevals.csv");
      interval_myfile.precision(15);

      iterations_myfile.open("pmcubes_iters.csv");
      iterations_myfile
        << "iter, estimate, errorest, chi_sq, iter_estimate, iter_errorest\n";
      iterations_myfile.precision(15);
    }
  }

  ~IterDataLogger()
  {
    if constexpr (DEBUG_MCUBES) {
      myfile_bin_bounds.close();
      myfile_randoms.close();
      myfile_funcevals.close();
      interval_myfile.close();
      iterations_myfile.close();
    }
  }

  void
  PrintIterResults(int iteration,
                   double estimate,
                   double errorest,
                   double chi_sq,
                   double iter_estimate,
                   double iter_errorest)
  {
    iterations_myfile << iteration << "," << estimate << "," << errorest << ","
                      << chi_sq << "," << iter_estimate << "," << iter_errorest
                      << "\n";
  }

  void
  PrintBins(int iter, double* xi, double* d, int ndim)
  {
    int ndmx1 = 501;   // Internal_Vegas_Params::get_NDMX_p1();
    int ndmx = 500;    // Internal_Vegas_Params::get_NDMX();
    int mxdim_p1 = 21; // Internal_Vegas_Params::get_MXDIM_p1();

    if (iter == 1) {
      myfile_bin_bounds
        << "iter, dim, bin, bin_length, left, right, contribution\n";
    }

    if (iter <= 2) {
      for (int dim = 1; dim <= ndim; dim++)
        for (int bin = 1; bin <= ndmx; bin++) {

          double bin_length = xi[dim * ndmx1 + bin] - xi[dim * ndmx1 + bin - 1];
          double left = xi[dim * ndmx1 + bin - 1];
          if (bin == 1)
            left = 0.;
          double right = xi[dim * ndmx1 + bin];
          double contribution = d[bin * mxdim_p1 + dim];
          myfile_bin_bounds << iter << "," << dim << "," << bin << ","
                            << bin_length << "," << left << "," << right << ","
                            << contribution << "\n";
        }
    }
  }

  void
  PrintRandomNums(int it, int ncubes, int npg, int ndim)
  {

    size_t nums_per_cube = npg * ndim;
    size_t nums_per_sample = ndim;

    if (it > 2)
      return;
    else {

      for (int cube = 0; cube < ncubes; cube++)
        for (int sample = 1; sample <= npg; sample++)
          for (int dim = 1; dim <= ndim; dim++) {

            size_t index =
              cube * nums_per_cube + nums_per_sample * (sample - 1) + dim - 1;

            myfile_randoms << it << "," << cube << "," << cube
                           << "," // same as chunk for single threaded
                           << sample << "," << dim << "," << randoms[index]
                           << "\n";
          }
    }
  }

  void
  PrintFuncEvals(int it, int ncubes, int npg, int ndim)
  {

    size_t nums_per_cube = npg * ndim;
    size_t nums_per_sample = ndim;

    if (it > 2)
      return;
    else {

      std::cout << "expecting total random numbers:" << ncubes * npg * ndim
                << "\n";
      for (int cube = 0; cube < ncubes; cube++)
        for (int sample = 1; sample <= npg; sample++) {

          size_t nums_evals_per_chunk = npg;
          size_t index = cube * nums_evals_per_chunk + (sample - 1);
          myfile_funcevals << it << "," << cube << "," << cube
                           << "," // same as chunk for single threaded
                           << sample << "," << funcevals[index] << "\n";
        }
    }
  }

  void
  PrintIntervals(int ndim,
                 int ng,
                 uint32_t totalNumThreads,
                 int chunkSize,
                 int it)
  {
    /*constexpr int mxdim_p1 = Internal_Vegas_Params::get_MXDIM_p1();

    if(it == 1)
        interval_myfile<<"m, kg[1], kg[2], kg[3], it\n";

    for(uint32_t m = 0; m < totalNumThreads; m++){
        uint32_t kg[mxdim_p1];
        get_indx(m , &kg[1], ndim, ng);

        interval_myfile<<m<<",";
            for(int ii = 1; ii<= ndim; ii++)
                interval_myfile<<kg[ii]<<",";
            interval_myfile<<it<<"\n";
    }*/
  }
};

#endif
