#ifndef GPUINTEGRATION_COMMON_INTEGRATION_RESULT_H
#define GPUINTEGRATION_COMMON_INTEGRATION_RESULT_H

#include <ostream>

namespace numint {

  // integration_result is the return type for all integration
  // routines. Not all members are filled by all algorithms.

  struct integration_result {
    double estimate = 0.;
    double errorest = 0.;
    size_t neval = 0;
    size_t nregions = 0;
    size_t nFinishedRegions = 0;
    int status = -1;
    int lastPhase = -1;
    double chi_sq = 0.;
    size_t iters = 0;
  };

  std::ostream& operator<<(std::ostream& os, integration_result const& res);
}

inline std::ostream&
numint::operator<<(std::ostream& os, numint::integration_result const& res)
{
  os << res.estimate << "," << res.errorest << "," << res.nregions << ","
     << res.chi_sq << "," << res.status;
  return os;
}

#endif
