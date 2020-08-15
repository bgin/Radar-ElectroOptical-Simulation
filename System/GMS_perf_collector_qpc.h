
#ifndef __GMS_PERF_COLLECTOR_QPC_H__
#define __GMS_PERF_COLLECTOR_QPC_H__



namespace file_info {
  const unsigned int gGMS_PERF_COLLECTOR_QPC_MAJOR = 1U;
  const unsigned int gGMS_PERF_COLLECTOR_QPC_MINOR = 0U;
  const unsigned int gGMS_PERF_COLLECTOR_QPC_MICRO = 0U;
  const unsigned int gGMS_PERF_COLLECTOR_QPC_FULLVER =
    1000U*gGMS_PERF_COLLECTOR_QPC_MAJOR + 100U*gGMS_PERF_COLLECTOR_QPC_MINOR + 10U*gGMS_PERF_COLLECTOR_QPC_MICRO;
  const char * const pgGMS_PERF_COLLECTOR_QPC_CREATE_DATE = "05-10-2019 12:12 +00200 (SAT 05 OCT 2019 GMT+2)";
  const char * const pgGMS_PERF_COLLECTOR_QPC_BUILD_DATE  = "00-00-0000 00:00";
  const char * const pgGMS_PERF_COLLECTOR_QPC_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_PERF_COLLECTOR_QPC_SYNOPSIS    = "Code section performance collector class based on QueryPerformanceCounter.";
}

#include <vector>
#include <string>
#include <cstdint>
#include <Windows.h>
#include "../GMS_config.h"
#if (USE_MKL) == 1 && defined (GMS_COMPILED_BY_ICC)
#include <mkl_dfti.h>
#endif
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>
#endif

void rdft(int, int, double *, int *, double *);

namespace gms {
	namespace system {

		struct PerfCollectorQPC {
				
				bool		 m_Iscleared;

				std::size_t  m_nsamples; // Counter of valid measurements.

				std::string  m_func_name;

				std::string  m_file_name;

				std::string  m_date;

				std::string  m_time;

				// this member describe successful computation of delta value (success == true, otherwise false)
				std::vector<bool>          m_nvalid_values;
				// profiling samples  start  values
				std::vector<LARGE_INTEGER> m_start_values;
				// profiling samples, stop values
				std::vector<LARGE_INTEGER> m_stop_values;
				// Time samples delta
				std::vector<double>        m_delta_values;
				//  Performance Counter values themselves 
				std::vector<LARGE_INTEGER> m_pc_values;

				constexpr static size_t lo_bound = 3Ui64;

				constexpr static size_t lowest = 1Ui64;

				constexpr static size_t zero = 0Ui64;

				PerfCollectorQPC();

				PerfCollectorQPC(_In_ const PerfCollectorQPC &);

				PerfCollectorQPC(_In_ PerfCollectorQPC &&) noexcept(true);

				~PerfCollectorQPC() noexcept(true) = default;

				PerfCollectorQPC &  operator=(_In_ const PerfCollectorQPC &);

				PerfCollectorQPC &  operator=(_In_ PerfCollectorQPC &&) noexcept(true);

				bool start() noexcept(false);

				bool stop()  noexcept(false);

				void clear_all();

				bool compute_delta() noexcept(false);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
				
				// This version relies on Intel libbfp754 and that means
				// probably slower execution speed (not confirmed yet).
				bool compute_stats(_Out_ double &,
								   _Out_ double &,
								   _Out_ double &,
								   _Out_ double &,
								   _Out_ double &);

#else
				bool compute_stats(_Out_ double &,
								   _Out_ double &,
								   _Out_ double &,
								   _Out_ double &,
								   _Out_ double &);

#endif
				
				/*
					If using MKL for FFT pass null pointers to
					correlation_set (those pointers are used to
					call Ouru FFT rdft function)
					Usage of MKL DFTI is set as a default.
				*/
				static bool correlation_set(PerfCollectorQPC &,
							    PerfCollectorQPC &,
							    double * __restrict,
							    int32_t * __restrict,
							    double * __restrict,
							    const int32_t );

				void print() const;

				bool check_data_size()const noexcept(true);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
				// This version relies on Intel libbfp754 and that means
				// probably slower execution speed (not confirmed yet).
				bool cvrt_to_double(double * __restrict,
						    const std::size_t) noexcept(true);
#else
				bool cvrt_to_double(double * __restrict,
					            const std::size_t) noexcept(true);
#endif		 



		};
	}
}




#endif /*__GMS_PERF_COLLECTOR_QPC_H__*/
