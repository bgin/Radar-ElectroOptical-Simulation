
#ifndef __GMS_PERF_COLLECTOR_GTC64_H__
#define __GMS_PERF_COLLECTOR_GTC64_H__




namespace file_info {
  const unsigned int gGMS_PERF_COLLECTOR_GTC64_MAJOR = 1U;
  const unsigned int gGMS_PERF_COLLECTOR_GTC64_MINOR = 0U;
  const unsigned int gGMS_PERF_COLLECTOR_GTC64_MICRO = 0U;
  const unsigned int gGMS_PERF_COLLECTOR_GTC64_FULLVER =
    1000U*gGMS_PERF_COLLECTOR_GTC64_MAJOR + 100U*gGMS_PERF_COLLECTOR_GTC64_MINOR + 10U*gGMS_PERF_COLLECTOR_GTC64_MICRO;
  const char * const pgGMS_PERF_COLLECTOR_GTC64_CREATION_DATE = "05-10-2019 14:25 +00200 (SAT 05 OCT 2019 GMT+2)";
  const char * const pgGMS_PERF_COLLECTOR_GTC64_BUILD_DATE    = "00-00-0000 00:00 +00200";
  const char * const pgGMS_PERF_COLLECTOR_GTC64_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_PERF_COLLECTOR_GTC64_SYNOPSIS      = "Code section performance collector class based on Win API GetTickCount64 function";
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

namespace lam {
	namespace system  {



		struct PerfCollectorGTC64 {

				bool  m_Iscleared;

				std::size_t  m_nsamples;

				std::string  m_funcname;

				std::string  m_filename;

				std::string  m_date;

				std::string  m_time;

				std::vector<bool>  m_nsuccess; // this member describe successful computation of delta value (success == true, otherwise false)

				std::vector<ULONGLONG>  m_start_values;

				std::vector<ULONGLONG>  m_stop_values;

				std::vector<ULONGLONG>  m_delta_values;

				constexpr static size_t  lo_bound = 3Ui64;

				constexpr static size_t  lowest = 1Ui64;

				constexpr static size_t  zero = 0Ui64;

				PerfCollectorGTC64();

				PerfCollectorGTC64(_In_ const PerfCollectorGTC64 &);

				PerfCollectorGTC64(_In_ PerfCollectorGTC64 &&) noexcept(true);

				~PerfCollectorGTC64() noexcept(true) = default;

				PerfCollectorGTC64 & operator=(_In_ const PerfCollectorGTC64 &);

				PerfCollectorGTC64 & operator=(_In_ PerfCollectorGTC64 &&) noexcept(true);

				void start() noexcept(false);

				void stop()  noexcept(false);

				void clear_all();

				bool compute_delta() noexcept(false);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
                // This version relies on Intel libbfp754 and that means
				// probably slower execution speed (not confirmed yet).
				bool  compute_stats(_Out_ double &,
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
				static bool correlation_set(_In_ PerfCollectorGTC64 &,
										    _In_ PerfCollectorGTC64 &,
											_Out_ double * __restrict,
											_In_ int32_t * __restrict,
											_In_ double * __restrict,
											_In_ const int32_t);

				void print() const;

		        bool check_data_size() noexcept(true);

		      void cvrt_to_double(_Out_ double * __restrict,
								   _In_ const std::size_t) noexcept(true);

		};
	}
}




#endif /*__GMS_PERF_COLLECTOR_GTC64_H__*/
