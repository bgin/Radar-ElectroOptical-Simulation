
#if defined _WIN64
    #ifndef __GMS_GTC64_TIMER_H__
    #define __GMS_GTC64_TIMER_H__

// File version
namespace file_info {
      const unsigned int gGMS_GTC64_TIMER_MAJOR = 1U;

      const unsigned int gGMS_GTC64_TIMER_MINOR = 0U;

      const unsigned int gGMS_GTC64_TIMER_MICRO = 0U;

      const unsigned int gGMS_GTC64_TIMER_FULLVER = 
		1000U*gGMS_GTC64_TIMER_MAJOR + 100U*gGMS_GTC64_TIMER_MINOR + 10U*gGMS_GTC64_TIMER_MICRO;

      const char * const pgGMS_GTC64_TIMER_CREATE_DATE = "26-08-2018 11:07 +00200 (SUN 26 AUG 2018 GMT+2) ";

      const char * const pgGMS_GTC64_TIMER_BUILD_DATE = "00-00-0000 00:00 ";

      const char * const pgGMS_GTC64_TIMER_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

      const char * const pgGMS_GTC64_TIMER_SYNOPSIS = "GetTickCount64 timer benchmarking based class.";
}




#include <vector>
#include <string>
#include <Windows.h>
#include "../LAM_config.h"

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>     // Accurate IEEE754 conforming arithmetics
#endif


namespace gms {
	namespace system {

			struct GTC64Timer {
				
				using fptr = void(*)();
				

				std::size_t  m_nruns;

				std::size_t  m_datum_size;

				fptr		 m_ptfunc;// Timed function called by the wrapper

				fptr		 m_ptwrap; // Empty wrapper body (expect prolog and epilog overhead)

				std::string  m_func_name;

				std::string  m_file_name;

				

				std::vector<ULONGLONG> m_timing_values;

				std::vector<ULONGLONG> m_overhead_values;

				std::vector<ULONGLONG> m_delta_values; // Delta between timing values and wrapper overhead calls.

				static constexpr ULONGLONG lo_bound = 3Ui64;

				static constexpr ULONGLONG zero = 0Ui64;

				GTC64Timer() = default;

				GTC64Timer(_In_ const std::size_t,
						   _In_ const std::size_t,
						   _In_ const fptr tf,
						   _In_ const fptr wf,
						   _In_ const char *,
						   _In_ const char *);

				GTC64Timer(_In_ const GTC64Timer &);

				GTC64Timer(_In_ GTC64Timer &&) noexcept(true);

				~GTC64Timer() noexcept(true);

				GTC64Timer & operator=(_In_ const GTC64Timer &);

				GTC64Timer & operator=(_In_ GTC64Timer &&) noexcept(true);

				inline std::size_t get_datum_size() const { return (m_nruns * m_datum_size); }

				void run_benchmark() noexcept(false);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
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

				void print() const;

				static bool   delta_values_eq(_Out_ std::vector<bool> &,
								       _In_  const GTC64Timer &,
									   _In_  const GTC64Timer &);

				static bool   delta_values_ineq(_Out_ std::vector<bool> &,
										 _In_ const GTC64Timer &,
										 _In_ const GTC64Timer &);

				static bool   delta_values_gt(_Out_ std::vector<bool> &,
									   _In_ const GTC64Timer &,
									   _In_ const GTC64Timer &);

				static bool   delta_values_lt(_Out_ std::vector<bool> &,
									   _In_ const GTC64Timer &,
									   _In_ const GTC64Timer &);


		};
	}
}


      #endif /*__GMS_GTC64_TIMER_H__*/
#endif
