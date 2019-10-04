
#ifndef __GMS_RDTSCP_TIMER_H__
#define __GMS_RDTSCP_TIMER_H__

// File major version
const unsigned int gGMS_RDTSCP_TIMER_MAJOR = 1U;

// File minor version
const unsigned int gGMS_RDTSCP_TIMER_MINOR = 0U;

// File micro version
const unsigned int gGMS_RDTSCP_TIMER_MICRO = 0U;

// File full version
const unsigned int gGMS_RDTSCP_TIMER_FULLVER = 
		1000U*gGMS_RDTSCP_TIMER_MAJOR + 100U*gGMS_RDTSCP_TIMER_MINOR + 10U*gGMS_RDTSCP_TIMER_MICRO;

// File creation date
const  char * const gpGMS_RDTSCP_TIMER_CREATE_DATE = "4-10-2019 17:47 +00200 (FRI 04 OCT  2019 GMT+2)";

// File build date (should be set after successful build of this translation unit.
const char * const gpGMS_RDTSCP_TIMER_BUILD_DATE = "00-00-0000 00:00";

// File author
const char * const gpGMS_RDTSCP_TIMER_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

// File short description
const char * const gpGMS_RDTSCP_TIMER_SYNOPSIS = "RDTSCP based performance measurement class.";

#include <vector>
#include <string>
#include <cstdint>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#endif
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>   // IEEE754_2008 accurate floating-point aruthmetics.
#endif

namespace gms {
	namespace system {

		struct RDTSCPTimer {

				using fptr = void(*)(void);

				std::size_t  m_nruns;

				std::size_t  m_datum_size;

				fptr         m_ptfunc; // Timed function called by the wrapper

				fptr         m_ptwrap; // Empty wrapper body (expect prolog and epilog overhead)

				std::string  m_func_name;

				std::string  m_file_name;

				std::vector<uint64_t> m_timing_values;

				std::vector<uint64_t> m_overhead_values;

				std::vector<uint64_t> m_delta_values; // Delta between timing values and wrapper overhead calls.

				std::vector<std::pair<uint32_t,uint32_t>> m_tscaux_values; // Executing core values (collected only for timing calls)

				static constexpr uint64_t lo_bound = 3Ui64;

				static constexpr uint64_t zero = 0Ui64;

				RDTSCPTimer() = default;

				RDTSCPTimer(const std::size_t,
					    const std::size_t,
					    const fptr,
					    const fptr,
					    const char *,
					    const char *);

				RDTSCPTimer(const RDTSCPTimer &);

				RDTSCPTimer(RDTSCPTimer &&) noexcept(true);

				inline ~RDTSCPTimer() noexcept(true);

				RDTSCPTimer & operator=(const RDTSCPTimer &);

				RDTSCPTimer & operator=(RDTSCPTimer &&);

				inline std::size_t get_datum_size() const { return (m_nruns*m_datum_size); }

				void run_benchmark() noexcept(false);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1

				bool  compute_stats(double &,
						    double &,
						    double &,
						    double &,
						    double & );
#else

				bool  compute_stats(double &,
						    double &,
						    double &,
						    double &,
						    double &);

#endif
				void print() const noexcept(true);

				static bool  delta_values_eq(std::vector<bool> &,
							     const RDTSCPTimer &,
							     const RDTSCPTimer &);

				static bool  delta_values_ineq(std::vector<bool> &,
							       const RDTSCPTimer &,
							       const RDTSCPTimer &);

				static bool  delta_values_gt(std::vector<bool> &,
							     const RDTSCPTimer &,
							     const RDTSCPTimer &);

				static bool  delta_values_lt(std::vector<bool> &,
							     const RDTSCPTimer &,
							     const RDTSCPTimer &);

		};
	}
}


#endif /*__GMS_RDTSCP_TIMER_H__*/
