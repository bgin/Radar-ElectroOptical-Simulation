

#ifndef __GMS_PERF_COLLECTOR_PMC_H__
#define __GMS_PERF_COLLECTOR_PMC_H__


namespace file_info {

    const unsigned int gGMS_PERF_COLLECTOR_PMC_MAJOR = 1;
    const unsigned int gGMS_PERF_COLLECTOR_PMC_MINOR = 0;
    const unsigned int gGMS_PERF_COLLECTOR_PMC_MICRO = 0;
    const unsigned int gGMS_PERF_COLLECTOR_PMC_FULLVER =
      1000U*gGMS_PERF_COLLECTOR_PMC_MAJOR+100U*gGMS_PERF_COLLECTOR_PMC_MINOR+10U*gGMS_PERF_COLLECTOR_PMC_MICRO;
    const char * const pgGMS_PERF_COLLECTOR_PMC_CREATE_DATE = "14-08-2020 16:32 +00200 (FRI 14 AUG 2020 GMT+2)";
    const char * const pgGMS_PERF_COLLECTOR_PMC_BUILD_DATE  = __DATE__ ":" __TIME__;
    const char * const pgGMS_PERF_COLLECTOR_PMC_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_PERF_COLLECTOR_PMC_SYNOPSIS    = "Code section performance collector based on PMC (libpfc)";
}


#include <vector>
#include <string>
#include <cstdint>
#include "GMS_config.h"
#include "libpfc.h"
#if (USE_MKL) == 1 && defined (GMS_COMPILED_BY_ICC)
#include <mkl_dfti.h>
#endif
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>
#endif

void rdft(int, int, double *, int *, double *);

namespace gms {

      namespace system {

                 struct PerfCollectorPMC {

                    std::size_t m_nsamples;

		    std::string m_func_name;

		    std::string m_file_name;

		    std::string m_date;

		    std::string m_time;

		    char * __restrict m_event1_name;

		    char * __restrict m_event2_name;

		    char * __restrict m_event3_name;

		    char * __restrict m_event4_name;

		    bool       m_iscleared;

		    std::vector<PFC_CNT> m_instr_issued;

		    std::vector<PFC_CNT> m_core_cycles;

		    std::vector<PFC_CNT> m_ref_cycles;

		    std::vector<PFC_CNT> m_hw_event1;

		    std::vector<PFC_CNT> m_hw_event2;

		    std::vector<PFC_CNT> m_hw_event3;

		    std::vector<PFC_CNT> m_hw_event4;

		    constexpr static size_t lo_bound = 3Ui64;

		    constexpr static size_t lowest = 1Ui64;

		    constexpr static size_t zero = 0Ui64;

		    PerfCollectorPMC() = default;

		    PerfCollectorPMC(const char * __restrict,
		                     const char * __restrict,
				     const char * __restrict,
				     const char * __restrict,
                                     const char * __restrict,
				     const char * __restrict,
				     const char * __restrict,
				     const char * __restrict,
				     const int32_t);

		    PerfCollectorPMC(const PerfCollectorPMC &) = delete;

		    PerfCollectorPMC(PerfCollectorPMC &&) = delete;

		    PerfCollectorPMC & operator=(const PerfCollectorPMC &) = delete;

		    PerfCollectorPMC & operator=(PerfCollectorPMC &&) = delete;

		    ~PerfCollectorPMC();

		    void start();

		    void stop();

		    void clear_all();
#if (USE_ACCURATE_IEEE754_2008_FP) == 1

                    static bool compute_stats(const std::vector<PFC_CNT> &,
		                              double &,
					      double &,
					      double &,
					      double &,
					      double &);
#else
                    
		    static bool compute_stats(const std::vector<PFC_CNT> &,
		                              double &,
					      double &,
					      double &,
					      double &,
					      double &);
#endif

                    /*
					If using MKL for FFT pass null pointers to 
					correlation_set (those pointers are used to
					                 call Ouru FFT rdft function)
					Usage of MKL DFTI is set as a default.
		    */

		    
	 };
                   
   } // system
 
}// gms

#endif /*__GMS_PERF_COLLECTOR_PMC_H__*/
