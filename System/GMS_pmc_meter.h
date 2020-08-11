
#ifndef __GMS_PMC_METER_H__
#define __GMS_PMC_METER_H__


namespace file_info {

     const unsigned int gGMS_PMC_METER_MAJOR = 1U;
     const unsigned int gGMS_PMC_METER_MINOR = 0U;
     const unsigned int gGMS_PMC_METER_MICRO = 0U;
     const unsigned int gGMS_PMC_METER_FULLVER =
         1000U*gGMS_PMC_METER_MAJOR+100U*gGMS_PMC_METER_MINOR+10U*gGMS_PMC_METER_MICRO;
     const char * const pgGMS_PMC_METER_CREATION_DATE = "09-08-2020 10:19 +00200 (SAT 09 AUG 2020 GMT+2)";
     const char * const pgGMS_PMC_METER_BUILD_DATE    = __DATE__ ":" __TIME__;
     const char * const pgGMS_PMC_METER_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
     const char * const pgGMS_PMC_METER_SYNOPSIS      = "PMC measurements (libpfc) based wrapper class";
}

#include <vector>
#include <string>
#include <cstdint>
#include "GMS_config.h"
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>   // IEEE754_2008 accurate floating-point aruthmetics.
#endif
#include "libpfc.h"

namespace gms {

            namespace system {

	              struct PMCMeter {

		            using fptr = void(*)(void);
			    
			    std::size_t        m_nruns;
			    
			    fptr               m_ptfunc;

			    std::string        m_func_name;

			    std::string        m_file_name;

			    char * __restrict  m_event1_name;

			    char * __restrict  m_event2_name;

			    char * __restrict  m_event3_name;

			    char * __restrict  m_event4_name;

			    bool         m_warmup;

			    std::vector<PFC_CNT> m_instr_issued;

			    std::vector<PFC_CNT> m_core_cycles;

			    std::vector<PFC_CNT> m_ref_cycles;

			    std::vector<PFC_CNT> m_hw_event1;

			    std::vector<PFC_CNT> m_hw_event2;

			    std::vector<PFC_CNT> m_hw_event3;

			    std::vector<PFC_CNT> m_hw_event4;

			    PMCMeter() = default;

			    PMCMeter(const std::size_t,
			             const fptr,
				     const char * __restrict,
				     const char * __restrict,
				     const char * __restrict,
				     const char * __restrict,
				     const char * __restrict,
				     const char * __restrict,
				     const bool);
				    

			    PMCMeter(const PMCMeter &) = delete;

			    PMCMeter(PMCMeter &&) noexcept(true) = delete;

			    ~PMCMeter();

			    PMCMeter & operator=(const PMCMeter &)  = delete;

			    PMCMeter & operator=(PMCMeter &&) noexcept(true) = delete;

			    int32_t run_benchmark();

#if (USE_ACCURATE_IEEE754_2008_FP) == 1

			    static bool       compute_stats(const std::vector<PFC_CNT> &,
						     double &,
			                             double &,
						     double &,
						     double &,
						     double &);
#else
			      static bool       compute_stats(const std::vector<PFC_CNT> &,
			                             double &,
			                             double &,
						     double &,
						     double &,
						     double &);
#endif

			    void       print() const;

			    static     bool delta_values_eq(std::vector<bool> &,
			                                    const PMCMeter &,
							    const PMCMeter &);

			    static     bool delta_values_ineq(std::vector<bool> &,
			                                      const PMCMeter &,
							      const PMCMeter &);

			    static     bool delta_values_gt(std::vector<bool> &,
			                                    const PMCMeter &,
							    const PMCMeter &);

			    static     bool delta_values_lt(std::vector<bool> &,
			                                    const PMCMeter &,
							    const PMCMeter &);
	     };

       } // system
  
} // gms


#endif /*__GMS_PMC_METER_H__*/
