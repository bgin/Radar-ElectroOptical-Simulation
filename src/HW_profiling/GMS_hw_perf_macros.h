

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_HW_PERF_MACROS_H__
#define __GMS_HW_PERF_MACROS_H__


#include <cstdint>
#include "GMS_fast_pmc_access.h"


/* Basic PMC analysis setup, i.e. cycle-analysis*/

#if !defined(HW_PMC_COLLECTION_PROLOGE_BODY)
      #define HW_PMC_COLLECTION_PROLOGE_BODY             \
        do {                                                        \
               uint64_t prog_counters_start[4] = {};                \
               uint64_t prog_counters_end[4]   = {};                \
               uint64_t tsc_start,tsc_stop;                         \
               uint64_t act_cyc_start,act_cyc_stop;                 \
               uint64_t ref_cyc_start,ref_cyc_stop;                 \
               uint64_t volatile dummy1,dummy2,dummy3,dummy4;       \
               int32_t core_counter_width;                          \
               double utilization,nom_ghz,avg_ghz;                  \
               dummy1 = rdtsc();                                    \
               dummy2 = rdtsc();                                    \
               dummy3 = rdpmc(0);                                   \
               core_counter_width = get_core_counter_width();       \
               for(int32_t i = 0; i != 4; ++i) {                    \
                   prog_counters_start[i] = rdpmc(i);               \
               }                                                    \
               act_cyc_start = rdpmc_actual_cycles();               \
               ref_cyc_start = rdpmc_reference_cycles();            \
               tsc_start = rdtsc();                                
	} while(0)
#endif

#if !defined(HW_PMC_COLLECTION_EPILOGE_BODY)
    #define  HW_PMC_COLLECTION_EPILOGE_BODY                                                                                                                         \
       do {                                                                                                                                                                   \
                for(int32_t i = 0; i != 4; ++i) {                                                                                                                             \
                     prog_counters_stop[i] = rdpmc(i);                                                                                                                        \
                }                                                                                                                                                             \                                                                                                                                                                                                       \
                act_cyc_stop = rdpmc_actual_cycles();                                                                                                                         \
                ref_cyc_stop = rdpmc_reference_cycles();                                                                                                                      \
                tsc_stop = rdtsc();                                                                                                                                           \
	            dummy4 = get_core_counter_width();                                                                                                                            \
	            nom_ghz = get_TSC_frequency()/1.0e9;                                                                                                                          \
	            utilization = (double)(ref_cyc_end-ref_cyc_start)/(double)(tsc_end-tsc_start-rdtscp_latency);                                                                 \
	            avg_ghz = (double)(act_cyc_end-act_cyc_start)/(double)(tsc_end-tsc_start-rdtscp_latency)*nom_ghz;                                                             \
	            syslog(LOG_INFO,"%-10s:\n", __PRETTY_FUNCTION__);                                                                                                             \
	            syslog(LOG_INFO, "*************** Hardware Counters -- Dump Begin **************");                                                                           \
	            syslog(LOG_INFO,"Core utilization                      : %f\n",utilization );                                                                                 \
	            syslog(LOG_INFO,"Core average frequency                : %f\n",avg_ghz);                                                                                      \
	            syslog(LOG_INFO,"Reference cycles                      : %20lld\n",ref_cyc_end-ref_cyc_start);                                                                \
	            syslog(LOG_INFO,"Actual cycles                         : %20lld\n",act_cyc_stop-act_cyc_start);                                                               \
	            syslog(LOG_INFO,"L2 cycles  at entry                   : %20lld\n",cycles_lvl2);                                                                              \
	            syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event1                    , corrected_pmc_delta(prog_counters_stop[0],prog_counters_start[0],core_counter_width);      \
	            syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event2                    , corrected_pmc_delta(prog_counters_stop[1],prog_counters_start[1],core_counter_width);      \
	            syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event3                    , corrected_pmc_delta(prog_counters_stop[2],prog_counters_start[2],core_counter_width);      \
	            syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event4                    , corrected_pmc_delta(prog_counters_stop[3],prog_counters_start[3],core_counter_width);      \    
	            syslog(LOG_INFO, "*************** Hardware Counters -- Dump End   **************");
       } while(0)
#endif











#endif /*__GMS_HW_PERF_MACROS_H__*/