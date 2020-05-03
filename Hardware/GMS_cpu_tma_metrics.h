
#ifndef __GMS_CPU_TMA_METRICS_H__
#define __GMS_CPU_TMA_METRICS_H__

namespace file_info {

  const unsigned int gGMS_CPU_TMA_METRICS_MAJOR = 1;
  const unsigned int gGMS_CPU_TMA_METRICS_MINOR = 0;
  const unsigned int gGMS_CPU_TMA_METRICS_MICRO = 0;
  const unsigned int gGMS_CPU_TMA_METRICS_FULLVER =
    1000U*gGMS_CPU_TMA_METRICS_MAJOR+100U*gGMS_CPU_TMA_METRICS_MINOR+10U*gGMS_CPU_TMA_METRICS_MICRO;
  const char * const pgGMS_CPU_TMA_METRICS_CREATE_DATE = "03-05-2020 10:34  +00200 (03 MAY 2020 10:34AM GMT+2)";
  const char * const pgGMS_CPU_TMA_METRICS_BUILD_DATE  = __DATE__ ":" __TIME__;
  const char * const pgGMS_CPU_TMA_METRICS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_CPU_TMA_METRICS_SYNOPSYS    = "CPU performance computation metrics (Intel only) based on TMA-Metrics";
 
}

#include <csdtint.h>

namespace gms {

  /*
       Front-End metrics
   */

   // Front-End bound Skylake/KabyLake
   // ICache Misses
   // ( ICACHE_16B.IFDATA_STALL + 2 * ICACHE_16B.IFDATA_STALL:c1:e1 ) / CLKS
       static inline
       float skl_kbl_icache_misses( const uint64_t ifdata_stall,
                                    const uint64_t ifdata_stall_c1e1,
				    const uint64_t clks) {
           return ((float)(ifdata_stall+2*ifdata_stall_c1e1)/clks);
       }

    //ITLB misses
    // ICACHE_64B.IFTAG_STALL / CLKS
    static inline
    float skl_kbl_itlb_misses( const uint64_t iftag_stall,
                               const uint64_t clks) {
           return ((float)iftag_stall/clks);
      }

    // Branch resteers
    // INT_MISC.CLEAR_RESTEER_CYCLES + #BAClear_Cost * BACLEARS.ANY ) / CLKS
    static inline
    float skl_kbl_branch_resteers( const uint64_t clear_resteer_cycles,
                                   const uint64_t ba_clear_cost,
				   const uint64_t baclears_any,
				   const uint64_t clks) {
            return ((float)(clear_resteer_cycles+ba_clear_cost*baclears_any)/clks);
      }

     // DTLB Load
     // ( #Mem_STLB_Hit_Cost * DTLB_LOAD_MISSES.STLB_HIT + DTLB_LOAD_MISSES.WALK_ACTIVE ) / CLKS
     static inline
     float skl_kbl_dtlb_load( const uint64_t mem_stlb_hit_cost,
                              const uint64_t stlb_hit,
			      const uint64_t walk_active,
			      const uint64_t clks) {
            return ((float)(mem_stlb_hit_cost*stlb_hit+walk_active)/clks);
      }

     //L2 Bound
     // ( 1 if FB_Full < 1.5 else #LOAD_L2_HIT / ( #LOAD_L2_HIT + L1D_PEND_MISS.FB_FULL:c1 ) ) * #L2_Bound_Ratio
     static inline
     float skl_kbl_l2_bound(const float fb_full,
                            const uint64_t load_l2_hit,
			    const uint64_t fb_full_c1,
			    const float l2_bound_ratio) {
            const uint64_t res =  fb_full < 1.5f ? 1:load_l2_hit;
	    return ((float)res/(load_l2_hit+fb_full_c1)*l2_bound_ratio);
     }

     //L3 Bound
     // ( CYCLE_ACTIVITY.STALLS_L2_MISS - CYCLE_ACTIVITY.STALLS_L3_MISS ) / CLKS
     static inline
     float skl_kbl_l3_bound(const uint64_t stalls_l2_miss,
                            const uint64_t stalls_l3_miss,
			    const uint64_t clks) {
            return ((float)(stalls_l2_miss-stalls_l3_miss)/clks);
     }

     // Divider
     // ARITH.DIVIDER_ACTIVE / CLKS
     static inline
     float skl_kbl_divider(const uint64_t divider_active,
                           const uint64_t clks) {
            return ((float)divider_active/clks);
     }

     // Ports Utilization
     // ( #Backend_Bound_Cycles - CYCLE_ACTIVITY.STALLS_MEM_ANY - EXE_ACTIVITY.BOUND_ON_STORES ) / CLKS if ( ARITH.DIVIDER_ACTIVE < EXE_ACTIVITY.EXE_BOUND_0_PORTS ) else ( #Backend_Bound_Cycles - CYCLE_     // ACTIVITY.STALLS_MEM_ANY - EXE_ACTIVITY.BOUND_ON_STORES - EXE_ACTIVITY.EXE_BOUND_0_PORTS ) / CLKS
     static inline
     float skl_kbl_ports_utilization(const uint64_t backend_bound_cycles,
                                     const uint64_t stalls_mem_any,
				     const uint64_t bound_on_stores,
				     const uint64_t clks,
				     const uint64_t divider_active,
				     const uint64_t exe_bound_0_ports) {
            float res - 0.0f;
	    res = divider_active<exe_bound_0_ports ?
	                  (float)(backend_bound_cycles-stalls_mem_any-bound_on_stores)/clks :
			  (float)(backend_bound_cycles-stalls_mem_any-bound_on_stores-exe_bound_0_ports)/clks;
	    return (res);
      }

      //Serializng operations
      // PARTIAL_RAT_STALLS.SCOREBOARD / CLKS
      static inline
      float skl_kbl_serializing_ops(const uint64_t scoreboard,
                                    const uint64_t clks) {
            return ((float)scoreboard/clks);
      }
      
#include <math.h>

      // IFetch line utilization
      // min( 1 , UOPS_ISSUED.ANY / ( UPI * 64 * ( ICACHE_64B.IFTAG_HIT + ICACHE_64B.IFTAG_MISS ) / 4.1 ) )
      static inline
      float skl_kbl_ifetch_line_utilization(const uint64_t uops_issued_any,
                                            const uint64_t upi,
					    const uint64_t iftag_hit,
					    const uint64_t iftag_miss) {
              return (std::min(1.0f,(float)uops_issued_any(upi*64*(iftag_hit+iftag_miss)/4.1f)));
      }

      // iPL
      // INST_RETIRED.ANY / MEM_INST_RETIRED.ALL_LOADS_PS
      static inline
      float skl_kbl_ipl(const uint64_t inst_retired_any,
                        const uint64_t all_load_ps) {
              return ((float)inst_retired_any/all_load_ps);
      }

      // iPS
      // INST_RETIRED.ANY / MEM_INST_RETIRED.ALL_STORES_PS
      static inline
      float skl_kbl_ips(const uint64_t inst_retired_any,
                        const uint64_t all_store_ps) {
              return ((float)inst_retired_any/all_store_ps);
      }

      // Load miss real latency
      // L1D_PEND_MISS.PENDING / ( MEM_LOAD_RETIRED.L1_MISS_PS + MEM_LOAD_RETIRED.FB_HIT_PS )
      static inline
      float skl_kbl_load_miss_real_latency(const uint64_t miss_pending,
                                           const uint64_t l1_miss_ps,
					   const uint64_t fb_hit_ps) {
              return ((float)miss_pending/(l1_miss_ps+fb_hit_ps));
      }

      // Page walk utilization
      /// ( ITLB_MISSES.WALK_PENDING + DTLB_LOAD_MISSES.WALK_PENDING + DTLB_STORE_MISSES.WALK_PENDING + EPT.WALK_PENDING ) / ( 2 * CORE_CLKS )
      static inline
      float skl_kbl_page_walk_utilization(const uint64_t itlb_misses_walk_pending,
                                          const uint64_t dtlb_load_walk_pending,
					  const uint64_t dtlb_store_miss_walk_pending,
					  const uint64_t ept_walk_pending,
					  const uint64_t core_clks) {
              const uint64_t t0 = itlb_misses_walk_pending+
	                          dtlb_load_walk_pending+
				  dtlb_store_miss_walk_pending+
				  ept_walk_pending;
	      return ((float)t0/(2ULL*core_clks));
      }

      // L3 cache access BW
      // 64 * OFFCORE_REQUESTS.ALL_REQUESTS / #OneBillion / Time
      static inline
      float skl_kbl_l3_cache_access_bw(const uint64_t offcore_request_all_requests,
                                       const uint64_t time_period) {
               return ((float)64ULL*offcore_requests_all_request/1000000000ULL/time_period);
      }

      // L1 MPKI
      // 1000 * MEM_LOAD_RETIRED.L1_MISS_PS / INST_RETIRED.ANY
      static inline
      float skl_kbl_l1_mpki(const uint64_t mem_load_retired_l1_miss_ps,
                            const uint64_t inst_retired_any) {
               return ((float)1000ULL*mem_load_retired_l1_miss_ps/inst_retired_any);
      }

      // L2 MPKI
      // 1000 * MEM_LOAD_RETIRED.L2_MISS_PS / INST_RETIRED.ANY
      static inline
      float skl_kbl_l2_mpki(const uint64_t mem_load_retired_l2_miss_ps,
                            const uint64_t inst_retired_any) {
               return ((float)1000ULL*mem_load_retired_l2_miss_ps/inst_retired_any);
      }

      // L3 MPKI
      // 1000 * MEM_LOAD_RETIRED.L3_MISS_PS / INST_RETIRED.ANY
      static inline
      float skl_kbl_l3_mpki(const uint64_t mem_load_retired_l3_miss_ps,
                            const uint64_t inst_retired_any) {
               return ((float)1000ULL*mem_load_retired_l3_miss_ps/inst_retired_any);
      }

      // DRAM latency
      // OneBillion * ( UNC_ARB_TRK_OCCUPANCY.DATA_READ / UNC_ARB_TRK_REQUESTS.DATA_READ ) / ( Socket_CLKS / Time )
      static inline
      float skl_kbl_dram_latency(const uint64_t unc_arb_trk_occupancy_data_read,
                                 const uint64_t unc_arb_trk_requests_data_read,
				 const uint64_t socket_clks,
				 const uint64_t time_period) {
                return ((float)1000000000ULL*(unc_arb_trk_occupancy_data_read/unc_arb_trk_requests_data_read)/
		                              (socket_clks/time_period));
      }

      // Data parallel reads
      // UNC_ARB_TRK_OCCUPANCY.DATA_READ / UNC_ARB_TRK_OCCUPANCY.DATA_READ:c1
      static inline
      float skl_kbl_data_parallel_reads(const uint64_t unc_arb_trk_occupancy_data_read,
                                        const uint64_t unc_arb_trk_occupancy_data_read_c1) {
                return ((float)unc_arb_trk_occupancy_data_read/
		               unc_arb_trk_occupancy_data_read_c1);
      }

      // IpFarBranch
      // INST_RETIRED.ANY / ( BR_INST_RETIRED.FAR_BRANCH_PS / 2 )
      static inline
      float skl_kbl_ip_far_branch(const uint64_t inst_retired_any,
                                  const uint64_t br_inst_retired_far_branch_ps) {
                return ((float)instr_retired_any/
		               (br_inst_retired_far_branch_ps/2ULL));
      }

      // Fetched uops
      // IDQ.DSB_UOPS + IDQ.MITE_UOPS + IDQ.MS_UOPS
      static inline
      uint64_t skl_kbl_fetched_uops(const uint64_t idq_dsb_uops,
                                    const uint64_t idq_mite_uops,
				    const uint64_t idq_ms_uops) {
                return (idq_dsb_uops+idq_mite_uops+idq_ms_uops);      
      }

      // Execute cycles
      // ( UOPS_EXECUTED.CORE_CYCLES_GE_1 / 2 ) if #SMT_on else UOPS_EXECUTED.CORE_CYCLES_GE_1
      static inline
      uint64_t skl_kbl_executed_cycles(const uint64_t uops_executed_core_ge_1,
                                       const uint64_t uops_executed_core_cycles_ge_1,
				       const bool is_smt_on) {
                return ( is_smt_on ? uops_executed_core_ge_1/ 2 :
		                     uops_executed_core_cycles_ge_1);
      }

      // ITLB miss cycles
      // ( 14 * ITLB_MISSES.STLB_HIT + ITLB_MISSES.WALK_ACTIVE )
      static inline
      uint64_t skl_kbl_itlb_miss_cycles(const uint64_t itlb_misses_stlb_hit,
                                        const uint64_t itlb_misses_walk_active) {
                return (14ULL * itlb_misses_stlb_hit +
		                itlb_misses_walk_active);
      }

      // Load L1 miss
      // MEM_LOAD_RETIRED.L2_HIT_PS + MEM_LOAD_RETIRED.L3_HIT_PS + MEM_LOAD_L3_HIT_RETIRED.XSNP_HIT_PS +
      // MEM_LOAD_L3_HIT_RETIRED.XSNP_HITM_PS + MEM_LOAD_L3_HIT_RETIRED.XSNP_MISS_PS
      static inline
      uint64_t skl_kbl_load_l1_miss(const uint64_t mem_load_retired_l2_hit_ps,
                                    const uint64_t mem_load_retired_l3_hit_ps,
				    const uint64_t mem_load_l3_hit_retired_xsnp_hit_ps,
				    const uint64_t mem_load_l3_hit_retired_xsnp_hitm_ps,
				    const uint64_t mem_load_l3_hit_retired_xsnp_miss_ps) {
                 return (mem_load_retired_l2_hit_p            +
		         mem_load_retired_l3_hit_ps           +
			 mem_load_l3_hit_retired_xsnp_hit_ps  +
			 mem_load_l3_hit_retired_xsnp_hitm_ps +
			 mem_load_l3_hit_retired_xsnp_miss_ps);
      }
}





#endif /*__GMS_CPU_TMA_METRICS_H__*/
