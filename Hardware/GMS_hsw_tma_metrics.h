
#ifndef __GMS_HSW_TMA_METRICS_H__
#define __GMS_HSW_TMA_METRICS_H__ 50520201732


namespace file_info {

  const unsigned int gGMS_HSW_TMA_METRICS_MAJOR = 1;
  const unsigned int gGMS_HSW_TMA_METRICS_MINOR = 0;
  const unsigned int gGMS_HSW_TMA_METRICS_MICRO = 0;
  const unsigned int gGMS_HSW_TMA_METRICS_FULLVER =
      1000U*gGMS_HSW_TMA_METRICS_MAJOR+
      100U*gGMS_HSW_TMA_METRICS_MINOR+
      10U*gGMS_HSW_TMA_METRICS_MICRO;
  const char * const pgGMS_HSW_TMA_METRICS_CREATION_DATE = "05-05-2020 17:36 +00200 (05 MAY 2020 5:36PM GMT+2)";
  const char * const pgGMS_HSW_TMA_METRICS_BUILD_DATE    = __DATE__ ":" __TIME__;
  const char * const pgGMS_HSW_TMA_METRICS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_HSW_TMA_METRICS_SYNOPSIS      = "Haswell client performance metrics based on TMA-Metrics (4.0)";

}

#include <cstdint>




   const static uint32_t Issue_Width = 4;
   const static uint32_t Mem_L2_Store_Cost = 9;
   const static uint32_t Mem_L3_Weight = 7;
   const static uint32_t Energy_Unit = 61;
   const static uint32_t BAClear_Cost = 12;
   const static uint32_t MS_Switch_Cost = 2;
   const static uint32_t Avg_Assist_Cost = 100;
   const static uint32_t Mem_L3_Weight = 7;
   const static uint32_t Mem_STLB_Hit_Cost = 8;
   const static uint32_t Mem_XSNP_HitM_Cost = 60;
   const static uint32_t Mem_XSNP_Hit_Cost = 43;
   const static uint32_t Mem_XSNP_None_Cost = 29;


   // TMA Metrics.

namespace gms {

   static inline
   uint64_t hsw_fetched_uops(const uint64_t idq_dsb_uops,
                         const uint64_t lsd_uops,
			 const uint64_t idq_mite_uops,
			 const uint64_t idq_ms_uops) {
        return (idq_dsb_uops +
	        lsd_uops     +
		idq_mite_uops +
		idq_ms_uops);
   }

   static inline
   uint64_t hsw_recovery_cycles(const uint64_t int_misc_recovery_cycles_any,
                            const uint64_t int_misc_recovery_cycles,
			    const bool is_ht_enabled) {
         return (is_ht_enabled ? int_misc_recovery_cycles_any/2ULL :
	                         int_misc_recovery_cycles);
   }

   static inline
   uint64_t hsw_execute_cycles(const uint64_t uops_executed_core_c1,
                           const bool is_ht_enabled) {
         return (is_ht_enabled ? uops_executed_core_c1 / 2ULL :
	                         uops_executed_core_c1);
   }

   static inline
   uint64_t hsw_sq_full_cycles(const uint64_t offcore_requests_buffer_sq_full,
                           const bool is_ht_enabled) {
         return (is_ht_enabled ? offcore_requests_buffer_sq_full / 2ULL :
	                         offcore_requests_buffer_sq_full);
   }

   static inline
   uint64_t hsw_itlb_miss_cycles(const uint64_t itlb_misses_stlb_hit,
                                 const uint64_t itlb_misses_walk_duration) {
         return (14ULL*itlb_misses_stlb_hit+itlb_misses_walk_duration);
   }

   static inline
   uint64_t hsw_frontend_rs_empty_cycles(const uint64_t rs_event_empty_cycles,
                                         const float frontend_latency) {
         return (frontend_latency<0.1f ? rs_event_empty_cycles :
	                                 0ULL);
   }

   static inline
   uint64_t hsw_cycles_0_ports_utilized(const uint64_t uops_executed_core_i1_c1,
                                        const uint64_t stalls_total,
					const uint64_t rs_event_empty_cycles,
					const float frontend_latency,
					const bool is_ht_enabled) {
         return (is_ht_enabled ? uops_executed_core_c1_i1/2ULL :
	                         stalls_total-hsw_frontend_rs_empty_cycles(rs_event_empty_cycles,frontend_latency));
   }

   static inline
   uint64_t hsw_cycles_1_ports_utilized(const uint64_t uops_executed_core_c1,
                                       const uint64_t uops_executed_core_c2,
				       const bool is_ht_enabled) {
         return (is_ht_enabled ? (uops_executed_core_c1-uops_executed_core_c2)/2ULL :
	                          uops_executed_core_c1-uops_executed_core_c2);
   }

   static inline
   uint64_t hsw_cycles_2_ports_utilized(const uint64_t uops_executed_core_c2,
                                            const uint64_t uops_executed_core_c3,
					    const bool is_ht_enabled) {
          return (is_ht_enabled ? (uops_executed_core_c2-uops_executed_core_c3)/2ULL :
	                           uops_executed_core_c2-uops_executed_core_c3);
   }

   static inline
   uint64_t hsw_cycles_3_ports_utilized(const uint64_t uops_executed_core_c3,
                                        const bool is_ht_enabled) {
          return (is_ht_enabled ? uops_executed_core_c3/2ULL :
	                          uops_executed_core_c3);
   }

#include <algorithm>

   static inline
   uint64_t hsw_frontend_latency_cycles(const uint64_t cpu_clk_unhalted_thread,
                                        const uint64_t idq_uops_not_delivered_cycles_0_uops_deliv_core) {
           return (std::min(cpu_clk_unhalted_thread,
	                    idq_uops_not_delivered_cycles_0_uops_deliv_core));
   }

   static inline
   uint64_t hsw_stalls_mem_any(const uint64_t cpu_clk_unhalted_thread,
                               const uint64_t cycles_activity_stalls_lm_pending) {
           return (std::min(cpu_clk_unhalted_thread,
	                    cycles_activity_stalls_lm_pending));
   }

   static inline
   uint64_t hsw_stalls_total(const uint64_t cpu_clk_unhalted_thread,
                             const uint64_t cycles_activity_cycles_no_execute) {
           return (std::min(cpu_clk_unhalted_thread,
	                    cycles_activity_cycles_no_execute));
   }

   static inline
   uint64_t hsw_oro_drd_any_cycles(const uint64_t cpu_clk_unhalted_thread,
                                   const uint64_t offcore_requests_oustanding_cycles_with_data_rd) {
           return (std::min(cpu_clk_unhalted_thread,
	                    offcore_requests_oustanding_cycles_with_data_rd ));
   }

   static inline
   uint64_t hsw_oro_drd_bw_cycles(const uint64_t cpu_clk_unhalted_thread,
                                  const uint64_t offcore_requests_outstanding_all_data_rd_c6) {
           return (std::min(cpu_clk_unhalted_thread,
	                    offcore_requests_outstanding_all_data_rd_c6));
   }

   static inline
   uint64_t hsw_oro_demand_rfo_c1(const uint64_t cpu_clk_unhalted_thread,
                                  const uint64_t offcore_requests_outstanding_cycles_with_demand_rfo) {
           return (std::min(cpu_clk_unhalted_thread,
	                    offcore_requests_outstanding_cycles_with_demand_rfo));
   }

   static inline
   float hsw_store_l2_hit_cycles(const uint64_t l2_rqsts_rfo_hit,
                                    const uint64_t mem_uops_retired_lock_loads,
				    const uint64_t mem_uops_retired_all_stores) {
           return ((float)l2_rqsts_rfo_hit*Mem_L2_Store_Cost*
	                (1.0f-hsw_mem_lock_st_fraction(mem_uops_retired_loack_loads,
	                                               mem_uops_retired_all_stores)));
   }

   static inline
   uint64_t hsw_load_l1_miss(const uint64_t mem_load_uops_retired_l2_hit,
                             const uint64_t mem_load_uops_retired_l3_hit,
			     const uint64_t mem_load_uops_l3_hit_retired_xsnp_hit,
			     const uint64_t mem_load_uops_l3_hit_retired_xsnp_hitm,
			     const uint64_t mem_load_uops_l3_hit_retired_xsnp_miss) {
           return (mem_load_uops_retired_l2_hit   +
	           mem_load_uops_retired_l3_hit   +
		   mem_load_uops_l3_hit_retired_xsnp_hit +
		   mem_load_uops_l3_hit_retired_xsnp_hitm +
		   mem_load_uops_l3_hit_retired_xsnp_miss);
   }

   static inline
   uint64_t hsw_load_l1_miss_net(const uint64_t mem_load_uops_retired_l3_miss,
                                const uint64_t mem_load_uops_retired_l2_hit,
                                const uint64_t mem_load_uops_retired_l3_hit,
			        const uint64_t mem_load_uops_l3_hit_retired_xsnp_hit,
			        const uint64_t mem_load_uops_l3_hit_retired_xsnp_hitm,
			        const uint64_t mem_load_uops_l3_hit_retired_xsnp_miss ) {
            return (hsw_load_l1_miss(mem_load_uops_retired_l2_hit,
	                             mem_load_uops_retired_l3_hit,
				     mem_load_uops_l3_hit_retired_xsnp_hit,
				     mem_load_uops_l3_hit_retired_xsnp_hitm,
				     mem_load_uops_l3_hit_retired_xsnp_miss) +
				     mem_load_uops_retired_l3_miss);
   }

   static inline
   float hsw_load_l3_hit(const uint64_t mem_load_uops_retired_l3_hit,
                         const uint64_t mem_load_uops_retired_hit_lfb,
			 const uint64_t mem_load_uops_retired_l2_hit,
			 const uint64_t mem_load_uops_l3_hit_retired_xsnp_hit,
			 const uint64_t mem_load_uops_l3_hit_retired_xsnp_hitm,
			 const uint64_t mem_load_uops_l3_hit_retired_xsnp_miss) {
             return ((float)mem_load_uops_retired_l3_hit*
	                    (1+mem_load_uops_retired_hit_lfb) /
			    hsw_load_l1_miss(mem_load_uops_retired_l2_hit,
	                                     mem_load_uops_retired_l3_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hitm,
				             mem_load_uops_l3_hit_retired_xsnp_miss) );
    }

    static inline
    float hsw_load_xsnp_hit(const uint64_t mem_load_uops_l3_hit_retired_xsnp_hit,
                            const uint64_t mem_load_uops_retired_hit_lfb,
			    const uint64_t mem_load_uops_retired_l2_hit,
			    const uint64_t mem_load_uops_retired_l3_hit,
			    const uint64_t mem_load_uops_l3_hit_retired_xsnp_hitm,
			    const uint64_t mem_load_uops_l3_hit_retired_xsnp_miss) {
              return ((float)mem_load_uops_l3_hit_retired_xsnp_hit *
	                     (1+mem_load_uops_retired_hit_lfb) /
			      hsw_load_l1_miss(mem_load_uops_retired_l2_hit,
	                                     mem_load_uops_retired_l3_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hitm,
				             mem_load_uops_l3_hit_retired_xsnp_miss) );
    }

    static inline
    float hsw_load_xsnp_hitm(
                             const uint64_t mem_load_uops_retired_hit_lfb,
			     const uint64_t mem_load_uops_retired_l2_hit,
	                     const uint64_t mem_load_uops_retired_l3_hit,
			     const uint64_t mem_load_uops_l3_hit_retired_xsnp_hit,
			     const uint64_t mem_load_uops_l3_hit_retired_xsnp_hitm,
			     const uint64_t mem_load_uops_l3_hit_hit_retired_xsnp_miss) {
                return ((float)mem_load_uops_l3_hit_retired_xsnp_hitm *
		               (1+mem_load_uops_retired_hit_lfb) /
			        hsw_load_l1_miss(mem_load_uops_retired_l2_hit,
	                                     mem_load_uops_retired_l3_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hitm,
				             mem_load_uops_l3_hit_retired_xsnp_miss) );
    }

    static inline
    float hsw_load_xsnp_miss( const uint64_t mem_load_uops_retired_hit_lfb,
			     const uint64_t mem_load_uops_retired_l2_hit,
	                     const uint64_t mem_load_uops_retired_l3_hit,
			     const uint64_t mem_load_uops_l3_hit_retired_xsnp_hit,
			     const uint64_t mem_load_uops_l3_hit_retired_xsnp_hitm,
			     const uint64_t mem_load_uops_l3_hit_hit_retired_xsnp_miss) {
                return ((float)mem_load_uops_l3_hit_retired_xsnp_miss*
		               (1+mem_load_uops_retired_hit_lfb) /
			        hsw_load_l1_miss(mem_load_uops_retired_l2_hit,
	                                     mem_load_uops_retired_l3_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hit,
				             mem_load_uops_l3_hit_retired_xsnp_hitm,
				             mem_load_uops_l3_hit_retired_xsnp_miss) );
    }

    static inline
    uint64_t hsw_few_uops_executed_threshold(const uint64_t uops_executed_core_c2,
                                             const uint64_t uops_executed_core_c3,
					     const float ipc) {
                return (ipc > 1.5f ? uops_executed_core_c3 :
		                     uops_executed_core_c2);
     }

     static inline
     float hsw_backend_bound_cycles(const uint64_t stalls_total,
                                    const uint64_t uops_executed_core_c1,
				    const uint64_t few_uops_executed_threshold,
				    const uint64_t frontend_rs_empty_cycles,
				    const uint64_t resource_stall_sb,
				    const bool is_ht_enabled) {
                return (is_ht_enabled ? (float)(stalls_total+uops_executed_core_c1-few_uops_executed_threshold)/2 -
		                               frontend_rs_empty_cycles+resource_stalls_sb :
					(float)(stalls_total+uops_executed_core_c1-few_uops_executed_threshold-
					        frontend_rs_empty_cycles+resource_stalls_sb));
     }

     static inline
     float hsw_memory_bound_fraction(const uint64_t stalls_mem_any,
                                     const uint64_t resource_stalls_sb,
				     const float backend_bound_cycles) {
                return ((float)(stalls_mem_any+resource_stalls_sb) /
		                backend_bound_cycles);
     }

     static inline
     float hsw_mem_l3_hit_fraction( const uint64_t mem_load_uops_retired_l3_hit,
                                    const uint64_t mem_load_uops_retired_l3_miss) {
                return ((float)mem_load_uops_retired_l3_hit / (mem_load_uops_retired_l3_hit+Mem_L3_Weight+
		                                               mem_load_uops_retired_l3_miss));
     }

     static inline
     float hsw_mem_lock_st_fraction( const uint64_t mem_uops_retired_lock_loads,
                                     const uint64_t mem_uops_retired_all_stores) {
                return ((float)mem_uops_retired_lock_loads/mem_uops_retired_all_stores);
     }

     static inline
     float hsw_mispred_clears_fraction( const uint64_t br_mispred_all_branches,
                                        const uint64_t machine_clears_count) {
                return ((float)br_mispred_all_branches/(br_mispred_all_branches+
		                                        machine_clear_count));
     }

     static inline
     float hsw_retire_fraction( const uint64_t uops_retired_retire_slots,
                                const uint64_t uops_issued_any) {
                return ((float)uops_retired_retire_slots/uops_issued_any);
     }

    // clks is CLOCK_UNHALTED.THREAD
    static inline
    float hsw_ipc( const uint64_t inst_retired_any,
                   const uint64_t clks) {
                return ((float)inst_retired_any/clks);
     }

     static inline
     float hsw_upi( const uint64_t uops_retired_retire_slots,
                    const uint64_t uops_retired_any) {
                return ((float)uops_retired_retire_slots/uops_retired_any);
     }

     // Instructions per taken branch
     static inline
     float hsw_iptb( const uint64_t instr_retired_any,
                     const uint64_t br_instr_retired_near_taken) {
                return ((float)instr_retired_any/br_instr_retired_near_taken);
     }

     static inline
     float hsw_cpi( const uint64_t instr_retired_any,
                    const uint64_t clks) {
                return (1.0f/hsw_ipc(istr_retired_any,clks));
     }

     static inline
     uint64_t hsw_issue_slots( const uint64_t core_clks) {
                return (Pipeline_Width*core_clks);
     }
     // Instructions per load
     static inline
     float hsw_ipload( const uint64_t instr_retired_any,
                       const uint64_t mem_uops_retired_all_loads) {
                return ((float)instr_retired_any/mem_uops_retired_all_loads);
     }
     // Instructions per store
     static inline
     float hsw_ipstore( const uint64_t instr_retired_any,
                        const uint64_t mem_uops_retired_all_stores) {
                return ((float)instr_retired_any/mem_uops_retired_all_stores);
     }

     static inline
     float hsw_ipbranch( const uint64_t instr_retired_any,
                         const uint64_t br_instr_retired_all_branches) {
                return ((float)instr_retired_any/br_instr_retired_all_branches);
     }
     
     static inline
     float hsw_ipcall( const uint64_t instr_retired_any,
                       const uint64_t br_instr_retired_near_call) {
                return ((float)instr_retired_any/br_instr_retired_near_call);
     }
     // Branch instuctions per taken branch
     static inline
     float hsw_biptb( const uint64_t br_inst_retired_all_branches,
                      const uint64_t br_inst_retired_near_taken) {
                return ((float)br_inst_retired_all_branches /
	                       br_inst_retired_near_taken);
     }

     static inline
     float hsw_dsb_coverage( const uint64_t idq_dsb_uops,
                             const uint64_t fetched_uops) {
                return ((float)idq_dsb_uops/fetched_uops);
     }

     static inline
     float hsw_ipbaclear( const uint64_t inst_retired_any,
                          const uint64_t baclears_any) {
                return ((float)inst_retired_any/baclears_any);
     }

     static inline
     float hsw_ipc_core( const uint64_t instr_retired_any,
                         const uint64_t core_clks) {
                return ((float)instr_retired_any/core_clks);
     }

     static inline
     float hsw_ilp( const uint64_t uops_executed_core,
                    const uint64_t execute_cycles,
		    const bool is_ht_enabled) {
                return (is_ht_enabled ? (float) uops_executed_core/2/execute_cycles :
		                        (float) uops_executed_core/execute_cycles);
     }

     static inline
     float hsw_ip_mispredict( const uint64_t inst_retired_any,
                              const uint64_t br_misp_retired_all_branches) {
                 return ((float)inst_retired_any/br_misp_retired_all_branches);
     }

     static inline
     uint64_t hsw_core_clks( const uint64_t cpu_clk_unhalted_thread,
                          const uint64_t cpu_clk_unhalted_one_thread_active,
			  const uint64_t cpu_clk_unhalted_ref_xclk,
			  const uint64_t cpu_clk_unhalted_thread_any,
			  const bool ebs_mode) {
                 const uint64_t t0 = cpu_clk_unhalted_thread/2ULL *
		                     ((1ULL+cpu_clk_unhalted_one_thread_active)/cpu_clk_unhalted_ref_xclk);
		 const uint64_t t1 = cpu_clk_unhalted_thread_any/2ULL;
		 return (ebs_mode ? t0 : t1);
     }

     static inline
     float hsw_load_miss_real_latency(const uint64_t l1d_pend_miss_pending,
                                      const uint64_t mem_load_uops_retired_l1_miss,
				      const uint64_t mem_load_uops_retired_hit_lfb) {
                 return ((float)l1d_pend_miss_pending /
		                (mem_load_uops_retired_l1_miss +
				 mem_load_uops_retired_hit_lfb));
     }

     static inline
     float hsw_mem_level_parallelism( const uint64_t l1d_pend_miss_pending,
                                      const uint64_t l1d_pend_miss_pending_cycles) {
                 return ((float)l1d_pend_miss_pending/
		                l1d_pend_miss_pending_cycles);
     }

     static inline
     float hsw_page_walk_util( const uint64_t itlb_misses_walk_duration,
                               const uint64_t dtlb_load_misses_walk_duration,
			       const uint64_t dtlb_store_misses_walk_duration,
			       const uint64_t core_clks) {
                  return ((float)(itlb_misses_walk_duration+dtlb_load_misses_walk_duration+
		                  dtlb_store_misses_walk_duration)/core_clks);
     }
     // Average cache fill bandwith Gigabytes/sec
     static inline
     float hsw_l1d_cache_fill_bw( const uint64_t l1d_replacement,
                                  const uint64_t time_iterval) {
                  return ((float)64ULL*l1d_replacement/1000000000ULL/time_interval);
     }

     static inline
     float hsw_l2_cache_fill_bw( const uint64_t l2_lines_in_all,
                                 const uint64_t time_interval) {
                  return ((float)64ULL*l2_lines_all/1000000000ULL/time_interval);
     }

     static inline
     float hsw_l3_cache_fill_bw( const uint64_t longest_lat_cache_miss,
                                 const uint64_t time_interval) {
                  return ((float)64ULL*longest_lat_cache_miss/1000000000ULL/time_interval);
     }

     static inline
     float hsw_l1mpki( const uint64_t mem_load_uops_retired_l1_miss,
                       const uint64_t inst_retired_any) {
                  return ((float)1000ULL*mem_load_uops_retired_l1_miss/
		                         inst_retired_any);
     }

     static inline
     float hsw_l2mpki( const uint64_t mem_load_uops_retired_l2_miss,
                       const uint64_t  inst_retired_any) {
                  return ((float)1000ULL*mem_load_uops_retired_l2_miss/
		                         inst_retired_any);
     }

     static inline
     float hsw_l2hpki(  const uint64_t mem_load_uops_retired_l2_miss,
                       const uint64_t inst_retired_any) {
                  return ((float)1000ULL*mem_load_uops_retired_l2_miss/
		                         inst_retired_any);
     }

     static inline
     float hsw_l3mpki( const uint64_t mem_load_uops_retired_l3_miss,
                       const uint64_t inst_retired_any) {
                 return ((float)1000ULL*mem_load_uops_retired_l3_miss/
		                         inst_retired_any);
     }

     static inline
     float hsw_ipfarbr( const uint64_t inst_retired_any,
                        const uint64_t br_inst_retired_far_branch) {
                 return ((float)inst_retired_any/
		                (br_inst_retired_far_branch));
     }

     static inline
     float hsw_assist( const uint64_t other_assists_any_wb_assists,
                       const uint64_t slots) {
                return ((float)Avg_Assist_Cost*other_assists_any_wb_assists / slots);
     }

     // Metrics based on 'pmu-tools'

     static inline
     float hsw_kernel_utilization(const uint64_t cpu_clk_unhalted_thread_sup,
                                  const uint64_t cpu_clk_unhalted_ref_tsc) {
                 return ((float)cpu_clk_unhalted_thread_sup/
		                cpu_clk_unhalted_ref_tsc);
     }

     static inline
     float hsw_dram_bw_use(const uint64_t unc_arb_trk_requests_all,
                           const uint64_t unc_arb_coh_trk_requests_all,
			   const uint64_t time_interval) {
                 return ((float)64ULL*(unc_arb_trk_requests_all+
		                       unc_arb_coh_trk_requests_all)/1000000000ULL/time_interval);
     }

     static inline
     float hsw_mem_requests_latency( const uint64_t unc_arb_trk_occupancy_all,
                                 const uint64_t unc_arb_trk_requests_all) {
                  return ((float)unc_arb_trk_occupancy_all/
		                 unc_arb_trk_requests_all);
     }

     static inline
     float hsw_mem_parallel_requests( const uint64_t unc_arb_trk_occupancy_all,
                                      const uint64_t unc_arb_trk_occupancy_cycles_with_any_request) {
                   return ((float)unc_arb_trk_occupancy_all/
		           unc_arb_trk_occupancy_cycles_with_any_requests);
     }

     static inline
     float hsw_ht_utilization( const uint64_t cpu_clk_thread_unhalted_one_thread_active,
                               const uint64_t cpu_clk_thread_unhalted_ref_xclk_any,
			       const bool is_ht_enabled) {
                   return (is_ht_enabled ? (float)cpu_clk_thread_unhalted_one_thread_active/
		                                  cpu_clk_thread_unhalted_ref_xclk_any :
						  0.0f);
     }

     static inline
     float hsw_frontend_bound( const uint64_t idq_uops_not_deliverd_core,
                               const uint64_t slots) {
                   return ((float)idq_uops_not_deliverd_core/slots);
     }
    // Domain: pipeline slots
     static inline
     float hsw_frontend_latency( const uint64_t frontend_latency_cycles,
                                 const uint64_t cycles) {
                   return ((float)(4ULL*frontend_latency_cycles)/slots);
     }

     static inline
     float hsw_icache_misses( const uint64_t icache_ifdata_stall,
                              const uint64_t clks) {
                   return ((float)icache_ifdata_stall/clks);
     }

     static inline
     float hsw_dsb_switches(const uint64_t dsb2mite_switches_penalty_cycles,
                            const uint64_t clks) {
                   return ((float)dsb2mite_switches_penalty_cycles/clks);
     }

     static inline
     float hsw_lcp(const uint64_t ild_stall_lcp,
                   const uint64_t clks) {
                   return ((float)ild_stall_lcp/clks);
     }

     static inline
     float hsw_ms_switches( const uint64_t idq_ms_switches,
                            const uint64_t clks) {
                   return ((float)(Ms_Switches_Cost*idq_ms_switches)/clks);
     }

     static inline
     float hsw_branch_resteers( const uint64_t br_misp_retired_all_branches,
                                const uint64_t machine_clears_count,
				const uint64_t baclears_any,
				const uint64_t clks) {
                    return ((float)BAClear_Cost*(br_misp_retired_all_branches+
		                                 machine_clears_count+baclears_any)/clks)
     }

     static inline
     float hsw_mite( const uint64_t idq_all_mite_cycles_any_uops,
                     const uint64_t idq_all_mite_cycles_4_uops,
		     const uint64_t core_clks) {
                    return ((float) (idq_all_mite_cycles_any_uops+
		                     idq_all_mite_cycles_4_uops)/core_clks);
     }

     static inline
     float hsw_dsb( const uint64_t idq_all_dsb_cycles_any_uops,
                    const uint64_t idq_all_dsb_cycles_4_uops,
		    const uint64_t core_clks) {
                   return ((float)(idq_all_dsb_cycles_any_uops+
		                   idq_all_dsb_cycles_4_uops)/core_clks);
     }

     static inline
     float hsw_l1_bound( const uint64_t stalls_mem_any,
                         const uint64_t cycles_activity_stalls_l1d_pending,
			 const uint64_t clks) {
                  return ((float)(stalls_mem_any-
		                  cycles_activity_stalls_l1d_pending)/clks);
     }

     static inline
     float hsw_dtlb_load(const uint64_t dtlb_load_misses_stlb_hit,
                         const uint64_t dtlb_load_misses_walk_duration,
			 const uint64_t clks) {
                  return ((float)(Mem_STLB_Hit_Cost*
		                  dtlb_load_misses_stlb_hit+
				  dtlb_load_misses_walk_duration)/clks);
     }

     static inline
     float hsw_store_fwd_blk( const uint64_t ld_blocks_store_forward,
                              const uint64_t clks) {
                   return ((float)(13ULL*ld_blocks_store_forward)/clks);
     }

     static inline
     float hsw_split_loads( const float load_miss_real_latency,
                            const uint64_t ld_blocks_no_sr,
			    const uint64_t clks) {
                   return ((float)(load_miss_real_latency*ld_blocks_no_sr)/clks);
     }
     

} // gms    
   










#endif /*__GMS_HSW_TMA_METRICS_H__*/
