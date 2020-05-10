
#include "GMS_hsw_tma_metrics_api.h"

// Used only with integral types
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)



  
   uint64_t hsw_fetched_uops(const uint64_t idq_dsb_uops,
                             const uint64_t lsd_uops,
			     const uint64_t idq_mite_uops,
			     const uint64_t idq_ms_uops) {
        return (idq_dsb_uops +
	        lsd_uops     +
		idq_mite_uops +
		idq_ms_uops);
   }

  
   uint64_t hsw_recovery_cycles(const uint64_t int_misc_recovery_cycles_any,
                            const uint64_t int_misc_recovery_cycles,
			    const bool is_ht_enabled) {
         return (is_ht_enabled ? int_misc_recovery_cycles_any/2ULL :
	                         int_misc_recovery_cycles);
   }

   
   uint64_t hsw_execute_cycles(const uint64_t uops_executed_core_c1,
                           const bool is_ht_enabled) {
         return (is_ht_enabled ? uops_executed_core_c1 / 2ULL :
	                         uops_executed_core_c1);
   }

  
   uint64_t hsw_sq_full_cycles(const uint64_t offcore_requests_buffer_sq_full,
                           const bool is_ht_enabled) {
         return (is_ht_enabled ? offcore_requests_buffer_sq_full / 2ULL :
	                         offcore_requests_buffer_sq_full);
   }

   
   uint64_t hsw_itlb_miss_cycles(const uint64_t itlb_misses_stlb_hit,
                                 const uint64_t itlb_misses_walk_duration) {
         return (14ULL*itlb_misses_stlb_hit+itlb_misses_walk_duration);
   }

  
   uint64_t hsw_frontend_rs_empty_cycles(const uint64_t rs_event_empty_cycles,
                                         const float frontend_latency) {
         return (frontend_latency<0.1f ? rs_event_empty_cycles :
	                                 0ULL);
   }

  
   uint64_t hsw_cycles_0_ports_utilized(const uint64_t uops_executed_core_i1_c1,
                                        const uint64_t stalls_total,
					const uint64_t rs_event_empty_cycles,
					const float frontend_latency,
					const bool is_ht_enabled) {
         return (is_ht_enabled ? uops_executed_core_c1_i1/2ULL :
	                         stalls_total-hsw_frontend_rs_empty_cycles(rs_event_empty_cycles,frontend_latency));
   }

  
   uint64_t hsw_cycles_1_ports_utilized(const uint64_t uops_executed_core_c1,
                                       const uint64_t uops_executed_core_c2,
				       const bool is_ht_enabled) {
         return (is_ht_enabled ? (uops_executed_core_c1-uops_executed_core_c2)/2ULL :
	                          uops_executed_core_c1-uops_executed_core_c2);
   }

  
   uint64_t hsw_cycles_2_ports_utilized(const uint64_t uops_executed_core_c2,
                                            const uint64_t uops_executed_core_c3,
					    const bool is_ht_enabled) {
          return (is_ht_enabled ? (uops_executed_core_c2-uops_executed_core_c3)/2ULL :
	                           uops_executed_core_c2-uops_executed_core_c3);
   }

  
   uint64_t hsw_cycles_3_ports_utilized(const uint64_t uops_executed_core_c3,
                                        const bool is_ht_enabled) {
          return (is_ht_enabled ? uops_executed_core_c3/2ULL :
	                          uops_executed_core_c3);
   }



  
   uint64_t hsw_frontend_latency_cycles(const uint64_t cpu_clk_unhalted_thread,
                                        const uint64_t idq_uops_not_delivered_cycles_0_uops_deliv_core) {
           return (MIN(cpu_clk_unhalted_thread,
	                    idq_uops_not_delivered_cycles_0_uops_deliv_core));
   }

  
   uint64_t hsw_stalls_mem_any(const uint64_t cpu_clk_unhalted_thread,
                               const uint64_t cycles_activity_stalls_lm_pending) {
           return (MIN(cpu_clk_unhalted_thread,
	                    cycles_activity_stalls_lm_pending));
   }

  
   uint64_t hsw_stalls_total(const uint64_t cpu_clk_unhalted_thread,
                             const uint64_t cycles_activity_cycles_no_execute) {
           return (MIN(cpu_clk_unhalted_thread,
	                    cycles_activity_cycles_no_execute));
   }

  
   uint64_t hsw_oro_drd_any_cycles(const uint64_t cpu_clk_unhalted_thread,
                                   const uint64_t offcore_requests_oustanding_cycles_with_data_rd) {
           return (MIN(cpu_clk_unhalted_thread,
	                    offcore_requests_oustanding_cycles_with_data_rd ));
   }

  
   uint64_t hsw_oro_drd_bw_cycles(const uint64_t cpu_clk_unhalted_thread,
                                  const uint64_t offcore_requests_outstanding_all_data_rd_c6) {
           return (MIN(cpu_clk_unhalted_thread,
	                    offcore_requests_outstanding_all_data_rd_c6));
   }

  
   uint64_t hsw_oro_demand_rfo_c1(const uint64_t cpu_clk_unhalted_thread,
                                  const uint64_t offcore_requests_outstanding_cycles_with_demand_rfo) {
           return (MIN(cpu_clk_unhalted_thread,
	                    offcore_requests_outstanding_cycles_with_demand_rfo));
   }

  
   float hsw_store_l2_hit_cycles(const uint64_t l2_rqsts_rfo_hit,
                                    const uint64_t mem_uops_retired_lock_loads,
				    const uint64_t mem_uops_retired_all_stores) {
           return ((float)l2_rqsts_rfo_hit*Mem_L2_Store_Cost*
	                (1.0f-hsw_mem_lock_st_fraction(mem_uops_retired_loack_loads,
	                                               mem_uops_retired_all_stores)));
   }

  
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

   
    uint64_t hsw_few_uops_executed_threshold(const uint64_t uops_executed_core_c2,
                                             const uint64_t uops_executed_core_c3,
					     const float ipc) {
                return (ipc > 1.5f ? uops_executed_core_c3 :
		                     uops_executed_core_c2);
     }

   
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

   
     float hsw_memory_bound_fraction(const uint64_t stalls_mem_any,
                                     const uint64_t resource_stalls_sb,
				     const float backend_bound_cycles) {
                return ((float)(stalls_mem_any+resource_stalls_sb) /
		                backend_bound_cycles);
     }

    
     float hsw_mem_l3_hit_fraction( const uint64_t mem_load_uops_retired_l3_hit,
                                    const uint64_t mem_load_uops_retired_l3_miss) {
                return ((float)mem_load_uops_retired_l3_hit / (mem_load_uops_retired_l3_hit+Mem_L3_Weight+
		                                               mem_load_uops_retired_l3_miss));
     }

   
     float hsw_mem_lock_st_fraction( const uint64_t mem_uops_retired_lock_loads,
                                     const uint64_t mem_uops_retired_all_stores) {
                return ((float)mem_uops_retired_lock_loads/mem_uops_retired_all_stores);
     }

    
     float hsw_mispred_clears_fraction( const uint64_t br_mispred_all_branches,
                                        const uint64_t machine_clears_count) {
                return ((float)br_mispred_all_branches/(br_mispred_all_branches+
		                                        machine_clear_count));
     }

   
     float hsw_retire_fraction( const uint64_t uops_retired_retire_slots,
                                const uint64_t uops_issued_any) {
                return ((float)uops_retired_retire_slots/uops_issued_any);
     }

    // clks is CLOCK_UNHALTED.THREAD
  
    float hsw_ipc( const uint64_t inst_retired_any,
                   const uint64_t clks) {
                return ((float)inst_retired_any/clks);
     }

    
     float hsw_upi( const uint64_t uops_retired_retire_slots,
                    const uint64_t uops_retired_any) {
                return ((float)uops_retired_retire_slots/uops_retired_any);
     }

     // Instructions per taken branch
    
     float hsw_iptb( const uint64_t instr_retired_any,
                     const uint64_t br_instr_retired_near_taken) {
                return ((float)instr_retired_any/br_instr_retired_near_taken);
     }

    
     float hsw_cpi( const uint64_t instr_retired_any,
                    const uint64_t clks) {
                return (1.0f/hsw_ipc(istr_retired_any,clks));
     }

   
     uint64_t hsw_issue_slots( const uint64_t core_clks) {
                return (Pipeline_Width*core_clks);
     }
     // Instructions per load
    
     float hsw_ipload( const uint64_t instr_retired_any,
                       const uint64_t mem_uops_retired_all_loads) {
                return ((float)instr_retired_any/mem_uops_retired_all_loads);
     }
     // Instructions per store
    
     float hsw_ipstore( const uint64_t instr_retired_any,
                        const uint64_t mem_uops_retired_all_stores) {
                return ((float)instr_retired_any/mem_uops_retired_all_stores);
     }

   
     float hsw_ipbranch( const uint64_t instr_retired_any,
                         const uint64_t br_instr_retired_all_branches) {
                return ((float)instr_retired_any/br_instr_retired_all_branches);
     }
     
   
     float hsw_ipcall( const uint64_t instr_retired_any,
                       const uint64_t br_instr_retired_near_call) {
                return ((float)instr_retired_any/br_instr_retired_near_call);
     }
     // Branch instuctions per taken branch
    
     float hsw_biptb( const uint64_t br_inst_retired_all_branches,
                      const uint64_t br_inst_retired_near_taken) {
                return ((float)br_inst_retired_all_branches /
	                       br_inst_retired_near_taken);
     }

   
     float hsw_dsb_coverage( const uint64_t idq_dsb_uops,
                             const uint64_t fetched_uops) {
                return ((float)idq_dsb_uops/fetched_uops);
     }

    
     float hsw_ipbaclear( const uint64_t inst_retired_any,
                          const uint64_t baclears_any) {
                return ((float)inst_retired_any/baclears_any);
     }

     
     float hsw_ipc_core( const uint64_t instr_retired_any,
                         const uint64_t core_clks) {
                return ((float)instr_retired_any/core_clks);
     }

    
     float hsw_ilp( const uint64_t uops_executed_core,
                    const uint64_t execute_cycles,
		    const bool is_ht_enabled) {
                return (is_ht_enabled ? (float) uops_executed_core/2/execute_cycles :
		                        (float) uops_executed_core/execute_cycles);
     }

    
     float hsw_ip_mispredict( const uint64_t inst_retired_any,
                              const uint64_t br_misp_retired_all_branches) {
                 return ((float)inst_retired_any/br_misp_retired_all_branches);
     }

    
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

    
     float hsw_load_miss_real_latency(const uint64_t l1d_pend_miss_pending,
                                      const uint64_t mem_load_uops_retired_l1_miss,
				      const uint64_t mem_load_uops_retired_hit_lfb) {
                 return ((float)l1d_pend_miss_pending /
		                (mem_load_uops_retired_l1_miss +
				 mem_load_uops_retired_hit_lfb));
     }

    
     float hsw_mem_level_parallelism( const uint64_t l1d_pend_miss_pending,
                                      const uint64_t l1d_pend_miss_pending_cycles) {
                 return ((float)l1d_pend_miss_pending/
		                l1d_pend_miss_pending_cycles);
     }

    
     float hsw_page_walk_util( const uint64_t itlb_misses_walk_duration,
                               const uint64_t dtlb_load_misses_walk_duration,
			       const uint64_t dtlb_store_misses_walk_duration,
			       const uint64_t core_clks) {
                  return ((float)(itlb_misses_walk_duration+dtlb_load_misses_walk_duration+
		                  dtlb_store_misses_walk_duration)/core_clks);
     }
     // Average cache fill bandwith Gigabytes/sec
    
     float hsw_l1d_cache_fill_bw( const uint64_t l1d_replacement,
                                  const uint64_t time_iterval) {
                  return ((float)64ULL*l1d_replacement/1000000000ULL/time_interval);
     }

   
     float hsw_l2_cache_fill_bw( const uint64_t l2_lines_in_all,
                                 const uint64_t time_interval) {
                  return ((float)64ULL*l2_lines_all/1000000000ULL/time_interval);
     }

    
     float hsw_l3_cache_fill_bw( const uint64_t longest_lat_cache_miss,
                                 const uint64_t time_interval) {
                  return ((float)64ULL*longest_lat_cache_miss/1000000000ULL/time_interval);
     }

    
     float hsw_l1mpki( const uint64_t mem_load_uops_retired_l1_miss,
                       const uint64_t inst_retired_any) {
                  return ((float)1000ULL*mem_load_uops_retired_l1_miss/
		                         inst_retired_any);
     }

    
     float hsw_l2mpki( const uint64_t mem_load_uops_retired_l2_miss,
                       const uint64_t  inst_retired_any) {
                  return ((float)1000ULL*mem_load_uops_retired_l2_miss/
		                         inst_retired_any);
     }

     
     float hsw_l2hpki(  const uint64_t mem_load_uops_retired_l2_miss,
                       const uint64_t inst_retired_any) {
                  return ((float)1000ULL*mem_load_uops_retired_l2_miss/
		                         inst_retired_any);
     }

   
     float hsw_l3mpki( const uint64_t mem_load_uops_retired_l3_miss,
                       const uint64_t inst_retired_any) {
                 return ((float)1000ULL*mem_load_uops_retired_l3_miss/
		                         inst_retired_any);
     }

     
     float hsw_ipfarbr( const uint64_t inst_retired_any,
                        const uint64_t br_inst_retired_far_branch) {
                 return ((float)inst_retired_any/
		                (br_inst_retired_far_branch));
     }

    
     float hsw_assist( const uint64_t other_assists_any_wb_assists,
                       const uint64_t slots) {
                return ((float)Avg_Assist_Cost*other_assists_any_wb_assists / slots);
     }


