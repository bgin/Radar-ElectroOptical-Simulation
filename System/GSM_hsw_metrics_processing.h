
#ifndef __GMS_HSW_METRICS_PROCESSING_H__
#define __GMS_HSW_METRICS_PROCESSING_H__


namespace file_info {

    const unsigned int gGMS_HSW_METRICS_PROCESSING_MAJOR = 1;
    const unsigned int gGMS_HSW_METRICS_PROCESSING_MINOR = 0;
    const unsigned int gGMS_HSW_METRICS_PROCESSING_MICRO = 0;
    const unsigned int gGMS_HSW_METRICS_PROCESSING_FULLVER =
         1000U*gGMS_HSW_METRICS_PROCESSING_MAJOR+
	 100U*gGMS_HSW_METRICS_PROCESSING_MINOR+
	 10U*gGMS_HSW_METRICS_PROCESSING_MICRO;
    const char * const pgGMS_HSW_METRICS_PROCESSING_CREATION_DATE = "15-08-2020 14:40 +00200 (15 AUG 2020 2:40PM GMT+2)";
    const char * const pgGMS_HSW_METRICS_PROCESSING_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_HSW_METRICS_PROCESSING_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_HSW_METRICS_PROCESSING_SYNOPSIS      = "Haswell TMA metrics processing (statistical mainly)";
}


#include <cstindt>



namespace gms {

        namespace  system {

	        void hsw_fetched_uops_data(const int64_t * __restrict,
		                           const int64_t * __restrict,
					   const int64_t * __restrict,
					   const int64_t * __restrict,
					   int64_t * __restrict,
					   const int32_t);

		void hsw_recovery_cycles_data(const int64_t * __restrict,
		                              const int64_t * __restrict,
					      int64_t * __restrict,
					      const int32_t,
					      const bool);
					      

		void hsw_execute_cycles_data(const int64_t * __restrict,
					     int64_t * __restrict,
					     const int32_t,
		                             const bool);
					     

		void hsw_sq_full_cycles_data(const int64_t * __restrict,
					     int64_t * __restrict,
					     const int32_t,
		                             const bool);
					    

		void hsw_itlb_miss_cycles_data(const int64_t * __restrict,
		                               const int64_t * __restrict,
					       int64_t * __restrict,
					       const int32_t);

		void hsw_frontend_rs_empty_cycles_data(const int64_t * __restrict,
		                                       const float * __restrict,
						       int64_t * __restrict,
						       const int32_t);

		void hsw_cycles_0_ports_utilized_data(const int64_t * __restrict,
		                                      const int64_t * __restrict,
						      const int64_t * __restrict,
						      const float * __restrict,
						      int64_t * __restrict,
						      const int32_t,
						      const bool);
						     

		void hsw_cycles_1_ports_utilized_data(const int64_t * __restrict,
		                                      const int64_t * __restrict,
						      int64_t * __restrict,
						      const int32_t,
						      const bool);

	        void hsw_cycles_2_ports_utilized_data(const int64_t * __restrict,
						      const int64_t * __restrict,
						      int64_t * __restrict,
						      const int32_t,
						      const bool);

		void hsw_cycles_3_ports_utilized_data(const int64_t * __restrict,
						      int64_t * __restrict,
						      const int32_t,
						      const bool);

		void hsw_frontend_lat_cycles_data(const int64_t * __restrict,
						  const int64_t * __restrict,
						  int64_t * __restrict,
						  const int32_t);

		void hsw_stalls_mem_any_data(const int64_t * __restrict,
					     const int64_t * __restrict,
					     int64_t * __restrict,
					     const int32_t);

		void hsw_stalls_total_data(const int64_t * __restrict,
					   const int64_t * __restrict,
					   int64_t * __restrict,
					   const int32_t);

		void hsw_oro_drd_any_cycles_data(const int64_t * __restrict,
						 const int64_t * __restrict,
						 int64_t * __restrict,
						 const int32_t);

		void hsw_oro_drd_bw_cycles_data(const int64_t * __restrict,
						 const int64_t * __restrict,
						 int64_t * __restrict,
						 const int32_t);

		void hsw_oro_demand_rfo_c1(const int64_t * __restrict,
					   const int64_t * __restrict,
					   int64_t * __restrict,
					   const int32_t);

		void hsw_store_l2_hit_cycles_data(const int64_t * __restrict,
						 const int64_t * __restrict,
						 float * __restrict,
						 const int32_t);

		void hsw_load_l1_miss_data(const int64_t * __restrict,
					   const int64_t * __restrict,
					   const int64_t * __restrict,
					   const int64_t * __restrict,
					   const int64_t * __restrict,
					   int64_t * __restrict,
					   const int32_t);

		void hsw_load_l1_miss_net_data(const int64_t * __restrict,
					       const int64_t * __restrict,
					       const int64_t * __restrict,
					       const int64_t * __restrict,
					       const int64_t * __restrict,
					       const int64_t * __restrict,
					       int64_t * __restrict,
					       const int32_t);

		void hsw_load_l3_hit_data(const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  float * __restrict,
					  const int32_t);

		void hsw_load_xsnp_hit_data(const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  float * __restrict,
					  const int32_t);

		void hsw_load_xsnp_hitm_data(const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  float * __restrict,
					  const int32_t);

		void hsw_load_xsnp_miss_data(const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  const int64_t * __restrict,
					  float * __restrict,
					  const int32_t);

		void hsw_few_uops_exec_thresh_data(const int64_t * __restrict,
						   const int64_t * __restrict,
						   const float * __restrict,
						   float * __restrict,
						   const int32_t);

		void hsw_backend_bound_cycles_data(const int64_t * __restrict,
					           const int64_t * __restrict,
					           const int64_t * __restrict,
					           const int64_t * __restrict,
					           const int64_t * __restrict,
						   float * __restrict,
						   const int32_t,
						   const bool);

		void hsw_mem_bound_fraction_data(const int64_t * __restrict,
						 const int64_t * __restrict,
						 const float * __restrict,
						 float * __restrict,
						 const int32_t);

		void hsw_mem_l3_hit_frac_data(const int64_t * __restrict,
					      const int64_t * __restrict,
					      float * __restrict,
					      const int32_t);

		void hsw_mem_lock_st_frac_data(const int64_t * __restrict,
					      const int64_t * __restrict,
					      float * __restrict,
					      const int32_t);

		void hsw_mispred_clear_frac_data(const int64_t * __restrict,
					         const int64_t * __restrict,
					         float * __restrict,
					         const int32_t);

		void hsw_retire_fraction_data(const int64_t * __restrict,
					      const int64_t * __restrict,
					      float * __restrict,
					      const int32_t);

		void hsw_ipc_data(const int64_t * __restrict,
				  const int64_t * __restrict,
				  float * __restrict,
				  const int32_t);

		void hsw_upi_data(const int64_t * __restrict,
				  const int64_t * __restrict,
				  float * __restrict,
				  const int32_t);

		void hsw_iptb_data(const int64_t * __restrict,
				   const int64_t * __restrict,
				   float * __restrict,
				   const int32_t);

		void hsw_cpi_data(const int64_t * __restrict,
				  const int64_t * __restrict,
				  float * __restrict,
				  const int32_t);

		void hsw_issue_slots_data(const int64_t * __restrict,
					  int64_t * __restrict,
					  const int32_t);

		void hsw_ipload_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_ipstore_data(const int64_t * __restrict,
				      const int64_t * __restrict,
				      float * __restrict,
				      const int32_t);

		void hsw_ipbranch_data(const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_ipcall_data( const int64_t * __restrict,
				      const int64_t * __restrict,
				      float * __restrict,
				      const int32_t);

		void hsw_biptb_data(  const int64_t * __restrict,
				      const int64_t * __restrict,
				      float * __restrict,
				      const int32_t);

		void hsw_dsb_coverage_data( const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_ipbaclear_data(  const int64_t * __restrict,
				          const int64_t * __restrict,
				          float * __restrict,
				          const int32_t);

		void hsw_ipc_core_data(   const int64_t * __restrict,
				          const int64_t * __restrict,
				          float * __restrict,
				          const int32_t);

		void hsw_ilp_data( const int64_t * __restrict,
				   const int64_t * __restrict,
				   float * __restrict,
				   const int32_t,
				   bool);

		void hsw_ip_mispredict_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_core_clks_data(const int64_t * __restrict,
					const int64_t * __restrict,
					const int64_t * __restrict,
					const int64_t * __restrict,
					int64_t * __restrict,
					const int32_t,
					const bool);

		void hsw_load_miss_real_lat_data(const int64_t * __restrict,
					         const int64_t * __restrict,
					         const int64_t * __restrict,
						 float * __restrict,
						 const int32_t);

		void hsw_mem_level_parallel_data(const int64_t * __restrict,
				                 const int64_t * __restrict,
				                 float * __restrict,
						 const int32_t);

		void hsw_page_walk_util_data(const int64_t * __restrict,
					     const int64_t * __restrict,
					     const int64_t * __restrict,
					     const int64_t * __restrict,
					     float * __restrict,
					     const int32_t);

		void hsw_l1d_cache_fill_bw_data(const int64_t * __restrict,
				                 const int64_t * __restrict,
				                 float * __restrict,
						 const int32_t);

		void hsw_l2_cache_fill_bw_data(const int64_t * __restrict,
				                 const int64_t * __restrict,
				                 float * __restrict,
						 const int32_t);

		void hsw_l3_cache_fill_bw_data(const int64_t * __restrict,
				                 const int64_t * __restrict,
				                 float * __restrict,
						 const int32_t);

		void hsw_l1mpki_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_l2mpki_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_l2hpki_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_l3mpki_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_ipfabr_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_assist_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_kernel_util_data(const int64_t * __restrict,
				          const int64_t * __restrict,
				          float * __restrict,
					  const int32_t);

		void hsw_dram_bw_use_data(const int64_t * __restrict,
				          const int64_t * __restrict,
					  const int64_t * __restrict,
				          float * __restrict,
					  const int32_t);

		void hsw_mem_rqsts_lat_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
					     const int32_t);

		void hsw_mem_parallel_rqsts_data(const int64_t * __restrict,
				                 const int64_t * __restrict,
				                 float * __restrict,
					         const int32_t);

		void hsw_ht_utilization_data(const int64_t * __restrict,
				             const int64_t * __restrict,
				             float * __restrict,
					     const int32_t,
					     const bool);

		void hsw_frontend_bound_data(const int64_t * __restrict,
				             const int64_t * __restrict,
				             float * __restrict,
					     const int32_t);

		void hsw_frontend_latency_data(const int64_t * __restrict,
				               const int64_t * __restrict,
				               float * __restrict,
					       const int32_t);

		void hsw_icache_misses_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
					    const int32_t);

		void hsw_dsb_switches_data(const int64_t * __restrict,
				           const int64_t * __restrict,
				           float * __restrict,
					   const int32_t);

		void hsw_lcp_data(const int64_t * __restrict,
				  const int64_t * __restrict,
				  float * __restrict,
				  const int32_t);

		void hsw_ms_switches_data(const int64_t * __restrict,
				           const int64_t * __restrict,
				           float * __restrict,
					   const int32_t);

		void hsw_branch_resteers_data(const int64_t * __restrict,
					     const int64_t * __restrict,
					     const int64_t * __restrict,
					     const int64_t * __restrict,
					     float * __restrict,
					     const int32_t);

		void hsw_mite_data(const int64_t * __restrict,
				   const int64_t * __restrict,
				   const int64_t * __restrict,
				   float * __restrict,
				   const int32_t);

		void hsw_dsb_data(const int64_t * __restrict,
				   const int64_t * __restrict,
				   const int64_t * __restrict,
				   float * __restrict,
				   const int32_t);

		void hsw_l1_bound_data(const int64_t * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_dtlb_load_data(const int64_t * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_store_fwd_blk_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_split_loads_data(const int64_t * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_single_mul_clks_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_single_mul_core_clks_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_single_mul_uops_any_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_single_mul_uops_ret_any_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_simd_mov_elim_not_elim_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_int_mov_elim_not_elim_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_uops_issued_any_mite_uops_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_single_mul_avx_inst_all_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

		void hsw_frontend_latency_data(const int64_t * __restrict,
				            const int64_t * __restrict,
				            float * __restrict,
				            const int32_t);

	        void hsw_frontend_bw_data(const float * __restrict,
				     const float * __restrict,
				     float * __restrict,
				     const int32_t);

		
                void hsw_mite_data(const int64_t * __restrict,
			      const int64_t * __restrict,
			      const int64_t * __restrict,
			      float * __restrict,
			      const int32_t);

		void hsw_store_fwd_blocked_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_lock_latency_data(const float * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_split_loads_data(const float * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_4k_aliasing_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

	        void hsw_fb_full_data(const float * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_l2_bound_data(const int64_t * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_l3_bound_data(const float * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_contested_accesses_data(const float * __restrict,
				       const float * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_data_sharing_data(const float * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_hsw_dram_bound_data(const float * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_mem_bw_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_mem_latency_data(const int64_t * __restrict,
		                          const int64_t * __restrict,
					  const float * __restrict,
					   float * __restrict,
				          const int32_t);

		void hsw_store_bound_data(const int64_t * __restrict,
				     const int64_t * __restrict,
				     float * __restrict,
				     const int32_t);

		void hsw_dtlb_bound_data(const int64_t * __restrict,
				       const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_l3_hit_latency_data(const float * __restrict,
		                             const int64_t * __restrict,
					      float * __restrict,
				             const int32_t);

		void hsw_false_sharing_data( const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		void hsw_split_stores_data( const int64_t * __restrict,
				       const int64_t * __restrict,
				       float * __restrict,
				       const int32_t);

		

		
						     
		

    } // system

} // gms




/*__GMS_HSW_METRICS_PROCESSING_H__*/
