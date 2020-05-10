

#ifndef __GMS_HSW_TMA_METRICS_API_H__
#define __GMS_HSW_TMA_METRICS_API_H__ 100520200913


// This implementation must be compiled as C version only.
// This is an API for Fortran callers.

// File info

  const unsigned int gGMS_HSW_CLIENT_METRICS_MAJOR = 1;
  const unsigned int gGMS_HSW_CLIENT_METRICS_MINOR = 0;
  const unsigned int gGMS_HSW_CLIENT_METRICS_MICRO = 0;
  const unsigned int gGMS_HSW_CLIENT_METRICS_FULLVER =
      1000U*gGMS_HSW_CLIENT_METRICS_MAJOR+
      100U*gGMS_HSW_CLIENT_METRICS_MINOR+
      10U*gGMS_HSW_CLIENT_METRICS_MICRO;
  const char * const pgGMS_HSW_CLIENT_METRICS_CREATION_DATE = "10-05-2020 09:15 +00200 (SUN 10 MAY 2020 09:36AM GMT+2)";
  const char * const pgGMS_HSW_CLIENT_METRICS_BUILD_DATE    = __DATE__ ":" __TIME__;
  const char * const pgGMS_HSW_CLIENT_METRICS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_HSW_CLIENT_METRICS_SYNOPSIS      = "Haswell client performance metrics based on TMA-Metrics (4.0)";


#include <stdint.h>

// Constants

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

    
    uint64_t hsw_fetched_uops(const uint64_t,
                              const uint64_t ,
			      const uint64_t ,
			      const uint64_t );
#include <stdbool.h>



    uint64_t hsw_recovery_cycles(const uint64_t ,
                                 const uint64_t ,
			         const bool);

    uint64_t hsw_execute_cycles(const uint64_t ,
                                const bool);

    uint64_t hsw_sq_full_cycles(const uint64_t,
                                const bool);

    uint64_t hsw_itlb_miss_cycles(const uint64_t,
                                  const uint64_t);

    uint64_t hsw_frontend_rs_empty_cycles(const uint64_t,
                                          const float);

    uint64_t hsw_cycles_0_ports_utilized(const uint64_t ,
                                         const uint64_t,
					 const uint64_t,
					 const float,
					 const bool);

    uint64_t hsw_cycles_1_ports_utilized(const uint64_t,
                                         const uint64_t ,
				         const bool);

    uint64_t hsw_cycles_2_ports_utilized(const uint64_t,
                                         const uint64_t,
					 const bool);

    uint64_t hsw_cycles_3_ports_utilized(const uint64_t ,
                                         const bool);

    uint64_t hsw_frontend_latency_cycles(const uint64_t,
                                         const uint64_t);

    uint64_t hsw_stalls_mem_any(const uint64_t,
                                const uint64_t);

    uint64_t hsw_stalls_total(const uint64_t,
                              const uint64_t);

    uint64_t  hsw_oro_drd_any_cycles(const uint64_t,
                                    const uint64_t);

    uint64_t hsw_oro_drd_bw_cycles(const uint64_t,
                                   const uint64_t);

    uint64_t hsw_oro_demand_rfo_c1(const uint64_t,
                                   const uint64_t);

    float hsw_store_l2_hit_cycles(const uint64_t,
                                  const uint64_t ,
				  const uint64_t);

    uint64_t hsw_load_l1_miss(const uint64_t,
                              const uint64_t,
			      const uint64_t,
			      const uint64_t,
			      const uint64_t);

    uint64_t hsw_load_l1_miss_net(const uint64_t,
                                  const uint64_t,
                                  const uint64_t,
			          const uint64_t,
			          const uint64_t,
			          const uint64_t);

    float hsw_load_l3_hit(const uint64_t,
                          const uint64_t,
			  const uint64_t,
			  const uint64_t,
			  const uint64_t,
			  const uint64_t);

    float hsw_load_xsnp_hit(const uint64_t,
                            const uint64_t,
			    const uint64_t,
			    const uint64_t,
			    const uint64_t,
			    const uint64_t);

    float hsw_load_xsnp_hitm(
                             const uint64_t,
			     const uint64_t,
	                     const uint64_t,
			     const uint64_t,
			     const uint64_t,
			     const uint64_t);

     float hsw_load_xsnp_miss( const uint64_t ,
			     const uint64_t,
	                     const uint64_t,
			     const uint64_t,
			     const uint64_t,
			     const uint64_t);

      uint64_t hsw_few_uops_executed_threshold(const uint64_t,
                                             const uint64_t,
					     const float);

      float hsw_backend_bound_cycles(const uint64_t,
                                    const uint64_t,
				    const uint64_t,
				    const uint64_t,
				    const uint64_t,
				    const bool);

       float hsw_memory_bound_fraction(const uint64_t,
                                       const uint64_t,
				       const float);

       float hsw_mem_l3_hit_fraction( const uint64_t,
                                      const uint64_t);

       float hsw_mem_lock_st_fraction( const uint64_t,
                                       const uint64_t);

       float hsw_mispred_clears_fraction( const uint64_t,
                                          const uint64_t);

       float hsw_retire_fraction( const uint64_t,
                                  const uint64_t);

       float hsw_ipc( const uint64_t,
                      const uint64_t);

       float hsw_upi( const uint64_t,
                      const uint64_t);

       float hsw_iptb( const uint64_t,
                       const uint64_t);

       float hsw_cpi( const uint64_t,
                      const uint64_t);

       uint64_t hsw_issue_slots( const uint64_t);

       float hsw_ipload( const uint64_t,
                       const uint64_t);

       float hsw_ipstore( const uint64_t,
                        const uint64_t);

       float hsw_ipbranch( const uint64_t,
                         const uint64_t);

       float hsw_ipcall( const uint64_t ,
                         const uint64_t);

       float hsw_biptb( const uint64_t,
                        const uint64_t);

       float hsw_dsb_coverage( const uint64_t,
                             const uint64_t);

       float hsw_ipbaclear( const uint64_t,
                            const uint64_t);

       float hsw_ipc_core( const uint64_t,
                         const uint64_t);

       float hsw_ilp( const uint64_t,
                      const uint64_t,
		      const bool);

       float hsw_ip_mispredict( const uint64_t,
                              const uint64_t);

       uint64_t hsw_core_clks( const uint64_t,
                          const uint64_t,
			  const uint64_t,
			  const uint64_t,
			  const bool);

       float hsw_load_miss_real_latency(const uint64_t,
                                      const uint64_t,
				      const uint64_t);

       float hsw_mem_level_parallelism( const uint64_t,
                                      const uint64_t);

       float hsw_page_walk_util( const uint64_t,
                               const uint64_t,
			       const uint64_t,
			       const uint64_t);

       float hsw_l1d_cache_fill_bw( const uint64_t,
                                  const uint64_t);

       float hsw_l2_cache_fill_bw( const uint64_t ,
                                 const uint64_t);

       float hsw_l3_cache_fill_bw( const uint64_t,
                                 const uint64_t);

       float hsw_l1mpki( const uint64_t,
                         const uint64_t);

       float hsw_l2mpki( const uint64_t,
                       const uint64_t);

       float hsw_l2hpki(  const uint64_t,
                          const uint64_t);

       float hsw_l3mpki( const uint64_t ,
                         const uint64_t);

       float hsw_ipfarbr( const uint64_t ,
                        const uint64_t);

        float hsw_assist( const uint64_t,
                       const uint64_t);

		       


#endif /*__GMS_HSW_TMA_METRICS_API_H__*/
