

#include "GMS_config.h"
#include "GSM_hsw_tma_metrics.h"

void
gms::system
::hsw_fetched_uops_data(const int64_t * __restrict idq_dsb_uops,
		        const int64_t * __restrict lsd_uops,
			const int64_t * __restrict idq_mite_uops,
			const int64_t * __restrict idq_ms_uops,
		        int64_t * __restrict samples,
			const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
          __assume_aligned(idq_dsb_uops,64);
	  __assume_aligned(lsd_uops,64);
	  __assume_aligned(idq_mite_uops,64);
          __assume_aligned(idq_ms_uops);
	  __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
          for(int32_t i = 0; i != len; ++i) {
              samples[i] = hsw_fetched_uops_data(idq_dsb_uops[i],
	                                         lsd_uops[i],
						 idq_mite_uops[i],
						 idq_ms_uops[i]);
	  }
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a  = (const int64_t*)__builtin_assume_aligned(idq_dsb_uops,64);
	 const int64_t * __restrict b  = (const int64_t*)__builtin_assume_aligned(lsd_uops,64);
	 const int64_t * __restrict c  = (const int64_t*)__builtin_assume_aligned(idq_mite_uops,64);
	 const int64_t * __restrict d  = (const int64_t*)__builtin_assume_aligned(idq_ms_uops,64);
	 int64_t *       __restrict s  = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
          for(int32_t i = 0; i != len; ++i) {
             s[i] = hsw_fetched_uops(a[i],
	                             b[i],
				     c[i],
				     d[i]);
	 }
#endif

}   



void
gms::system
::hsw_recovery_cycles_data(const int64_t * __restrict int_misc_recovery_cycles_any,
		           const int64_t * __restrict int_misc_recovery_cycles,
			   int64_t * __restrict samples,
			   const int32_t len,
			   const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(int_misc_recovery_cycles_any,64);
	 __assume_aligned(int_misc_recovery_cycles,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
          for(int32_t i = 0; i != len; ++i) {
              samples[i] = hsw_recovery_cycles( int_misc_recovery_cycles_any[i],
	                                         int_misc_recovery_cycles[i],
						 is_ht_enabled);
	  }
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(int_misc_recovery_cycles_any,64);
	 const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(int_misc_recovery_cycles,64);
	 int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
         for(int32_t i = 0; i != len; ++i) {
             s[i] = hsw_recovery_cycles(a[i],
	                                b[i],
					is_ht_enabled);
	 }
#endif

}
       
void
gms::system
::hsw_execute_cycles_data(const int64_t * __restrict uops_executed_core_c1,
                          int64_t * __restrict samples,
			  const int32_t len,
			  const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(uops_executed_core_c1,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_execute_cycles(uops_executed_core_c[i],
	                                    is_ht_enabled);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c1,64);
	 int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
          for(int32_t i = 0; i != len; ++i) {
              s[i] = hsw_execute_cycles(a[i],
	                                is_ht_enabled);
       }
#endif

}

      


void
gms_system
::hsw_sq_full_cycles_data(const int64_t * __restrict offcore_requests_buffer_sq_full,
                          int64_t * __restrict samples,
			  const int32_t len,
			  const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(offcore_requests_buffer_sq_full,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_sq_full_cycles(offcore_requests_buffer_sq_full[i],
	                                     is_ht_enabled);
	 }
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(offcore_requests_buffer_sq_full,64);
	 int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
         for(int32_t i = 0; i != len; ++i) {
             s[i] = hsw_sq_full_cycles(a[i],
	                               is_ht_enabled);
	 }
#endif

}


void
gms::system
::hsw_itlb_miss_cycles_data(const int64_t * __restrict  itlb_misses_stlb_hit ,
		            const int64_t * __restrict  itlb_misses_walk_duration,
			    int64_t * __restrict samples,
			    const int32_t len ) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(itlb_misses_stlb_hit,64);
	 __assume_aligned(itlb_misses_walk_duration,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_itlb_miss_cycles(itlb_misses_stlb_hit[i],
	                                      itlb_misses_walk_duration[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(itlb_misses_stlb_hit,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(itlb_misses_walk_duration,64)
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_itlb_misses_cycles(a[i],
	                                  b[i]);
	}
#endif
}


void
gms::system
::hsw_frontend_rs_empty_cycles_data(const int64_t * __restrict rs_event_empty_cycles ,
		                    const double * __restrict frontend_latency ,
				    int64_t * __restrict samples,
				    const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(rs_event_empty_cycles,64);
	 __assume_aligned(frontend_latency,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_frontend_rs_empty_cycles(rs_event_empty_cycles[i],
	                                               frontend_latency[i]);
	 }
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(rs_event_empty_cycles,64);
	 const double  * __restrict b = (const double*)__builtin_assume_aligned(frontend_latency,64);
	 int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
         for(int32_t i = 0; i != len; ++i) {
             s[i] = hsw_frontend_rs_empty_cycles(a[i],
	                                         b[i]);
	 }
#endif
}


void
gms::system
::hsw_cycles_0_ports_utilized_data(const int64_t * __restrict uops_executed_core_i1_c1,
		                   const int64_t * __restrict stalls_total,
				   const int64_t * __restrict rs_event_empty_cycles,
				   const double * __restrict  frontend_latency,
				   int64_t * __restrict samples,
				   const int32_t len,
				   const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(uops_executed_core_i1_c1,64);
	 __assume_aligned(stalls_total,64);
	 __assume_aligned(rs_event_empty_cycles,64);
	 __assume_aligned(frontend_latency,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_cycles_0_ports_utilized(uops_executed_core_i1_c1[i],
	                                              stalls_total[i],
						      rs_event_empty_cycles[i],
						      frontend_latency[i],
						      is_ht_enabled);
	 }
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
       const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_executed_core_i1_c1,64);
       const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(stalls_total,64);
       const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(rs_event_empty_cycles,64);
       const double  * __restrict d = (const double*)__builtin_assume_aligned(frontend_latency,64);
       int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
       for(int32_t i = 0; i != len; ++i) {
           s[i] = hsw_cycles_0_ports_utilized(a[i],
	                                      b[i],
					      c[i],
					      d[i],
					      is_ht_enabled);
       }
#endif
}


void
gms::system
::hsw_cycles_1_ports_utilized_data(const int64_t * __restrict uops_executed_core_c1,
		                   const int64_t * __restrict uops_executed_core_c2,
				   int64_t * __restrict samples,
				   const int32_t len,
				   const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(uops_executed_core_c1,64);
	 __assume_aligned(uops_executed_core_c2,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_cycles_1_ports_utilized(uops_executed_core_c1[i],
	                                              uops_executed_core_c2[i],
						      is_ht_enabled);
	 }
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
       const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c1,64);
       const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c2,64);
       int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
       for(int32_t i = 0; i != len; ++i) {
           s[i] = hsw_cycles_1_ports_utilized(a[i],
	                                      b[i],
					      is_ht_enabled);
       }
#endif
}


void
gms::system
::hsw_cycles_2_ports_utilized_data(const int64_t * __restrict uops_executed_core_c2,
		                   const int64_t * __restrict uops_executed_core_c3,
				   int64_t * __restrict samples,
				   const int32_t len,
				   const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(uops_executed_core_c2,64);
	 __assume_aligned(uops_executed_core_c3,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_cycles_2_ports_utilized(uops_executed_core_c2[i],
	                                              uops_executed_core_c3[i],
						      is_ht_enabled);
	 }
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
       const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c2,64);
       const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c3,64);
       int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
       for(int32_t i = 0; i != len; ++i) {
           s[i] = hsw_cycles_2_ports_utilized(a[i],
	                                      b[i],
					      is_ht_enabled);
       }
#endif
}


void
gms::system
::hsw_cycles_3_ports_utilized(const int64_t * __restrict uops_executed_core_c3,
                              int64_t * __restrict samples,
			      const int32_t len,
			      const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(uops_executed_core_c3,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_cycles_3_ports_utilized(uops_executed_core_c3[i],
	                                              is_ht_enabled);
	 }
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c3,64);
	 int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
         for(int32_t i = 0; i != len; ++i) {
             s[i] = hsw_cycles_3_ports_utilized(a[i],
	                                        is_ht_enabled);
	 }
#endif
}


void
gms::system
::hsw_frontend_lat_cycles_data(const int64_t * __restrict cpu_clk_unhalted_thread,
			       const int64_t * __restrict idq_uops_not_delivered_cycles_0_uops_deliv_core ,
			       int64_t * __restrict samples,
			       const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(cpu_clk_unhalted_thread,64);
	 __assume_aligned(idq_uops_not_delivered_cycles_0_uops_deliv_core,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_frontend_latency_cycles(cpu_clk_unhalted_thread[i],
	                                              idq_uops_not_delivered_cycles_0_uops_deliv_core[i]);
	                                              
	 }
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(cpu_clk_unhalted_thread,64);
	 const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(idq_uops_not_delivered_cycles_0_uops_deliv_core,64);
	 int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
         for(int32_t i = 0; i != len; ++i) {
             s[i] = hsw_frontend_latency_cycles(a[i],
	                                        b[i]);
	 }
#endif
}


void
gms::system
::hsw_stalls_mem_any_data(const int64_t * __restrict cpu_clk_unhalted_thread,
			  const int64_t * __restrict cycles_activity_stalls_lm_pending,
			  int64_t * __restrict samples,
			  const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(cpu_clk_unhalted_thread,64);
	__assume_aligned(cycles_activity_stalls_lm_pending,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_stalls_mem_any(cpu_clk_unhalted_thread[i],
	                                    cycles_activity_stalls_lm_pending[i]);
	}
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(cpu_clk_unhalted_thread,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(cycles_activity_stalls_lm_pending,64);
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_stalls_mem_any(a[i],
	                              b[i]);
	}
#endif
}



void
gms::system
::hsw_stalls_total_data(const int64_t * __restrict cpu_clk_unhalted_thread,
			  const int64_t * __restrict  cycles_activity_cycles_no_execute,
			  int64_t * __restrict samples,
			  const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(cpu_clk_unhalted_thread,64);
	__assume_aligned(cycles_activity_cycles_no_execute,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_stalls_mem_any(cpu_clk_unhalted_thread[i],
	                                    cycles_activity_cycles_no_execute[i]);
	}
#elif  defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(cpu_clk_unhalted_thread,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(cycles_activity_cycles_no_execute,64);
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_stalls_mem_any(a[i],
	                              b[i]);
	}
#endif
}


void
gms::system
::hsw_oro_drd_any_cycles_data(const int64_t * __restrict cpu_clk_unhalted_thread ,
			      const int64_t * __restrict offcore_requests_oustanding_cycles_with_data_rd,
			      int64_t * __restrict samples,
			      const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(cpu_clk_unhalted_thread,64);
	__assume_aligned(offcore_requests_oustanding_cycles_with_data_rd,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_oro_drd_any_cycles(cpu_clk_unhalted_thread[i],
	                                        offcore_requests_oustanding_cycles_with_data_rd[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(cpu_clk_unhalted_thread,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned( offcore_requests_oustanding_cycles_with_data_rd,64);
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_oro_drd_any_cycles(a[i],
	                                  b[i]);
	}
#endif
}


void
gms::system
::hsw_oro_drd_bw_cycles_data(const int64_t * __restrict cpu_clk_unhalted_thread ,
			      const int64_t * __restrict offcore_requests_oustanding_all_data_rd_c6,
			      int64_t * __restrict samples,
			      const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(cpu_clk_unhalted_thread,64);
	__assume_aligned(offcore_requests_oustanding_all_data_rd_c6,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_oro_drd_bw_cycles(cpu_clk_unhalted_thread[i],
	                                        offcore_requests_oustanding_all_data_rd_c6[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(cpu_clk_unhalted_thread,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned( offcore_requests_oustanding_all_data_rd_c6,64);
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_oro_drd_bw_cycles(a[i],
	                                  b[i]);
	}
#endif
}


void
gms::system
::hsw_oro_demand_rfo_c1_data(const int64_t * __restrict cpu_clk_unhalted_thread ,
			      const int64_t * __restrict offcore_requests_oustanding_cycles_with_demand_rfo,
			      int64_t * __restrict samples,
			      const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(cpu_clk_unhalted_thread,64);
	__assume_aligned(offcore_requests_oustanding_cycles_with_demand_rfo,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_oro_demand_rf_c1(cpu_clk_unhalted_thread[i],
	                                      offcore_requests_oustanding_cycles_with_demand_rfo[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(cpu_clk_unhalted_thread,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned( offcore_requests_oustanding_cycles_with_demand_rfo,64);
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_oro_drd_bw_cycles(a[i],
	                                  b[i]);
	}
#endif
}


void
gms::system
::hsw_store_l2_hit_cycles_data(const int64_t * __restrict l2_rqsts_rfo_hit ,
			       const int64_t * __restrict mem_uops_retired_lock_loads,
			       const int64_t * __restrict mem_uops_retired_all_stores,
			       double * __restrict samples,
			       const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(l2_rqsts_rfo_hit,64);
	__assume_aligned(mem_uops_retired_lock_loads,64);
	__assume_aligned(mem_uops_retired_all_stores,64)
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_store_l2_hit_cycles(l2_rqsts_rfo_hit[i],
	                                         mem_uops_retired_lock_loads[i],
						 mem_uops_retired_all_stores[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(l2_rqsts_rfo_hit,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_uops_retired_lock_loads,64);
	const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(mem_uops_retired_all_stores,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_store_l2_hit_cycles(a[i],
	                                   b[i],
					   c[i]);
	}
#endif
}


void
gms::system
::hsw_load_l1_miss_data(const int64_t * __restrict mem_load_uops_retired_l2_hit,
			const int64_t * __restrict mem_load_uops_retired_l3_hit,
			const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hit ,
			const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hitm,
		        const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_miss,
		        int64_t * __restrict samples,
			const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_load_uops_retired_l2_hit,64);
	__assume_aligned(mem_load_uops_retired_l3_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_load_l1_miss(mem_load_uops_retired_l2_hit[i],
	                                  mem_load_uops_retired_l3_hit[i],
					  mem_load_uops_l3_hit_retired_xsnp_hit[i],
					  mem_load_uops_l3_hit_retired_xsnp_hitm[i],
					  mem_load_uops_l3_hit_retired_xsnp_miss[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l2_hit,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_hit,64);
	const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	const int64_t * __restrict d = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	const int64_t * __restrict e = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_load_l1_miss(a[i],
	                            b[i],
				    c[i],
				    d[i],
				    e[i]);
	}
#endif
}
