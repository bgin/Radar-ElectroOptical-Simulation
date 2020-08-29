

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


void
gms::system
::hsw_load_l1_miss_net_data(const int64_t * __restrict mem_load_uops_retired_l3_miss,
			    const int64_t * __restrict mem_load_uops_retired_l2_hit,
			    const int64_t * __restrict mem_load_uops_retired_l3_hit,
			    const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hit,
			    const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hitm,
			    const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_miss,
			    int64_t * __restrict samples,
			    const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_load_uops_retired_l3_miss,64);
        __assume_aligned(mem_load_uops_retired_l2_hit,64);
	__assume_aligned(mem_load_uops_retired_l3_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_load_l1_miss_net(mem_load_uops_retired_l3_miss[i],
	                                      mem_load_uops_retired_l2_hit[i],
	                                      mem_load_uops_retired_l3_hit[i],
					      mem_load_uops_l3_hit_retired_xsnp_hit[i],
					      mem_load_uops_l3_hit_retired_xsnp_hitm[i],
					      mem_load_uops_l3_hit_retired_xsnp_miss[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_miss,64);
        const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l2_hit,64);
	const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_hit,64);
	const int64_t * __restrict d = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	const int64_t * __restrict e = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	const int64_t * __restrict f = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragms GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_load_l1_miss_net(a[i],
	                                b[i],
					c[i],
					d[i],
					e[i],
					f[i]);
	}
#endif
}


void
gms::system
::hsw_load_l3_hit_data(const int64_t * __restrict mem_load_uops_retired_l3_hit,
		       const int64_t * __restrict mem_load_uops_retired_hit_lfb,
		       const int64_t * __restrict mem_load_uops_retired_l2_hit,
		       const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hit,
		       const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hitm,
		       const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_miss,
		       double * __restrict samples,
		       const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_load_uops_retired_l3_miss,64);
        __assume_aligned(mem_load_uops_retired_hit_lfb,64);
	__assume_aligned(mem_load_uops_retired_l2_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_load_l3_hit(mem_load_uops_retired_l3_hit[i],
	                                 mem_load_uops_retired_hit_lfb[i],
	                                 mem_load_uops_retired_l2_hit[i],
					 mem_load_uops_l3_hit_retired_xsnp_hit[i],
					 mem_load_uops_l3_hit_retired_xsnp_hitm[i],
					 mem_load_uops_l3_hit_retired_xsnp_miss[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_hit,64);
        const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_hit_lfb,64);
	const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l2_hit,64);
	const int64_t * __restrict d = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	const int64_t * __restrict e = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	const int64_t * __restrict f = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	double* __restrict         s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_load_l3_hit(a[i],
	                           b[i],
				   c[i],
				   d[i],
				   e[i],
				   f[i]);
	}
#endif
        
}


void
gms::system
::hsw_load_xsnp_hit_data(const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hit,
			 const int64_t * __restrict mem_load_uops_retired_hit_lfb,
			 const int64_t * __restrict mem_load_uops_retired_l2_hit,
			 const int64_t * __restrict mem_load_uops_retired_l3_hit,
			 const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hitm,
			 const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hit,
			 double * __restrict samples,
			 const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
        __assume_aligned(mem_load_uops_retired_hit_lfb,64);
	__assume_aligned(mem_load_uops_retired_l2_hit,64);
	__assume_aligned(mem_load_uops_retired_l3_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_load_xsnp_hit(mem_load_uops_l3_hit_retired_xsnp_hit[i],
	                                 mem_load_uops_retired_hit_lfb[i],
	                                 mem_load_uops_retired_l2_hit[i],
					 mem_load_uops_retired_l3_hit[i],
					 mem_load_uops_l3_hit_retired_xsnp_hitm[i],
					 mem_load_uops_l3_hit_retired_xsnp_miss[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
        const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_hit_lfb,64);
	const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l2_hit,64);
	const int64_t * __restrict d = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_hit,64);
	const int64_t * __restrict e = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	const int64_t * __restrict f = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	double* __restrict         s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_load_xsnp_hit(a[i],
	                           b[i],
				   c[i],
				   d[i],
				   e[i],
				   f[i]);
	}
#endif
}


void
gms::system
::hsw_load_xsnp_hitm_data(const int64_t * __restrict mem_load_uops_retired_hit_lfb,
			  const int64_t * __restrict mem_load_uops_retired_l2_hit,
			  const int64_t * __restrict mem_load_uops_retired_l3_hit,
			  const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hit,
			  const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hitm,
			  const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_miss,
			  double * __restrict samples,
			  const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_load_uops_retired_hit_lfb,64);
        __assume_aligned(mem_load_uops_retired_l2_hit,64);
	__assume_aligned(mem_load_uops_retired_l3_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_load_xsnp_hitm(mem_load_uops_retired_hit_lfb[i],
	                                 mem_load_uops_retired_l2_hit[i],
	                                 mem_load_uops_retired_l3_hit[i],
					 mem_load_uops_l3_hit_retired_xsnp_hit[i],
					 mem_load_uops_l3_hit_retired_xsnp_hitm[i],
					 mem_load_uops_l3_hit_retired_xsnp_miss[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_hit_lfb,64);
        const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l2_hit,64);
	const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_hit,64);
	const int64_t * __restrict d = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	const int64_t * __restrict e = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	const int64_t * __restrict f = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	double* __restrict         s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_load_xsnp_hit(a[i],
	                           b[i],
				   c[i],
				   d[i],
				   e[i],
				   f[i]);
	}
#endif
}


void
gms::system
::hsw_load_xsnp_miss_data(const int64_t * __restrict mem_load_uops_retired_hit_lfb,
			  const int64_t * __restrict mem_load_uops_retired_l2_hit,
			  const int64_t * __restrict mem_load_uops_retired_l3_hit,
			  const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hit,
			  const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_hitm,
			  const int64_t * __restrict mem_load_uops_l3_hit_retired_xsnp_miss,
			  double * __restrict samples,
			  const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_load_uops_retired_hit_lfb,64);
        __assume_aligned(mem_load_uops_retired_l2_hit,64);
	__assume_aligned(mem_load_uops_retired_l3_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	__assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_load_xsnp_miss(mem_load_uops_retired_hit_lfb[i],
	                                 mem_load_uops_retired_l2_hit[i],
	                                 mem_load_uops_retired_l3_hit[i],
					 mem_load_uops_l3_hit_retired_xsnp_hit[i],
					 mem_load_uops_l3_hit_retired_xsnp_hitm[i],
					 mem_load_uops_l3_hit_retired_xsnp_miss[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_hit_lfb,64);
        const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l2_hit,64);
	const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_hit,64);
	const int64_t * __restrict d = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hit,64);
	const int64_t * __restrict e = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_hitm,64);
	const int64_t * __restrict f = (const int64_t*)__builtin_assume_aligned(mem_load_uops_l3_hit_retired_xsnp_miss,64);
	double* __restrict         s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_load_xsnp_miss(a[i],
	                           b[i],
				   c[i],
				   d[i],
				   e[i],
				   f[i]);
	}
#endif
}


void
gms::system
::hsw_few_uops_exec_thresh_data(const int64_t * __restrict uops_executed_core_c2,
                                const int64_t * __restrict uops_executed_core_c3,
				const double  * __restrict ipc,
				int64_t * __restrict samples,
				const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(uops_executed_core_c2,64);
	__assume_aligned(uops_executed_core_c3,64);
	__assume_aligned(ipc,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_few_uops_executed_threshold(uops_executed_core_c2[i],
	                                                 uops_executed_core_c3[i],
							 ipc[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
       const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c2,64);
       const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c3,64);
       const double  * __restrict c = (const double*)__builtin_assume_aligned(ipc,64);
       int64_t * __restrict       s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
       for(int32_t i = 0; i != len; ++i) {
           s[i] = hsw_few_uops_executed_threshold(a[i],
	                                          b[i],
						  c[i]);
       }
#endif
}


void
gms::system
::hsw_backend_bound_cycles_data(const int64_t * __restrict stalls_total,
                                const int64_t * __restrict uops_executed_core_c1,
				const int64_t * __restrict few_uops_executed_threshold,
				const int64_t * __restrict frontend_rs_empty_cycles,
				const int64_t * __restrict resource_stalls_sb,
				double * __restrict samples,
				const int32_t len,
				const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(stalls_total,64);
	 __assume_aligned(uops_executed_core_c1,64);
	 __assume_aligned(few_uops_executed_threshold,64);
	 __assume_aligned(frontend_rs_empty_cycles,64);
	 __assume_aligned(resource_stalls_sb,64);
	 __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
         for(int32_t i = 0; i != len; ++i) {
             samples[i] = hsw_backend_bound_cycles(stalls_total[i],
	                                           uops_executed_core_c1[i],
						   few_uops_executed_threshold[i],
						   frontend_rs_empty_cycles[i],
						   resource_stalls_sb[i],
						   is_ht_enabled);
	 }
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
         const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(stalls_total,64);
	 const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(uops_executed_core_c1,64);
	 const int64_t * __restrict c = (const int64_t*)__builtin_assume_aligned(frontend_rs_empty_cycles,64);
	 const int64_t * __restrict d = (const int64_t*)__builtin_assume_aligned(few_uops_executed_threshold,64);
	 const int64_t * __restrict e = (const int64_t*)__builtin_assume_aligned(resource_stalls_sb,64);
	 double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
         for(int32_t i = 0; i != len; ++i) {
             s[i] = hsw_backend_bound_cycles(a[i],
	                                     b[i],
					     c[i],
					     d[i],
					     e[i],
					     is_ht_enabled);
					     
	 }
#endif
}


void
gms::system
::hsw_mem_bound_fraction_data(const int64_t * __restrict stalls_mem_any,
                              const int64_t * __restrict resource_stalls_sb,
			      const double *  __restrict backend_bound_cycles,
			      double * __restrict samples,
			      const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
       __assume_aligned(stalls_mem_any,64);
       __assume_aligned(resource_stalls_sb,64);
       __assume_aligned(backend_bound_cycles,64);
       __assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
       for(int32_t i = 0; i != len; ++i) {
           samples[i] = hsw_mem_bound_fraction(stalls_mem_any[i],
	                                       resource_stalls_sb[i],
					       backend_bound_cycles[i]);
       }
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
       const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(stalls_mem_any,64);
       const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(resource_stalls_sb,64);
       const double  * __restrict c = (const double*)__builtin_assume_aligned(backend_bound_cycles,64);
       double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
       for(int32_t i = 0; i != len; ++i) {
           s[i] = hsw_mem_bound_fraction(a[i],
	                                 b[i],
					 c[i]);
       }
#endif
}


void
gms::system
::hsw_mem_l3_hit_frac_data(const int64_t * __restrict mem_load_uops_retired_l3_hit,
                           const int64_t * __restrict mem_load_uops_retired_l3_miss,
			   double * __restrict samples,
			   const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_load_uops_retired_l3_hit,64);
	__assume_aligned(mem_load_uops_retired_l3_miss,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_mem_l3_hit_fraction(mem_load_uops_retired_l3_hit[i],
	                                         mem_load_uops_retired_l3_miss[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_hit,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_load_uops_retired_l3_miss,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_mem_l3_hit_fraction(a[i],
	                                   b[i]);
	}
#endif
}


void
gms::system
::hsw_mem_lock_st_frac_data(const int64_t * __restrict mem_uops_retired_lock_loads,
                            const int64_t * __restrict mem_uops_retired_all_stores,
			    double * __restric samples,
			    const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(mem_uops_retired_lock_loads,64);
	__assume_aligned(mem_uops_retired_all_stores,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_mem_lock_st_fraction(mem_uops_retired_lock_loads[i],
	                                          mem_uops_retired_all_stores[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(mem_uops_retired_lock_loads,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_uops_retired_all_stores,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_mem_lock_st_fraction(a[i],
	                                    b[i]);
	}
#endif
}


void
gms::system
::hsw_mispred_clear_frac_data(const int64_t * __restrict br_mispred_all_branches,
                              const int64_t * __restrict machines_clear_counts,
			      double * __restrict samples,
			      const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(br_mispred_all_branches,64);
	__assume_aligned(machines_clear_count,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_mispred_clear_fraction(br_mispred_all_branches[i],
	                                            machines_clear_counts[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(br_mispred_all_branches,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(machines_clear_count,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_mispred_clear_fraction(a[i],
	                                      b[i]);
	}
#endif
}


void
gms::system
::hsw_retire_fraction_data(const int64_t * __restrict uops_retired_retire_slots,
                           const int64_t * __restrict uops_issued_any,
			   double * __restrict samples,
			   const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(uops_retired_retire_slots,64);
	__assume_aligned(uops_issued_any,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_retire_fraction(uops_retired_retire_slots[i],
	                                     uops_issued_any[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_retired_retire_slots,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(uops_issued_any,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_retire_fraction(a[i],
	                               b[i]);
	}
#endif
}


void
gms::system
::hsw_ipc_data(const int64_t * __restrict instr_retired_any,
               const int64_t * __restrict clks,
	       double * __restrict samples,
	       const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(clks,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ipc(instr_retired_any[i],
	                         clks[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(clks,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ipc(a[i],
	                   b[i]);
	}
#endif
}


void
gms::system
::hsw_upi_data(const int64_t * __restrict uops_retired_retire_slots,
               const int64_t * __restrict uops_retired_any,
	       double * __restrict samples,
	       const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(uops_retired_retire_slots,64);
	__assume_aligned(uops_retired_any,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_upi(uops_retired_retire_slots[i],
	                         uops_retired_any[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_retired_retire_slots,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(uops_retired_any,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_upi(a[i],
	                   b[i]);
	}
#endif
}


void
gms::system
::hsw_iptb_data(const int64_t * __restrict instr_retired_any,
                const int64_t * __restrict br_instr_retired_near_taken,
		double * __restrict samples,
		const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assumed_aligned(instr_retired_any,64);
	__assume_aligned(br_instr_retired_near_taken,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t  i = 0; i != len; ++i) {
            samples[i] = hsw_iptb(instr_retired_any[i],
	                          br_instr_retired_near_taken[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(br_instr_retired_near_taken,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_iptb(a[i],
	                    b[i]);
	}
#endif
}


void
gms::system
::hsw_cpi_data(const int64_t * __restrict instr_retired_any,
               const int64_t * __restrict clks,
	       double * __restrict samples,
	       const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(clks,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_cpi(instr_retired_any[i],
	                         clks[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(clks,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_cpi(a[i],
	                   b[i]);
	}
#endif
}


void
gms::system
::hsw_issue_slots_data(const int64_t * __restrict core_clks,
                       int64_t * __restrict samples,
		       const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(core_clks,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_issue_slots(core_clks[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
       const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(core_clks,64);
       int64_t * __restrict       s = (int64_t*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
       for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_issue_slots(a[i]);
       }
#endif
}


void
gms::system
::hsw_ipload_data(const int64_t * __restrict instr_retired_any,
                  const int64_t * __restrict mem_uops_retired_all_loads,
		  double * __restrict samples,
		  const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(mem_uops_retired_all_loads,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ipload(instr_retired_any[i],
	                            mem_uops_retired_all_loads[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_uops_retired_all_loads,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ipload(a[i],
	                      b[i]);
	}
#endif
}


void
gms::system
::hsw_ipstore_data(const int64_t * __restrict instr_retired_any,
                   const int64_t * __restrict mem_uops_retired_all_stores,
		   double * __restrict samples,
		   const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(mem_uops_retired_all_stores,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ipstore(instr_retired_any[i],
	                            mem_uops_retired_all_stores[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(mem_uops_retired_all_stores,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ipstore(a[i],
	                      b[i]);
	}
#endif
}


void
gms::system
::hsw_ipbranch_data(const int64_t * __restrict instr_retired_any,
                    const int64_t * __restrict br_instr_retired_all_branches,
		    double * __restrict samples,
		    const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(br_instr_retired_all_branches,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ipbranch(instr_retired_any[i],
	                              br_instr_retired_all_branches[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(br_instr_retired_all_branches,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ipbranch(a[i],
	                      b[i]);
	}
#endif
}


void
gms::system
::hsw_ipcall_data(const int64_t * __restrict instr_retired_any,
                  const int64_t * __restrict br_instr_retired_near_call,
		  double * __restrict samples,
		  const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(br_instr_retired_near_call,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ipcall(instr_retired_any[i],
	                              br_instr_retired_near_call[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(br_instr_retired_near_call,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ipcall(a[i],
	                      b[i]);
	}
#endif
}


void
gms::system
::hsw_biptb_data( const int64_t * __restrict br_instr_retired_all_branches,
                  const int64_t * __restrict br_instr_retired_near_taken,
		  double * __restrict samples,
		  const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(br_instr_retired_all_branches,64);
	__assume_aligned(br_instr_retired_near_call,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_biptb(instr_retired_any[i],
	                              br_instr_retired_near_call[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(br_instr_retired_near_call,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_biptb(a[i],
	                      b[i]);
	}
#endif
}


void
gms::system
::hsw_dsb_coverage_data( const int64_t * __restrict idq_dsb_uops,
                    const int64_t * __restrict fetched_uops,
		    double * __restrict samples,
		    const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(idq_dsb_uops,64);
	__assume_aligned(fetched_uops,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_dsb_coverage(idq_dsb_uops[i],
	                            fetched_uops[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(idq_dsb_uops,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(fetched_uops,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_dsb_coverage(a[i],
	                            b[i]);
	}
#endif
}


void
gms::system
::hsw_ipbaclear_data( const int64_t * __restrict instr_retired_any,
                      const int64_t * __restrict baclears_any,
		      double * __restrict samples,
		      const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(baclears_any,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ipbaclear(instr_retired_any[i],
	                                  baclears_any[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(baclears_any,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ipbaclear(a[i],
	                         b[i]);
	}
#endif
}


void
gms::system
::hsw_ipc_core_data(  const int64_t * __restrict instr_retired_any,
                      const int64_t * __restrict core_clks,
		      double * __restrict samples,
		      const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(core_clks,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ipc_core(instr_retired_any[i],
	                                  core_clks[i]);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(core_clks,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ipc_core(a[i],
	                         b[i]);
	}
#endif
}


void
gms::system
::hsw_ilp_data( const int64_t * __restrict uops_executed_core,
                const int64_t * __restrict execute_cycles,
		double * __restrict samples,
		const int32_t len,
		const bool is_ht_enabled) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(uops_executed_core,64);
	__assume_aligned(execute_cycles,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ilp(uops_execute_core[i],
	                              execute_cycles[i],
				      is_ht_enabled);
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(uops_execute_core,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(execute_cycles,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ilp(a[i],
	                   b[i],
			   is_ht_enabled);
	}
#endif
}


void
gms::system
::hsw_ip_mispredict_data(const int64_t * __restrict instr_retired_any,
                         const int64_t * __restrict br_misp_retired_all_branches,
		         double * __restrict samples,
		          const int32_t len){
		         
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(instr_retired_any,64);
	__assume_aligned(br_misp_retired_all_branches,64);
	__assume_aligned(samples,64);
#pragma simd vectorlengthfor(8)
        for(int32_t i = 0; i != len; ++i) {
            samples[i] = hsw_ip_mispredict(instr_retired_any[i],
	                         br_misp_retired_all_branches[i]);
				     
	}
#elif defined __GNUC__ && (!defined __ICC || !defined __INTEL_COMPILER)
        const int64_t * __restrict a = (const int64_t*)__builtin_assume_aligned(instr_retired_any,64);
	const int64_t * __restrict b = (const int64_t*)__builtin_assume_aligned(br_misp_retired_all_branches,64);
	double * __restrict        s = (double*)__builtin_assume_aligned(samples,64);
#pragma GCC ivdep
        for(int32_t i = 0; i != len; ++i) {
            s[i] = hsw_ip_mispredict(a[i],
	                   b[i])
			  
	}
#endif
}
