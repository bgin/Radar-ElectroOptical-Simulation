

#ifndef __GMS_POSTPROCESS_HW_EVENTS_METRICS_H__
#define __GMS_POSTPROCESS_HW_EVENTS_METRICS_H__

/*
     Statistical post-processing of HW events precomputed metrics.
     The functions arguments:
     const double * -- HW metric samples
     const char *   -- File name (the results will be written to the file)
     const int32_t  -- The length of HW metric samples
     Every function calls the TIMSAC "CANARM" subroutine
*/



// HW metric: CPU Operating Frequency

bool
postprocess_cpu_freq_canarm(const double * __restrict,
                            const char * __restrict,
                            const int32_t) __attribute__((hot))
				           __attribute__((aligned(32)));


// HW Metric: CPU Utilization (all cores)
bool
postprocess_cpu_util_canarm(const double * __restrict,
                            const char * __restrict,
			    const int32_t)  __attribute__((hot))
				            __attribute__((aligned(32)));


// HW Metric: CPU Utilization (all cores, kernel mode only)
bool
postprocess_cpu_util_km_canarm(const double * __restrict,
                               const char * __restrict,
			       const int32_t) __attribute__((hot))
			                      __attribute__((aligned(32)));


// HW Metric: Cycles Per Instruction (CPI)
bool
postprocess_cpi_canarm(const double * __restrict,
                       const char * __restrict,
		       const int32_t) __attribute__((hot))
			              __attribute__((aligned(32)));


// HW Metric: Cycles Per Instruction (CPI) (kernel mode only)
bool
postprocess_cpi_km_canarm(const double * __restrict,
                          const char * __restrict,
		          const int32_t) __attribute__((hot))
			                 __attribute__((aligned(32)));


// HW Metric: Branch Misprediction Ratio
bool
postprocess_br_misp_canarm(const double * __restrict,
                           const char * __restrict,
		           const int32_t) __attribute__((hot))
			                  __attribute__((aligned(32)));


// HW Metric: Loads Per Instruction
bool
postprocess_lpi_canarm(const double * __restrict,
                       const char * __restrict,
		       const int32_t) __attribute__((hot))
			              __attribute__((aligned(32)));


// HW Metric: Stores Per Instruction
bool
postprocess_spi_canarm(const double * __restrict,
                       const char * __restrict,
		       const int32_t) __attribute__((hot))
			              __attribute__((aligned(32)));


// HW Metric: Memory operations per instructions
bool
postprocess_mpi_canarm(const double * __restrict,
                       const char * __restrict,
		       const int32_t) __attribute__((hot))
			              __attribute__((aligned(32)));


// HW Metric: Locks retired per instruction
bool
postprocess_ret_locks_canarm(const double * __restrict,
                             const char * __restrict,
		             const int32_t) __attribute__((hot))
			                    __attribute__((aligned(32)));


// HW Metric: Uncacheable reads per instruction.
bool
postprocess_uncache_rpi_canarm(const double * __restrict,
                               const char * __restrict,
		               const int32_t) __attribute__((hot))
			                      __attribute__((aligned(32)));


// HW Metric: Streaming-stores (full line) per instruction.
bool
postprocess_stream_stores_fl_canarm(const double * __restrict,
                                    const char * __restrict,
		                    const int32_t) __attribute__((hot))
			                           __attribute__((aligned(32)));


// HW Metric:   Streaming-stores (partial line) per instruction.
bool
postprocess_stream_stores_pl_canarm(const double * __restrict,
                                    const char * __restrict,
		                    const int32_t) __attribute__((hot))
			                           __attribute__((aligned(32)));


// HW Metric:   L1D$ Misses per instruction
bool
postprocess_L1D_misses_canarm(const double * __restrict,
                              const char * __restrict,
		              const int32_t) __attribute__((hot))
			                     __attribute__((aligned(32)));


// HW Metric:  L1D$ Demand data read hit per instruction.
bool
postprocess_L1D_read_hits_canarm(  const double * __restrict,
                                   const char * __restrict,
		                    const int32_t) __attribute__((hot))
			                           __attribute__((aligned(32)));


// HW Metric:   L1I$ code read misses (including prefetches) per instruction.
bool
postprocess_L1I_read_misses_canarm( const double * __restrict,
                                    const char * __restrict,
		                    const int32_t) __attribute__((hot))
			                           __attribute__((aligned(32)));


// HW Metric:  L2 Demand data read hits per instruction.
bool
postprocess_L2_data_read_hits_canarm( const double * __restrict,
                                      const char * __restrict,
		                      const int32_t) __attribute__((hot))
			                           __attribute__((aligned(32)));


// HW Metric:   L2 Misses per instruction including
bool
postprocess_L2_all_misses_canarm(const double * __restrict,
                             const char * __restrict,
		             const int32_t) __attribute__((hot))
			                    __attribute__((aligned(32)));


// HW Metric:  L2 demand data read misses per instruction
bool
postprocess_L2_data_read_misses_canarm(const double * __restrict,
                                       const char * __restrict,
		                       const int32_t) __attribute__((hot))
			                              __attribute__((aligned(32)));


// HW Metric:  L2 demand data code misses per instruction
bool
postprocess_L2_code_misses_canarm(const double * __restrict,
                                  const char * __restrict,
		                  const int32_t) __attribute__((hot))
			                         __attribute__((aligned(32)));


// HW Metric:  L2 Any local request that HITM in a sibling core (per instruction).
bool
postprocess_L2_req_hitm_score_canarm(const double * __restrict,
                                     const char * __restrict,
		                     const int32_t) __attribute__((hot))
			                         __attribute__((aligned(32)));


// HW Metric:  L2 (percentage) of all lines evicted that are unused prefetches.
bool
postprocess_L2_evict_lines_up_canarm(const double * __restrict,
                                     const char * __restrict,
		                     const int32_t) __attribute__((hot))
			                         __attribute__((aligned(32)));


// HW Metric:  L2 percentage of L2 evictions that are allocated into L3$.
bool
postprocess_L2_evict_L3_alloc_canarm(const double * __restrict,
                                     const char * __restrict,
		                     const int32_t) __attribute__((hot))
			                         __attribute__((aligned(32)));


// HW Metric:  L2 percentage of L2 evictions that are NOT allocated into L3$.
bool
postprocess_L2_evict_L3_notalloc_canarm(const double * __restrict,
                                        const char * __restrict,
		                        const int32_t) __attribute__((hot))
			                         __attribute__((aligned(32)));


// HW Metric:  LLC code references per instruction (L3 prefetch excluded).
bool
postprocess_LLC_code_ref_canarm(const double * __restrict,
                                const char * __restrict,
		                const int32_t) __attribute__((hot))
			                       __attribute__((aligned(32)));


// HW Metric:  LLC data read references per instruction (L3 prefetch excluded).
bool
postprocess_LLC_data_ref_canarm(const double * __restrict,
                                const char * __restrict,
		                const int32_t) __attribute__((hot))
			                       __attribute__((aligned(32)));


// HW Metric:  LLC RFO references per instructions (L3 prefetch excluded).
bool
postprocess_LLC_rfo_ref_canarm(const double * __restrict,
                                const char * __restrict,
		                const int32_t) __attribute__((hot))
			                       __attribute__((aligned(32)));


// HW Metric:  LLC Misses Per Instruction (includes code and data and rfo with prefetches)
bool
postprocess_LLC_all_misses_canarm(  const double * __restrict,
                                    const char * __restrict,
		                    const int32_t) __attribute__((hot))
			                       __attribute__((aligned(32)));


// HW Metric:  LLC data read Misses Per Instruction (demand and prefetch).
bool
postprocess_LLC_data_read_misses_canarm(  const double * __restrict,
                                    const char * __restrict,
		                    const int32_t) __attribute__((hot))
			                       __attribute__((aligned(32)));

// HW Metric:  LLC RFO read Misses Per Instruction (demand and prefetch).  
bool
postprocess_LLC_rfo_read_misses_canarm(const double * __restrict,
                                    const char * __restrict,
		                    const int32_t) __attribute__((hot))
			                       __attribute__((aligned(32)));


// HW Metric:   LLC code read Misses Per Instruction (demand and prefetch).
bool
postprocess_LLC_code_read_misses_canarm(const double * __restrict,
                                        const char * __restrict,
		                        const int32_t) __attribute__((hot))
			                               __attribute__((aligned(32)));


// HW Metrics: LLC total 'HITM' per instruction
bool
postprocessing_LLC_total_hitm_canarm(   const double * __restrict,
                                        const char * __restrict,
		                        const int32_t) __attribute__((hot))
			                               __attribute__((aligned(32)));


// HW Metric:   LLC total HIT clean line forwards (per instr) (excludes LLC prefetches)
bool
postprocess_LLC_tot_hit_clean_lines_canarm(const double * __restrict,
                                           const char * __restrict,
		                           const int32_t) __attribute__((hot))
			                               __attribute__((aligned(32)));


// HW Metric:  Average LLC data read (demand+prefetch) miss latency (in ns).
bool
postprocess_L3_data_rd_miss_ns_canarm(const double * __restrict,
                                      const char * __restrict,
		                      const int32_t) __attribute__((hot))
			                             __attribute__((aligned(32)));


// HW Metric:    Average LLC data read (demand and prefetch) miss latency (in UNCORE clk).
bool
postprocess_L3_data_rd_miss_clk_canarm(const double * __restrict,
                                       const char * __restrict,
		                       const int32_t) __attribute__((hot))
			                              __attribute__((aligned(32)));


// HW Metric:  Average LLC data read (demand+prefetch) miss latency for LOCAL requests (in ns)
bool
postprocess_L3_data_rd_loc_miss_ns_canarm(const double * __restrict,
                                          const char * __restrict,
		                          const int32_t) __attribute__((hot))
			                                 __attribute__((aligned(32)));


// HW Metric: Average LLC data read (demand+prefetch) miss latency for REMOTE requests (in ns).
bool
postprocess_L3_data_rd_rem_miss_ns_canarm(const double * __restrict,
                                          const char * __restrict,
		                          const int32_t) __attribute__((hot))
			                                 __attribute__((aligned(32)));


// HW Metric:   Average LLC data read (demand+prefetch) miss latency  for REMOTE requests (in UNCORE clk)
bool
postprocess_L3_data_rd_rem_miss_clk_canarm(const double * __restrict,
                                          const char * __restrict,
		                          const int32_t) __attribute__((hot))
			                                 __attribute__((aligned(32)));


//HW Metric: ITLB MPI
bool
postprocess_ITLB_mpi_canarm(const double * __restrict,
                            const char * __restrict,
		            const int32_t) __attribute__((hot))
			                   __attribute__((aligned(32)));


// HW Metric:  ITLB large page MPI
bool
postprocess_ITLB_2M4M_mpi_canarm(const double * __restrict,
                                 const char * __restrict,
		                 const int32_t) __attribute__((hot))
			                   __attribute__((aligned(32)));


#endif /*__GMS_POSTPROCESS_HW_EVENTS_METRICS_H__*/
