

#ifndef __GMS_PREPROCESS_SKX_HW_METRICS_H__
#define __GMS_PREPROCESS_SKX_HW_METRICS_H__


#include <cstdint>


/*
    CPU operating frequency (in GHz)
*/
void
skx_cpu_operating_freq_samples(const double * __restrict,
			    const double * __restrict,
			    const double * __restrict,
			    double * __restrict,
			    const int32_t) __attribute__((hot))
                                           __attribute__((aligned(32)));

					   
/*
     CPU utilization (percentage) all cores.
*/
void
skx_cpu_utilization_samples(const double * __restrict,
                         const double * __restrict,
			 double * __restrict,
			 const int32_t)  __attribute__((hot))
                                         __attribute__((aligned(32)));


/*
    CPU utilization (percentage) in kernel mode (all cores).
*/
void
skx_cpu_utlization_kernel_samples(const double * __restrict,
                               const double * __restrict,
			       double * __restrict,
			       const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));


/*
    Cycles Per Instruction (CPI).
*/
void
skx_cycles_per_instr_samples(const double * __restrict,
                          const double * __restrict,
			  double * __restrict,
			  const int32_t)  __attribute__((hot))
                                          __attribute__((aligned(32)));


/*
    Cycles Per Instruction (CPI) kernel mode.
*/
void
skx_cycles_per_instr_kernel_samples(const double * __restrict,
                                 const double * __restrict,
			         double * __restrict,
			         const int32_t)  __attribute__((hot))
                                                 __attribute__((aligned(32)));


/*
    EMON event multiplexing reliability (>95% --- means satisfying ratio).
*/
void
skx_mux_reliability_samples(const double * __restrict,
                         const double * __restrict,
			 double * __restrict,
			 const int32_t)  __attribute__((hot))
                                         __attribute__((aligned(32)));


/*
      Branch mispredict ratio.
*/
void
skx_branch_mispred_ratio_samples(const double * __restrict,
                              const double * __restrict,
			      double * __restrict,
			      const int32_t)  __attribute__((hot))
                                              __attribute__((aligned(32)));


/*
    Loads per instruction.
*/
void
skx_loads_per_instr_samples(const double * __restrict,
                         const double * __restrict,
			 double * __restrict,
			 const int32_t)  __attribute__((hot))
                                         __attribute__((aligned(32)));


/*
    Stores per instruction.
*/
void
skx_stores_per_instr_samples(const double * __restrict,
                         const double * __restrict,
			 double * __restrict,
			 const int32_t)  __attribute__((hot))
                                         __attribute__((aligned(32)));


/*
   Memory operations per instruction.
*/
void
skx_mem_ops_instr_samples(const double * __restrict,
                       const double * __restrict,
		       const double * __restrict,
		       double * __restrict,
		       const int32_t) __attribute__((hot))
                                      __attribute__((aligned(32)));


/*
    Locks retired per instruction.
*/
void
skx_lock_per_instr_samples(const double * __restrict,
                        const double * __restrict,
			double * __restrict,
			const int32_t)  __attribute__((hot))
                                        __attribute__((aligned(32)));


/*
   Uncacheable reads per instruction.
*/
void
skx_uncacheable_reads_instr_samples(const double * __restrict,
                                 const double * __restrict,
			         double * __restrict,
			         const int32_t)  __attribute__((hot))
                                                 __attribute__((aligned(32)));


/*
    Streaming-stores (full line) per instruction.
*/
void
skx_stream_stores_fl_instr_samples(const double * __restrict,
                                const double * __restrict,
			        double * __restrict,
			        const int32_t)  __attribute__((hot))
                                                __attribute__((aligned(32)));


/*
     Streaming-stores (partial line) per instruction.
*/
void
skx_stream_stores_pl_instr_samples(const double * __restrict,
                                const double * __restrict,
			        double * __restrict,
			        const int32_t)  __attribute__((hot))
                                                __attribute__((aligned(32)));


/*
     L1D$ Misses per instruction including:
     -- data
     -- reads for ownership (rfo) with prefetches
*/
void
skx_L1D_misses_instr_samples(const double * __restrict,
                          const double * __restrict,
			  double * __restrict,
			  const int32_t)  __attribute__((hot))
                                          __attribute__((aligned(32)));


/*
     L1D$ Demand data read hit per instruction.
*/
void
skx_L1D_hits_instr_samples(const double * __restrict,
                        const double * __restrict,
			double * __restrict,
			const int32_t)  __attribute__((hot))
                                        __attribute__((aligned(32)));


/*
    L1$ code read misses (including prefetches) per instruction.
*/
void
skx_L1I_read_misses_instr_samples(const double * __restrict,
                               const double * __restrict,
			       double * __restrict,
			       const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));


/*
     L2D Demand data read hits per instruction.
*/
void
skx_L2_data_read_misses_instr_samples(const double * __restrict,
                                   const double * __restrict,
			           double * __restrict,
			           const int32_t)  __attribute__((hot))
                                                   __attribute__((aligned(32)));


/*
      L2 Misses per instruction including:
     -- code
     -- data
     -- read for ownership and prefetches
*/
void
skx_L2_all_misses_instr_samples(const double * __restrict,
                             const double * __restrict,
			     double * __restrict,
			     const int32_t)  __attribute__((hot))
                                             __attribute__((aligned(32)));


/*
       L2 demand data read misses per instruction.
*/
void
skx_L2_demand_data_read_mpi_samples(const double * __restrict,
                                 const double * __restrict,
			         double * __restrict,
			         const int32_t)  __attribute__((hot))
                                                 __attribute__((aligned(32)));


/*
     L2 demand data code misses per instruction.
*/
void
skx_L2_demand_code_mpi_samples(const double * __restrict,
                            const double * __restrict,
			    double * __restrict,
			    const int32_t)  __attribute__((hot))
                                            __attribute__((aligned(32)));


/*
    L2 Any local request that HITM in a sibling core (per instruction).
*/
void
skx_L2_req_hitm_sibling_core_samples(const double * __restrict,
                                  const double * __restrict,
			          double * __restrict,
			          const int32_t)  __attribute__((hot))
                                                  __attribute__((aligned(32)));


/*
     L2 (percentage) of all lines evicted that are unused prefetches.
*/
void
skx_L2_lines_evict_unused_prefetch_samples(const double * __restrict,
                                        const double * __restrict,
					const double * __restrict,
					const double * __restrict,
					double * __restrict,
					const int32_t)  __attribute__((hot))
                                                        __attribute__((aligned(32)));


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
void
skx_L2_evict_L3_alloc_samples(const double * __restrict,
                           const double * __restrict,
			   double * __restrict,
			   const int32_t)  __attribute__((hot))
                                           __attribute__((aligned(32)));


/*
      L2 percentage of L2 evictions that are NOT allocated into L3$.
*/
void
skx_L2_evict_L3_no_alloc_data(const double * __restrict,
                              const double * __restrict,
			      double * __restrict,
			      const int32_t)  __attribute__((hot))
                                              __attribute__((aligned(32)));


/*
     LLC code references per instruction (L3 prefetch excluded).
*/
void
skx_L3_code_references_instr_data(const double * __restrict,
                                  const double * __restrict,
			          double * __restrict,
			          const int32_t)  __attribute__((hot))
                                                  __attribute__((aligned(32)));


/*
    LLC data read references per instruction (L3 prefetch excluded).
*/
void
skx_L3_data_references_instr_data(const double * __restrict,
                                  const double * __restrict,
			          double * __restrict,
			          const int32_t)  __attribute__((hot))
                                                  __attribute__((aligned(32)));


/*
      LLC RFO references per instrctions (L3 prefetch excluded).
*/
void
skx_L3_rfo_references_instr_data(const double * __restrict,
                                  const double * __restrict,
			          double * __restrict,
			          const int32_t)  __attribute__((hot))
                                                  __attribute__((aligned(32)));


/*
     LLC Misses Per Instruction (includes code and data and rfo with prefetches).
*/
void
skx_L3_all_mpi(const double * __restrict,
               const double * __restrict,
	       const double * __restrict,
               const double * __restrict,
	       double * __restrict,
	       const int32_t)  __attribute__((hot))
                               __attribute__((aligned(32)));


/*
     LLC data read Misses Per Instruction (demand and prefetch).
*/
void
skx_L3_data_read_mpi_data(const double * __restrict,
                          const double * __restrict,
			  double * __restrict,
			  const int32_t)  __attribute__((hot))
                                          __attribute__((aligned(32)));


/*
      LLC RFO read Misses Per Instruction (demand and prefetch).  
*/
void
skx_L3_rfo_read_mpi_data(const double * __restrict,
                          const double * __restrict,
			  double * __restrict,
			  const int32_t)  __attribute__((hot))
                                          __attribute__((aligned(32)));


/*
     LLC code read MPI (demand and prefetch).
*/
void
skx_L3_total_hitm_instr_data(const double * __restrict,
                             const double * __restrict,
			     double * __restrict,
			     const int32_t)  __attribute__((hot))
                                             __attribute__((aligned(32)));


/*
     LLC total HIT clean line forwards (per instr) (excludes LLC prefetches)
*/
void
skx_L3_total_hitm_clean_lines_instr_data(const double * __restrict,
                                         const double * __restrict,
			                 double * __restrict,
			                 const int32_t)  __attribute__((hot))
                                                         __attribute__((aligned(32)));


/*
    Average LLC data read (demand+prefetch) miss latency (in ns).  
*/
void
skx_L3_avg_data_read_ns_data(const double * __restrict,
                             const double * __restrict,
			     const double * __restrict,
                             const double * __restrict,
			     double * __restrict,
			     const int32_t)  __attribute__((hot))
                                             __attribute__((aligned(32)));


/*
      Average LLC data read (demand and prefetch) miss latency (in UNCORE clk).
*/
void
skx_L3_avg_data_read_unc_clk_data(const double * __restrict,
                                  const double * __restrict,
			          double * __restrict,
			          const int32_t)  __attribute__((hot))
                                                  __attribute__((aligned(32)));


/*
      Average LLC data read (demand+prefetch) miss latency for LOCAL requests (in ns)
*/
void
skx_L3_avg_data_read_loc_ns_data(const double * __restrict,
                                 const double * __restrict,
			         const double * __restrict,
                                 const double * __restrict,
			         double * __restrict,
			         const int32_t)  __attribute__((hot))
                                             __attribute__((aligned(32)));


/*
    Average LLC data read (demand+prefetch) miss latency for REMOTE requests (in ns).
*/
void
skx_L3_avg_data_read_rem_ns_data(const double * __restrict,
                                 const double * __restrict,
			         const double * __restrict,
                                 const double * __restrict,
			         double * __restrict,
			         const int32_t)  __attribute__((hot))
                                             __attribute__((aligned(32)));


/*
      Average LLC data read (demand+prefetch) miss latency  for LOCAL requests (in UNCORE clk)
*/
void
skx_L3_avg_data_read_loc_unc_clk_data(const double * __restrict,
                                  const double * __restrict,
			          double * __restrict,
			          const int32_t)  __attribute__((hot))
                                                  __attribute__((aligned(32)));


/*
     Average LLC data read (demand+prefetch) miss latency  for REMOTE requests (in UNCORE clk)
*/
void
skx_L3_avg_data_read_rem_unc_clk(const double * __restrict,
                                  const double * __restrict,
			          double * __restrict,
			          const int32_t)  __attribute__((hot))
                                                  __attribute__((aligned(32)));


/*
    ITLB MPI
*/
void
skx_ITLB_mpi_data(const double * __restrict,
                  const double * __restrict,
		  double * __restrict,
		  const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
   ITLB large page MPI
*/
void
skx_ITLB_2M4M_mpi_data(const double * __restrict,
                  const double * __restrict,
		  double * __restrict,
		  const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
     DTLB load MPI.
*/
void
skx_DTLB_load_mpi_data(const double * __restrict,
                  const double * __restrict,
		  double * __restrict,
		  const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
    DTLB large page load MPI.
*/
void
skx_DTLB_2M4M_mpi_data(const double * __restrict,
                  const double * __restrict,
		  double * __restrict,
		  const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
   DTLB store MPI.
*/
void
skx_DTLB_store_mpi_data(const double * __restrict,
                  const double * __restrict,
		  double * __restrict,
		  const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
    DTLB load miss latency (in core clks)
*/
void
skx_DTLB_load_miss_clks_data(const double * __restrict,
                             const double * __restrict,
		             double * __restrict,
		             const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
     DTLB store miss latency (in core clks).
*/
void
skx_DTLB_store_miss_clks_data(const double * __restrict,
                             const double * __restrict,
		             double * __restrict,
		             const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
      ITLB miss latency (in core clks).
*/
void
skx_ITLB_miss_latency_clks_data(const double * __restrict,
                             const double * __restrict,
		             double * __restrict,
		             const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
    NUMA percentage of Reads addressed to local DRAM.
*/
void
skx_numa_reads_local_dram_samples(const double * __restrict,
                             const double * __restrict,
		             double * __restrict,
		             const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
    NUMA percentage of Reads addressed to remote  DRAM.
*/
void
skx_numa_reads_remote_dram_samples(const double * __restrict,
                             const double * __restrict,
		             double * __restrict,
		             const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
     Uncore Frequency Ghz.
*/
void
skx_uncore_frequency_samples(const double * __restrict,
                          const double * __restrict,
			  const double * __restrict,
		          double * __restrict,
		          const int32_t)  __attribute__((hot))
                                  __attribute__((aligned(32)));


/*
    UPI speed - GT/s.
*/
void
skx_UPI_speed_samples(const double * __restrict,
              const double * __restrict,
	      double * __restrict,
	      const int32_t)  __attribute__((hot))
                              __attribute__((aligned(32)));


/*
    UPI Data transmit BW (MB/sec) (only data)
*/
void
skx_UPI_data_bw_samples(const double * __restrict,
                        const double * __restrict,
	                double * __restrict,
	                const int32_t)  __attribute__((hot))
                              __attribute__((aligned(32)));


/*
     UPI Total transmit BW (MB/sec) (includes control)
*/
void
skx_UPI_tx_total_bw_samples(const double * __restrict,
                            const double * __restrict,
			    const double * __restrict,
	                    double * __restrict,
	                    const int32_t)  __attribute__((hot))
                              __attribute__((aligned(32)));


/*
     UPI Transmit utilization percentage (includes control).
   Percentage of time the processor is communicating with sibling processors.
*/
void
skx_UPI_tx_utlization_samples(const double * __restrict,
                              const double * __restrict,
			      const double * __restrict,
			      const double * __restrict,
			      const double * __restrict,
			      const double * __restrict,
	                      double * __restrict,
	                      const int32_t)  __attribute__((hot))
                              __attribute__((aligned(32)));


/*
   UPI percentage of  cycles transmit link is half-width (L0p) 
*/
void
skx_UPI_half_width_link_tx_cycles_samples( const double * __restrict,
			                   const double * __restrict,
			                   const double * __restrict,
	                                   double * __restrict,
	                                   const int32_t)  __attribute__((hot))
                                                           __attribute__((aligned(32)));


/*
    UPI percentage of  cycles receive link is half-width (L0p)
*/
void
skx_UPI_half_width_link_rx_cycles_samples( const double * __restrict,
			                   const double * __restrict,
			                   const double * __restrict,
	                                   double * __restrict,
	                                   const int32_t)  __attribute__((hot))
                                                           __attribute__((aligned(32)));


/*
    HA - Reads vs. all requests
*/
void
skx_HA_reads_vs_all_reqs_samples(const double * __restrict,
                                 const double * __restrict,
	                         double * __restrict,
	                         const int32_t)  __attribute__((hot))
                                                 __attribute__((aligned(32)));


/*
    HA - Writes vs. all requests
*/
void
skx_HA_writes_vs_all_reqs_samples(const double * __restrict,
                                 const double * __restrict,
	                         double * __restrict,
	                         const int32_t)  __attribute__((hot))
                                                 __attribute__((aligned(32)));


/*
      HA percentage of all reads that are local.
*/
void
skx_HA_all_reads_local_samples(const double * __restrict,
                               const double * __restrict,
	                       double * __restrict,
	                       const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));


/*
     HA percentage of all writes that are local.
*/
void
skx_HA_all_writes_local_samples(const double * __restrict,
                               const double * __restrict,
	                       double * __restrict,
	                       const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));


/*
     HA conflict responses per instruction.
*/
void
skx_HA_conflict_resp_samples(const double * __restrict,
                               const double * __restrict,
	                       double * __restrict,
	                       const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));


/*
    HA directory lookups that spawned a snoop (per instruction)
*/
void
skx_HA_dir_lookup_snoop_samples(const double * __restrict,
                               const double * __restrict,
	                       double * __restrict,
	                       const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));


/*
    HA directory lookups that did not spawn a snoop (per instruction).
*/
void
skx_HA_dir_lookup_no_snoop_samples(const double * __restrict,
                               const double * __restrict,
	                       double * __restrict,
	                       const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));


/*
    M2M directory updates (per instruction). 
*/
void
skx_M2M_dir_updates_samples(const double * __restrict,
                            const double * __restrict,
			    const double * __restrict,
			    const double * __restrict,
	                    double * __restrict,
	                    const int32_t)  __attribute__((hot))
                                               __attribute__((aligned(32)));

#endif /*__GMS_PREPROCESS_SKX_HW_METRICS_H__*/
