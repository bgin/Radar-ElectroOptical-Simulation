

#ifndef __GMS_SKX_HW_EVENTS_METRICS_HPP__
#define __GMS_SKX_HW_EVENTS_METRICS_HPP__


#include <omp.h>


/*
   CPU operating frequency (in GHz)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cpu_operating_freq(const double CPU_CLK_UNHALTED_THREAD,
                              const double CPU_CLK_UNHALTED_REF_TSC,
                              const double TSC_FREQUENCY) {

  return ((CPU_CLK_UNHALTED_THREAD/
                           CPU_CLK_UNHALTED_REF_TSC*TSC_FREQUENCY)/
                                                      1000000000ULL);
}

/*
   CPU utilization (percentage) all cores.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cpu_utilization(const double CPU_CLK_UNHALTED_REF_TSC,
                           const double TIME_STAMP_CYCLES) {

        return (100.0*CPU_CLK_UNHALTED_REF_TSC/
                                      TIME_STAMP_CYCLES);
}


/*
    CPU utilization (percentage) in kernel mode (all cores).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cpu_utilization_kernel(const double CPU_CLK_UNHALTED_REF_TSC_SUP,
                                  const double TIME_STAMP_CYCLES) {

        return (100.0*CPU_CLK_UNHALTED_REF_TSC_SUP/
                                        TIME_STAMP_CYCLES);
}


/*
    Cycles Per Instruction (CPI).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cycles_per_instr(const double CPU_CLK_UNHALTED_THREAD,
                            const double INST_RETIRED_ANY) {

        return (CPU_CLK_UNHALTED_THREAD/
                                    INST_RETIRED_ANY);
}


/*
     Cycles Per Instruction (CPI) kernel mode.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cycles_per_instr_kernel(const double CPU_CLK_UNHALTED_THREAD_SUP,
                                   const double INST_RETIRED_ANY_SUP) {

        return (CPU_CLK_UNHALTED_THREAD_SUP/
                                      INST_RETIRED_ANY_SUP);
}


/*
    EMON event multiplexing reliability (>95% --- means satisfying ratio).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_emon_mux_reliability(const double CPU_CLK_UNHALTED_THREAD_P,
                                const double CPU_CLK_UNHALTED_THREAD) {

        return (CPU_CLK_UNHALTED_THREAD_P-
                            CPU_CLK_UNHALTED_THREAD < 0.0 ? 
                                               CPU_CLK_UNHALTED_THREAD_P/
                                                             CPU_CLK_UNHALTED_THREAD : 
                                                                   CPU_CLK_UNHALTED_THREAD/
                                                                                CPU_CLK_UNHALTED_THREAD_P);
}


/*
    Branch mispredict ratio.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_branch_mispredict_ratio(const double BR_MISP_RETIRED_ALL_BRANCHES,
                                   const double BR_INST_RETIRED_ALL_BRANCHES) {

        return (BR_MISP_RETIRED_ALL_BRANCHES/
                                  BR_INST_RETIRED_ALL_BRANCHES);
}


/*
    Loads per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_loads_per_instr(const double MEM_INST_RETIRED_ALL_LOADS,
                           const double INST_RETIRED_ANY) {

        return (MEM_INST_RETIRED_ALL_LOADS/
                                  INST_RETIRED_ANY);
}


/*
    Stores per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_stores_per_instr(const double MEM_INST_RETIRED_ALL_STORES,
                            const double INST_RETIRED_ANY) {

        return (MEM_INST_RETIRED_ALL_STORES/
                                   INST_RETIRED_ANY);
}


/*
    Memory operations per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_ops_per_instr(const double MEM_INST_RETIRED_ALL_LOADS,
                             const double MEM_INST_RETIRED_ALL_STORES,
                             const double INST_RETIRED_ANY) {

        const double mem_ops = MEM_INST_RETIRED_ALL_LOADS+
                                 MEM_INST_RETIRED_ALL_STORES;
        return (mem_ops/INST_RETIRED_ANY);
}


/*
    Locks retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_locks_per_instr(const double MEM_INST_RETIRED_LOCK_LOADS,
                           const double INST_RETIRED_ANY) {

        return (MEM_INST_RETIRED_LOCK_LOADS/
                                       INST_RETIRED_ANY);
}


/*
    Uncacheable reads per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_uncacheable_reads_per_instr(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40e33,
                                       const double INST_RETIRED_ANY) {
        
        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40e33/
                                                     INST_RETIRED_ANY);
}


/*
    Streaming-stores (full line) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_streaming_stores_fl_per_instr(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41833,
                                         const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41833/
                                                      INST_RETIRED_ANY);
}


/*
     Streaming-stores (partial line) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_streaming_stores_pl_per_instr(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41a33,
                                         const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41a33/
                                                        INST_RETIRED_ANY);
}


/*
     L1D$ Misses per instruction including:
     -- data
     -- reads for ownership (rfo) with prefetches
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L1D_misses_per_instr(const double L1D_REPLACEMENT,
                                const double INST_RETIRED_ANY) {

        return (L1D_REPLACEMENT/INST_RETIRED_ANY);
}


/*
    L1D$ Demand data read hit per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L1D_hits_per_instr(const double MEM_LOAD_RETIRED_L1_HIT,
                              const double INST_RETIRED_ANY) {

        return (MEM_LOAD_RETIRED_L1_HIT/INST_RETIRED_ANY);
}


/*
     L1$ code read misses (including prefetches) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L1I_read_misses_per_instr(const double L2_RQSTS_ALL_CODE_RD,
                                     const double INST_RETIRED_ANY) {

        return (L2_RQSTS_ALL_CODE_RD/INST_RETIRED_ANY);
}


/*
     L2D Demand data read hits per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_data_read_hits_per_instr(const double MEM_LOAD_RETIRED_L2_HIT,
                                       const double INST_RETIRED_ANY) {

        return (MEM_LOAD_RETIRED_L2_HIT/INST_RETIRED_ANY);
}


/*
     L2 Misses per instruction including:
     -- code
     -- data
     -- read for ownership and prefetches
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_all_misses_per_instr(const double L2_LINES_IN_ALL,
                                   const double INST_RETIRED_ANY) {

        return (L2_LINES_IN_ALL/INST_RETIRED_ANY);
}


/*
    L2 demand data read misses per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_demand_data_read_mpi(const double MEM_LOAD_RETIRED_L2_MISS,
                                   const double INST_RETIRED_ANY) {

        return (MEM_LOAD_RETIRED_L2_MISS/INST_RETIRED_ANY);
}


/*
    L2 demand data code misses per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_demand_code_mpi(const double L2_RQSTS_CODE_RD_MISS,
                              const double INST_RETIRED_ANY) {

        return (L2_RQSTS_CODE_RD_MISS/INST_RETIRED_ANY);
}


/*
    L2 Any local request that HITM in a sibling core (per instruction).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_request_hitm_sibling_core(const double OFFCORE_RESPONSE_request_ALL_READS_response_L3_HIT_HITM_OTHER_CORE,
                                        const double INST_RETIRED_ANY) {

        return (OFFCORE_RESPONSE_request_ALL_READS_response_L3_HIT_HITM_OTHER_CORE/
                                                                    INST_RETIRED_ANY);
}


/*
    L2 (percentage) of all lines evicted that are unused prefetches.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_lines_evict_unused_prefetch(const double L2_LINES_OUT_USELESS_HWPF,
                                          const double LINES_OUT_NON_SILENT,
                                          const double LINES_OUT_SILENT,
                                          const double HW_Thread_Count) {

        const double term1 = LINES_OUT_SILENT/HW_Thread_Count;
        return (100.0*L2_LINES_OUT_USELESS_HWPF/
                                          (LINES_OUT_NON_SILENT+term1));
}


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double  skx_L2_evict_L3_allocated(const double L2_LINES_OUT_NON_SILENT,
                                  const double IDI_MISC_WB_DOWNGRADE) {

        const double term1 = L2_LINES_OUT_NON_SILENT-
                                                 IDI_MISC_WB_DOWNGRADE;
        return (100.0*term1/L2_LINES_OUT_NON_SILENT);
}


/*
    L2 percentage of L2 evictions that are NOT allocated into L3$.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_evict_L3_not_allocated(const double L2_LINES_OUT_NON_SILENT,
                                     const double IDI_MISC_WB_DOWNGRADE) {

        return (100.0*IDI_MISC_WB_DOWNGRADE/
                                        L2_LINES_OUT_NON_SILENT);
}


/*
    LLC code references per instruction (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double  skx_L3_code_references_per_instr(const double UNC_CHA_TOR_INSERTS_IA_filter1_0x40433,
                                         const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_filter1_0x40433/
                                              INST_RETIRED_ANY);
}


/*
    LLC data read references per instruction (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_data_references_per_instr(const double UNC_CHA_TOR_INSERTS_IA_filter1_0x40433,
                                        const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_filter1_0x40433/
                                                  INST_RETIRED_ANY);
}


/*
    LLC RFO references per instrctions (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_rf_references_per_instr(const double UNC_CHA_TOR_INSERTS_IA_filter1_0x40033,
                                      const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_filter1_0x40033/
                                                     INST_RETIRED_ANY);
}


/*
    LLC Misses Per Instruction (includes code and data and rfo with prefetches).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_all_mpi(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433,
                   const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233,
                   const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033,
                   const double INST_RETIRED_ANY) {

        const double sum = UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433+
                           UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233+
                           UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033;
        return (sum/INST_RETIRED_ANY);
}


/*
    LLC data read Misses Per Instruction (demand and prefetch).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_data_read_mpi(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433,
                            const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433/
                                                      INST_RETIRED_ANY);
}


/*
    LLC RFO read Misses Per Instruction (demand and prefetch).  
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_rfo_read_mpi(const uint64_t UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033,
                           const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033/
                                                        INST_RETIRED_ANY));
}


/*
     LLC code read Misses Per Instruction (demand and prefetch).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_code_read_mpi(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233,
                            const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233/
                                                           INST_RETIRED_ANY);
}


/*
      LLC code read MPI (demand and prefetch).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3total_hitm_per_instr(const double OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HITM,
                                const double INST_RETIRED_ANY) {

       return (OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HITM/
                                                                   INST_RETIRED_ANY);
}


/*
        LLC total HIT clean line forwards (per instr) (excludes LLC prefetches)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_total_hitm_clean_lines_per_instr(const double OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HIT_FORWARD,
                                            const double INST_RETIRED_ANY) {

       return (OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HIT_FORWARD/
                                                                          INST_RETIRED_ANY);
}


/*
        Average LLC data read (demand+prefetch) miss latency (in ns). 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_ns(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433,
                               const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433,
			       const double UNC_CHA_CLOCKTICKS,
			       const double TIME_INTERVAL_SEC) {

        const double term1 = UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433/
	                     UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433;
	return (1000000000.0*TIME_INTERVAL_SEC*term1/UNC_CHA_CLOCKTICKS);
}


/*
      Average LLC data read (demand and prefetch) miss latency (in UNCORE clk).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_unc_clk(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433,
                                    const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433) {

       return (UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433/
               UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433);
}


/*
    Average LLC data read (demand+prefetch) miss latency for LOCAL requests (in ns)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_loc_req_ns(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432,
                                       const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432,
				       const double UNC_CHA_CLOCKTICKS,
			               const double TIME_INTERVAL_SEC) {
    
        const double term1 = UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432/
	                     UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432;
	return (1000000000.0*TIME_INTERVAL_SEC*term1/UNC_CHA_CLOCKTICKS);
}


/*
     Average LLC data read (demand+prefetch) miss latency for REMOTE requests (in ns).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_rem_req_ns(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40431,
                                       const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431,
				       const double UNC_CHA_CLOCKTICKS,
			               const double TIME_INTERVAL_SEC) {

        const double term1 = UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40431/
	                     UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431;
	return (1000000000.0*TIME_INTERVAL_SEC*term1/UNC_CHA_CLOCKTICKS);
}


/*
      Average LLC data read (demand+prefetch) miss latency  for LOCAL requests (in UNCORE clk)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_loc_req_unc_clk(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432,
                                            const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432) {

       return ( UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432/
                UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432);
}


/*
      Average LLC data read (demand+prefetch) miss latency  for REMOTE requests (in UNCORE clk)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_data_read_rem_req_unc_clk(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432,
                                        const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432) {

       return (UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432/
               UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432);
}


/*
   ITLB MPI. 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ITLB_mpi(const double ITLB_MISSES_WALK_COMPLETED,
                    const double INT_RETIRED_ANY) {

       return (ILTB_MISSES_WALK_COMPLETED/INST_RETIRED_ANY);
}


/*
    ITLB large page MPI
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ITLB_2M_4M_mpi(const double ITLB_MISSES_WALK_COMPLETED_2M_4M,
                          const double INST_RETIRED_ANY) {

       return (ITLB_MISSES_WALK_COMPLETED_2M_4M/INST_RETIRED_ANY);
}


/*
    DTLB load MPI.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_load_mpi(const double DTLB_LOAD_MISSES_WALK_COMPLETED,
                         const double INST_RETIRED_ANY) {

       return (DTLB_LOAD_MISSES_WALK_COMPLETED/INST_RETIRED_ANY);
}


/*
    DTLB large page load MPI.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_2M_4M_mpi(const double DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
                          const double INST_RETIRED_ANY) {

       return (DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M/
                                        INST_RETIRED_ANY);
}


/*
    DTLB store MPI.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_store_mpi(const double DTLB_STORE_MISSES_WALK_COMPLETED,
                          const double INST_RETIRED_ANY) {

       return (DTLB_STORE_MISSES_WALK_COMPLETED/
                                INST_RETIRED_ANY);
}


/*
    DTLB load miss latency (in core clks)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_load_miss_clks(const double DTLB_LOAD_MISSES_WALK_ACTIVE,
                               const double DTLB_LOAD_MISSES_WALK_COMPLETED) {

       return (DTLB_LOAD_MISSES_WALK_ACTIVE/
                DTLB_LOAD_MISSES_WALK_COMPLETED);
}


/*
    DTLB store miss latency (in core clks).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_store_miss_clks(const double DTLB_STORE_MISSES_WALK_ACTIVE,
                                const double DTLB_STORE_MISSES_WALK_COMPLETED) {

       return ( DTLB_STORE_MISSES_WALK_ACTIVE/
                DTLB_STORE_MISSES_WALK_COMPLETED);
}          


/*
    ITLB miss latency (in core clks).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ITLB_miss_latency_clks(const double ITLB_MISSES_WALK_ACTIVE,
                                  const double ITLB_MISSES_WALK_COMPLETED) {

       return (ITLB_MISSES_WALK_ACTIVE/
                 ITLB_MISSES_WALK_COMPLETED);
}


/*
     NUMA percentage of Reads addressed to local DRAM.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_numa_reads_local_dram(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432,
                                 const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431) {

       const double sum = UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432+
                          UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431;
       return (100.0*UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432/sum);
}


/*
     NUMA percentage of Reads addressed to remote  DRAM.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_numa_reads_remote_dram(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432,
                                  const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431){

       const double sum = UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432+
                          UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431;
       return (100.0*UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431/sum);
}


/*
      Uncore Frequency Ghz.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_uncore_frequency_ghz(const double UNC_CHA_CLOCKTICKS,
                                const double HW_CORE_COUNT,
				const double TIME_INTERVAL_SECONDS) {
                              
       const double term1 = UNC_CHA_CLOCKTICKS/HW_CORE_COUNT;
       const double term2 = 1000000000.0/TIME_INTERVAL_SECONDS;
       return (term1/term2);
}


/*
    UPI speed - GT/s.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_speed(const double UNC_UPI_CLOCKTICKS,
                     const double TIME_INTERVAL_SECONDS) {

       const double term1 = 8.0*UNC_UPI_CLOCKTICKS;
       return (term1/1000000000.0/TIME_INTERVAL_SECONDS);
}


/*
    UPI Data transmit BW (MB/sec) (only data)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_data_bw(const double UNC_UPI_TxL_FLITS_ALL_DATA,
                       const double TIME_INTERVAL_SECONDS) {

       const double term1 = UNC_UPI_TxL_FLITS_ALL_DATA*
                            7.1111111111111111111111;
       const double term2 = 1000000/TIME_INTERVAL_SECONDS;
       return (term1/term2);
}


/*
    UPI Total transmit BW (MB/sec) (includes control)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_transmit_total_bw(const double UNC_UPI_TxL_FLITS_ALL_DATA,
                                 const double UNC_UPI_TxL_FLITS_NON_DATA,
				 const double TIME_INTERVAL_SECONDS) {
       
       const double term1 =  UNC_UPI_TxL_FLITS_ALL_DATA+
                             UNC_UPI_TxL_FLITS_NON_DATA*
			     7.1111111111111111111111;
       const double term2 = 1000000/TIME_INTERVAL_SECONDS;
       return (term1/term2);
}


/*
   UPI Transmit utilization percentage (includes control).
   Percentage of time the processor is communicating with sibling processors.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_transmit_utilization(const double UNC_UPI_CLOCKTICKS,
                                    const double UNC_UPI_L1_POWER_CYCLES,
				    const double FREERUN_PKG_C6_RESIDENCY,
				    const double UNC_UPI_TxL_FLITS_ALL_DATA,
				    const double UNC_UPI_TxL_FLITS_NON_DATA,
				    const double TIME_INTERVAL_SECONDS) {

        const double term1 = UNC_UPI_TxL_FLITS_ALL_DATA+
                             UNC_UPI_TxL_FLITS_NON_DATA*
			     0.333333333333333333333333;
	const double term2 = TIME_INTERVAL_SECONDS/(TIME_INTERVAL_SECONDS-
	                                            FREERUN_PKG_C6_RESIDENCY);
	const double term3 = UNC_UPI_CLOCKTICKS-UNC_UPI_L1_POWER_CYCLES;
	return (100.0*term1/term2*term3*0.8333333333333333333333);
}


/*
   UPI percentage of  cycles transmit link is half-width (L0p) 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_half_width_link_tx_cycles(const double UNC_UPI_TxL0P_POWER_CYCLES,
                                     const double UNC_UPI_CLOCKTICKS,
				     const double UNC_UPI_L1_POWER_CYCLES) {

        const double term1 = UNC_UPI_TxL0P_POWER_CYCLES/
	                     (UNC_UPI_CLOCKTICKS-UNC_UPI_1_POWER_CYCLES);
	return (100.0*term1);
}


/*
    UPI percentage of  cycles receive link is half-width (L0p)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_half_width_link_rx_cycles(const double UNC_UPI_RxL0P_POWER_CYCLES,
                                     const double UNC_UPI_CLOCKTICKS,
				     const double UNC_UPI_L1_POWER_CYCLES) {

        const double term1 = UNC_UPI_RxL0P_POWER_CYCLES/
	                     (UNC_UPI_CLOCKTICKS-UNC_UPI_1_POWER_CYCLES);
	return (100.0*term1);
}


/*
     HA - Reads vs. all requests
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_reads_vs_all_req(const double UNC_CHA_REQUESTS_READS,
                               const double UNC_CHA_REQUESTS_WRITES) {

       return (UNC_CHA_REQUESTS_READS/
              (UNC_CHA_REQUESTS_READS+UNC_CHA_REQUESTS_WRITES));
}




#endif /*__GMS_SKX_HW_EVENTS_METRICS_HPP__*/
