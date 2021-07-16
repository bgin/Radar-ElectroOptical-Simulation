

#ifndef __SKX_HW_EVENTS_METRICS_HPP__
#define __SKX_HW_EVENTS_METRICS_HPP__


#include <cstdint>


/*
   CPU operating frequency (in GHz)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_cpu_operating_freq(const uint64_t CPU_CLK_UNHALTED_THREAD,
                              const uint64_t CPU_CLK_UNHALTED_REF_TSC,
                              const uint64_t TSC_FREQUENCY) {

        return ((double)(CPU_CLK_UNHALTED_THREAD/
                           CPU_CLK_UNHALTED_REF_TSC*TSC_FREQUENCY)/
                                                      1000000000ULL);
}

/*
   CPU utilization (percentage) all cores.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_cpu_utilization(const uint64_t CPU_CLK_UNHALTED_REF_TSC,
                           const uint64_t TIME_STAMP_CYCLES) {

        return ((double)(100ULL*CPU_CLK_UNHALTED_REF_TSC/
                                            TIME_STAMP_CYCLES));
}


/*
    CPU utilization (percentage) in kernel mode (all cores).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_cpu_utilization_kernel(const uint64_t CPU_CLK_UNHALTED_REF_TSC_SUP,
                                  const uint64_t TIME_STAMP_CYCLES) {

        return ((double)(100ULL*CPU_CLK_UNHALTED_REF_TSC_SUP/
                                                  TIME_STAMP_CYCLES));
}


/*
    Cycles Per Instruction (CPI).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_cycles_per_instr(const uint64_t CPU_CLK_UNHALTED_THREAD,
                            const uint64_t INST_RETIRED_ANY) {

        return ((double)(CPU_CLK_UNHALTED_THREAD/
                                    INST_RETIRED_ANY));
}


/*
     Cycles Per Instruction (CPI) kernel mode.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_cycles_per_instr_kernel(const uint64_t CPU_CLK_UNHALTED_THREAD_SUP,
                                   const uint64_t INST_RETIRED_ANY_SUP) {

        return ((double)(CPU_CLK_UNHALTED_THREAD_SUP/
                                      INST_RETIRED_ANY_SUP));
}


/*
    EMON event multiplexing reliability (>95% --- means satisfying ratio).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_emon_mux_reliability(const uint64_t CPU_CLK_UNHALTED_THREAD_P,
                                const uint64_t CPU_CLK_UNHALTED_THREAD) {

        return ((CPU_CLK_UNHALTED_THREAD_P-
                            CPU_CLK_UNHALTED_THREAD < 0ULL) ? 
                                               (double)(CPU_CLK_UNHALTED_THREAD_P/
                                                             CPU_CLK_UNHALTED_THREAD) : 
                                                                   (double)(CPU_CLK_UNHALTED_THREAD/
                                                                                CPU_CLK_UNHALTED_THREAD_P));
}


/*
    Branch mispredict ratio.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_branch_mispredict_ratio(const uint64_t BR_MISP_RETIRED_ALL_BRANCHES,
                                   const uint64_t BR_INST_RETIRED_ALL_BRANCHES) {

        return ((double)(BR_MISP_RETIRED_ALL_BRANCHES/
                                  BR_INST_RETIRED_ALL_BRANCHES));
}


/*
    Loads per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_loads_per_instr(const uint64_t MEM_INST_RETIRED_ALL_LOADS,
                           const uint64_t INST_RETIRED_ANY) {

        return ((double)(MEM_INST_RETIRED_ALL_LOADS/
                                  INST_RETIRED_ANY));
}


/*
    Stores per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_stores_per_instr(const uint64_t MEM_INST_RETIRED_ALL_STORES,
                            const uint64_t INST_RETIRED_ANY) {

        return ((double)(MEM_INST_RETIRED_ALL_STORES/
                                   INST_RETIRED_ANY));
}


/*
    Memory operations per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_mem_ops_per_instr(const uint64_t MEM_INST_RETIRED_ALL_LOADS,
                             const uint64_t MEM_INST_RETIRED_ALL_STORES,
                             const uint64_t INST_RETIRED_ANY) {

        const uint64_t mem_ops = MEM_INST_RETIRED_ALL_LOADS+
                                 MEM_INST_RETIRED_ALL_STORES;
        return ((double)(mem_ops/INST_RETIRED_ANY));
}


/*
    Locks retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_locks_per_instr(const uint64_t MEM_INST_RETIRED_LOCK_LOADS,
                           const uint64_t INST_RETIRED_ANY) {

        return ((double)(MEM_INST_RETIRED_LOCK_LOADS/
                                       INST_RETIRED_ANY));
}


/*
    Uncacheable reads per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_uncacheable_reads_per_instr(const uint64_t UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40e33,
                                       const uint64_t INST_RETIRED_ANY) {
        
        return ((double)(UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40e33/
                                                     INST_RETIRED_ANY));
}


/*
    Streaming-stores (full line) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_streaming_stores_fl_per_instr(const uint64_t UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41833,
                                      const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41833/
                                                      INST_RETIRED_ANY));
}


/*
     Streaming-stores (partial line) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_streaming_stores_pl_per_instr(const uint64_t UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41a33,
                                         const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41a33/
                                                        INST_RETIRED_ANY));
}


/*
     L1D$ Misses per instruction including:
     -- data
     -- reads for ownership (rfo) with prefetches
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L1D_misses_per_instr(const uint64_t L1D_REPLACEMENT,
                                const uint64_t INST_RETIRED_ANY) {

        return ((double)(L1D_REPLACEMENT/INST_RETIRED_ANY));
}


/*
    L1D$ Demand data read hit per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L1D_hits_per_instr(const uint64_t MEM_LOAD_RETIRED_L1_HIT,
                              const uint64_t INST_RETIRED_ANY) {

        return ((double)(MEM_LOAD_RETIRED_L1_HIT/INST_RETIRED_ANY));
}


/*
     L1$ code read misses (including prefetches) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L1I_read_misses_per_instr(const uint64_t L2_RQSTS_ALL_CODE_RD,
                                     const uint64_t INST_RETIRED_ANY) {

        return ((double)(L2_RQSTS_ALL_CODE_RD/INST_RETIRED_ANY));
}


/*
     L2D Demand data read hits per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L2_data_read_hits_per_instr(const uint64_t MEM_LOAD_RETIRED_L2_HIT,
                                       const uint64_t INST_RETIRED_ANY) {

        return ((double)(MEM_LOAD_RETIRED_L2_HIT/INST_RETIRED_ANY));
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
static inline
double skx_L2_all_misses_per_instr(const uint64_t L2_LINES_IN_ALL,
                                   const uint64_t INST_RETIRED_ANY) {

        return ((double)(L2_LINES_IN_ALL/INST_RETIRED_ANY));
}


/*
    L2 demand data read misses per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L2_demand_data_read_mpi(const uint64_t MEM_LOAD_RETIRED_L2_MISS,
                                   const uint64_t INST_RETIRED_ANY) {

        return ((double)(MEM_LOAD_RETIRED_L2_MISS/INST_RETIRED_ANY));
}


/*
    L2 demand data code misses per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L2_demand_code_mpi(const uint64_t L2_RQSTS_CODE_RD_MISS,
                              const uint64_t INST_RETIRED_ANY) {

        return ((double)(L2_RQSTS_CODE_RD_MISS/INST_RETIRED_ANY));
}


/*
    L2 Any local request that HITM in a sibling core (per instruction).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L2_request_hitm_sibling_core(const uint64_t OFFCORE_RESPONSE_request_ALL_READS_response_L3_HIT_HITM_OTHER_CORE,
                                        const uint64_t INST_RETIRED_ANY) {

        return ((double)(OFFCORE_RESPONSE_request_ALL_READS_response_L3_HIT_HITM_OTHER_CORE/
                                                                                 INST_RETIRED_ANY));
}


/*
    L2 (percentage) of all lines evicted that are unused prefetches.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L2_lines_evict_unused_prefetch(const uint64_t L2_LINES_OUT_USELESS_HWPF,
                                          const uint64_t LINES_OUT_NON_SILENT,
                                          const uint64_t LINES_OUT_SILENT,
                                          const uint64_t HW_Thread_Count) {

        const double term1 = (double)(LINES_OUT_SILENT/HW_Thread_Count);
        return ((double)(100ULL*L2_LINES_OUT_USELESS_HWPF/
                                          (LINES_OUT_NON_SILENT+term1)));
}


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double  skx_L2_evict_L3_allocated(const uint64_t L2_LINES_OUT_NON_SILENT,
                                  const uint64_t IDI_MISC_WB_DOWNGRADE) {

        const double term1 = (double)(L2_LINES_OUT_NON_SILENT-
                                                 IDI_MISC_WB_DOWNGRADE);
        return ((double)(100ULL*term1/L2_LINES_OUT_NON_SILENT));
}


/*
    L2 percentage of L2 evictions that are NOT allocated into L3$.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L2_evict_L3_not_allocated(const uint64_t L2_LINES_OUT_NON_SILENT,
                                     const uint64_t IDI_MISC_WB_DOWNGRADE) {

        return ((double)(100ULL*IDI_MISC_WB_DOWNGRADE/
                                        L2_LINES_OUT_NON_SILENT));
}


/*
    LLC code references per instruction (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double  skx_L3_code_references_per_instr(const uint64_t UNC_CHA_TOR_INSERTS_IA_filter1_0x40433,
                                         const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_filter1_0x40433/
                                              INST_RETIRED_ANY));
}


/*
    LLC data read references per instruction (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L3_data_references_per_instr(const uint64_t UNC_CHA_TOR_INSERTS_IA_filter1_0x40433,
                                        const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_filter1_0x40433/
                                                  INST_RETIRED_ANY));
}


/*
    LLC RFO references per instrctions (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
static inline
double skx_L3_rf_references_per_instr(const uint64_t UNC_CHA_TOR_INSERTS_IA_filter1_0x40033,
                                      const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_filter1_0x40033/
                                                     INST_RETIRED_ANY));
}


/*
    LLC Misses Per Instruction (includes code and data and rfo with prefetches).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
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
static inline
double skx_L3_data_read_mpi(const uint64_t UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433,
                            const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433/
                                                      INST_RETIRED_ANY));
}


/*
    LLC RFO read Misses Per Instruction (demand and prefetch).  
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
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
static inline
double skx_L3_code_read_mpi(const uint64_t UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233,
                            const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233/
                                                           INST_RETIRED_ANY));
}



#endif /*__SKX_HW_EVENTS_METRICS_HPP__*/