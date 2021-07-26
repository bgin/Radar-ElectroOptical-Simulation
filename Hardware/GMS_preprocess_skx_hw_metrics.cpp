

#include <omp.h>
#include "GMS_skx_hw_events_metrics.hpp"



/*
   CPU operating frequency (in GHz)
*/
void
skx_cpu_operating_freq_samples(const double * __restrict __attribute__((aligned(64))) a,
			       const double * __restrict __attribute__((aligned(64))) b,
			       const double * __restrict __attribute__((aligned(64))) c,
			       double * __restrict __attribute__((aligned(64))) d,
			       const int32_t data_len) {

       
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,d,data_len) private(i) \
        aligned(a:64,b,c,d) linear(i:1) unroll partial(10) if(data_len>=10000)
           for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_cpu_operating_freq(a[i],
	                                     b[i],
					     c[i]);
	   }       
}


/*
   CPU utilization (percentage) all cores.
*/
void
skx_cpu_utilization_samples(const double * __restrict __attribute__((aligned(64))) a,
                            const double * __restrict __attribute__((aligned(64))) b,
			    double * __restrict __attribute__((aligned(64))) c,
			    const int32_t data_len) {

      
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cpu_utilization(a[i],
	                                  b[i]);
					
          }
}


/*
    CPU utilization (percentage) in kernel mode (all cores).
*/
void
skx_cpu_utlization_kernel_samples(const double * __restrict __attribute__((aligned(64))) a,
                                  const double * __restrict __attribute__((aligned(64))) b,
			          double * __restrict c,
			          const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif        
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cpu_utilization_kernel(a[i],
	                                         b[i]);
					
          }
}


/*
    Cycles Per Instruction (CPI).
*/
void
skx_cycles_per_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                             const double * __restrict __attribute__((aligned(64))) b,
			     double * __restrict __attribute__((aligned(64))) c,
			     const int32_t data_len) {

 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif      
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cycles_per_instr(a[i],
	                                   b[i]);
					
          }	 
}


/*
    Cycles Per Instruction (CPI) kernel mode.
*/
void
skx_cycles_per_instr_kernel_samples(const double * __restrict __attribute__((aligned(64))) a,
                                    const double * __restrict __attribute__((aligned(64))) b,
			            double * __restrict __attribute__((aligned(64))) c,
			            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif       
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cycles_per_instr_kernel(a[i],
	                                          b[i]);
					
          }	 
}


/*
    EMON event multiplexing reliability (>95% --- means satisfying ratio).
*/
void
skx_mux_reliability_samples(const double * __restrict  __attribute__((aligned(64))) a,
                            const double * __restrict  __attribute__((aligned(64))) b,
			    double * __restrict  __attribute__((aligned(64))) c,
			    const int32_t data_len) {
			    
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_emon_mux_reliability(a[i],
	                                       b[i]);
					
          }	       
}


/*
      Branch mispredict ratio.
*/
void
skx_branch_mispred_ratio_samples(const double * __restrict __attribute__((aligned(64))) a,
                                 const double * __restrict __attribute__((aligned(64))) b,
			         double * __restrict __attribute__((aligned(64))) c,
			         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_branch_mispredict_ratio(a[i],
	                                          b[i]);
					
          }	
}


/*
    Loads per instruction.
*/
void
skx_loads_per_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                            const double * __restrict  __attribute__((aligned(64))) b,
			    double * __restrict  __attribute__((aligned(64))) c,
			    const int32_t data_len data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_loads_per_instr(a[i],
	                                  b[i]);
					
          }
}


/*
    Stores per instruction.
*/
void
skx_stores_per_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                             const double * __restrict  __attribute__((aligned(64))) b,
			     double * __restrict  __attribute__((aligned(64))) c,
			     const int32_t data_len data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_stores_per_instr(a[i],
	                                  b[i]);
					
          }
}


/*
   Memory operations per instruction.
*/
void
skx_mem_ops_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                          const double * __restrict  __attribute__((aligned(64))) b,
		          const double * __restrict  __attribute__((aligned(64))) c,
		          double * __restrict  __attribute__((aligned(64))) d,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,d,data_len) private(i) \
        aligned(a:64,b,c,d) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_mem_ops_per_instr(a[i],
	                                    b[i],
					    c[i]);
					
          }
}


/*
    Locks retired per instruction.
*/
void
skx_lock_per_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                           const double * __restrict  __attribute__((aligned(64))) b,
			   double * __restrict  __attribute__((aligned(64))) c,
			   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_locks_per_instr(a[i],
	                                  b[i]);
					
          }
}


/*
   Uncacheable reads per instruction.
*/
void
skx_uncacheable_reads_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                    const double * __restrict __attribute__((aligned(64))) b,
			            double * __restrict __attribute__((aligned(64))) c,
			            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_uncacheable_reads_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
    Streaming-stores (full line) per instruction.
*/
void
skx_stream_stores_fl_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                   const double * __restrict __attribute__((aligned(64))) b,
			           double * __restrict __attribute__((aligned(64))) c,
			           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_streaming_stores_fl_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
     Streaming-stores (partial line) per instruction.
*/
void
skx_stream_stores_pl_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                   const double * __restrict __attribute__((aligned(64))) b,
			           double * __restrict __attribute__((aligned(64))) c,
			           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_streaming_stores_pl_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
     L1D$ Misses per instruction including:
     -- data
     -- reads for ownership (rfo) with prefetches
*/
void
skx_L1D_misses_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                             const double * __restrict __attribute__((aligned(64))) b,
			     double * __restrict __attribute__((aligned(64))) c,
			     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L1D_misses_per_instr(a[i],
	                                       b[i]);
					
          }
}


/*
     L1D$ Demand data read hit per instruction.
*/
void
skx_L1D_hits_instr_samples(  const double * __restrict __attribute__((aligned(64))) a,
                             const double * __restrict __attribute__((aligned(64))) b,
			     double * __restrict __attribute__((aligned(64))) c,
			     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L1D_hits_per_instr(a[i],
	                                     b[i]);
					
          }
}


/*
    L1$ code read misses (including prefetches) per instruction.
*/
void
skx_L1I_read_misses_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                  const double * __restrict __attribute__((aligned(64))) b,
			          double * __restrict __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L1I_read_misses_per_instr(a[i],
	                                            b[i]);
					
          }
}


/*
     L2 Demand data read hits per instruction.
*/
void
skx_L2_data_read_misses_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                      const double * __restrict __attribute__((aligned(64))) b,
			              double * __restrict __attribute__((aligned(64))) c,
			              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_data_read_hits_per_instr(a[i],
	                                            b[i]);
					
          }
}


/*
      L2 Misses per instruction including:
     -- code
     -- data
     -- read for ownership and prefetches
*/
void
skx_L2_all_misses_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                const double * __restrict __attribute__((aligned(64))) b,
			        double * __restrict __attribute__((aligned(64))) c,
			        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_all_misses_per_instr(a[i],
	                                          b[i]);
					
          }
}


/*
       L2 demand data read misses per instruction.
*/
void
skx_L2_demand_data_read_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                                    const double * __restrict __attribute__((aligned(64))) b,
			            double * __restrict __attribute__((aligned(64))) c,
			            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_demand_data_read_mpi(a[i],
	                                          b[i]);
					
          }
}


/*
     L2 demand data code misses per instruction.
*/
void
skx_L2_demand_code_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                               const double * __restrict __attribute__((aligned(64))) b,
			       double * __restrict __attribute__((aligned(64))) c,
			       const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_demand_code_mpi(a[i],
	                                     b[i]);
					
          }
}


/*
    L2 Any local request that HITM in a sibling core (per instruction).
*/
void
skx_L2_req_hitm_sibling_core_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                     const double * __restrict  __attribute__((aligned(64))) b,
			             double * __restrict  __attribute__((aligned(64))) c,
			             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_request_hitm_sibling_core(a[i],
	                                               b[i]);
					
          }
}


/*
     L2 (percentage) of all lines evicted that are unused prefetches.
*/
void
skx_L2_lines_evict_unused_prefetch_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                           const double * __restrict  __attribute__((aligned(64))) b,
					   const double * __restrict  __attribute__((aligned(64))) c,
					   const double * __restrict  __attribute__((aligned(64))) d,
					   double * __restrict  __attribute__((aligned(64))) e,
					   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L2_lines_evict_unused_prefetch(a[i],
	                                                 b[i],
							 c[i],
							 d[i]);
					
          }
}


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
void
skx_L2_evict_L3_alloc_samples(const double * __restrict  __attribute__((aligned(64))) a,
                              const double * __restrict  __attribute__((aligned(64))) b,
			      double * __restrict  __attribute__((aligned(64))) c,
			      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_evict_L3_allocated(a[i],
	                                        b[i]);
					
          }
}


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
void
skx_L2_evict_L3_no_alloc_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                 const double * __restrict  __attribute__((aligned(64))) b,
			         double * __restrict  __attribute__((aligned(64))) c,
			         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_evict_L3_not_allocated(a[i],
	                                            b[i]);
					
          }
}


/*
     LLC code references per instruction (L3 prefetch excluded).
*/
void
skx_L3_code_references_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          double * __restrict  __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_code_references_per_instr(a[i],
	                                               b[i]);
					
          }
}


/*
     LLC data read references per instruction (L3 prefetch excluded).
*/
void
skx_L3_data_references_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          double * __restrict  __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_data_references_per_instr(a[i],
	                                               b[i]);
					
          }
}


/*
      LLC RFO references per instrctions (L3 prefetch excluded).
*/
void
skx_L3_rfo_references_instr_samples( const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          double * __restrict  __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_rfo_references_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
     LLC Misses Per Instruction (includes code and data and rfo with prefetches).
*/
void
skx_L3_all_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
               const double * __restrict __attribute__((aligned(64))) b,
	       const double * __restrict __attribute__((aligned(64))) c,
               const double * __restrict __attribute__((aligned(64))) d,
	       double * __restrict __attribute__((aligned(64))) e,
	       const int32_t data_len)  {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_all_mpi(a[i],
	                             b[i],
				     c[i],
				     d[i]);
					
          }
}


/*
     LLC data read Misses Per Instruction (demand and prefetch).
*/
void
skx_L3_data_read_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
			  double * __restrict  __attribute__((aligned(64))) c,
			  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_data_read_mpi(a[i],
	                                   b[i]);
					
          }
}


/*
      LLC RFO read Misses Per Instruction (demand and prefetch).  
*/
void
skx_L3_rfo_data_read_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                              const double * __restrict __attribute__((aligned(64))) b,
			      double * __restrict  __attribute__((aligned(64))) c,
			      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_rfo_read_mpi(a[i],
	                                  b[i]);
					
          }
}





