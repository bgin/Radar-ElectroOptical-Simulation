

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


/*
     LLC code read MPI (demand and prefetch).
*/
void
skx_L3_total_hitm_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                const double * __restrict  __attribute__((aligned(64))) b,
			        double * __restrict   __attribute__((aligned(64))) c,
			        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_total_hitm_per_instr(a[i],
	                                          b[i]);
					
          }
}


/*
     LLC total HIT clean line forwards (per instr) (excludes LLC prefetches)
*/
void
skx_L3_total_hitm_clean_lines_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                            const double * __restrict __attribute__((aligned(64))) b,
			                    double * __restrict __attribute__((aligned(64))) c,
			                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_total_hitm_clean_lines_per_instr(a[i],
	                                                      b[i]);
					
          }
}


/*
    Average LLC data read (demand+prefetch) miss latency (in ns).  
*/
void
skx_L3_avg_data_read_ns_samples(const double * __restrict __attribute__((aligned(64))) a,
                                const double * __restrict __attribute__((aligned(64))) b,
			        const double * __restrict __attribute__((aligned(64))) c,
                                const double * __restrict  __attribute__((aligned(64))) d,
			        double * __restrict __attribute__((aligned(64))) e,
			        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_avg_data_read_ns(a[i],
	                                      b[i],
				              c[i],
				              d[i]);
					
          }
}


/*
      Average LLC data read (demand and prefetch) miss latency (in UNCORE clk).
*/
void
skx_L3_avg_data_read_unc_clk_samples(const double * __restrict __attribute__((aligned(64))) a,
                                     const double * __restrict __attribute__((aligned(64))) b,
			             double * __restrict __attribute__((aligned(64))) c,
			             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_avg_data_read_unc_clk(a[i],
	                                           b[i]);
					
          }
}


/*
      Average LLC data read (demand+prefetch) miss latency for LOCAL requests (in ns)
*/
void
skx_L3_avg_data_read_loc_ns_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                    const double * __restrict  __attribute__((aligned(64))) b,
			            const double * __restrict  __attribute__((aligned(64))) c,
                                    const double * __restrict  __attribute__((aligned(64))) d,
			            double * __restrict  __attribute__((aligned(64))) e,
			            const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_avg_data_read_loc_req_ns(a[i],
	                                              b[i],
				                      c[i],
				                      d[i]);
					
          }
}


/*
    Average LLC data read (demand+prefetch) miss latency for REMOTE requests (in ns).
*/
void
skx_L3_avg_data_read_rem_ns_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                    const double * __restrict  __attribute__((aligned(64))) b,
			            const double * __restrict  __attribute__((aligned(64))) c,
                                    const double * __restrict  __attribute__((aligned(64))) d,
			            double * __restrict  __attribute__((aligned(64))) e,
			            const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_avg_data_read_rem_req_ns(a[i],
	                                              b[i],
				                      c[i],
				                      d[i]);
					
          }
}


/*
      Average LLC data read (demand+prefetch) miss latency  for LOCAL requests (in UNCORE clk)
*/
void
skx_L3_avg_data_read_loc_unc_clk_samples(const double * __restrict __attribute__((aligned(64))) a,
                                         const double * __restrict __attribute__((aligned(64))) b,
			                 double * __restrict __attribute__((aligned(64))) c,
			                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_avg_data_read_loc_unc_clk(a[i],
	                                               b[i]);
					
          }
}


/*
     Average LLC data read (demand+prefetch) miss latency  for REMOTE requests (in UNCORE clk)
*/
void
skx_L3_avg_data_read_rem_unc_clk_samples(const double * __restrict __attribute__((aligned(64))) a,
                                         const double * __restrict __attribute__((aligned(64))) b,
			                 double * __restrict __attribute__((aligned(64))) c,
			                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_avg_data_read_rem_unc_clk(a[i],
	                                               b[i]);
					
          }
}


/*
    ITLB MPI
*/
void
skx_ITLB_mpi_data_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_ITLB_mpi(a[i],
	                           b[i]);
					
          }
}


/*
   ITLB large page MPI
*/
void
skx_ITLB_2M4M_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_ITLB_2M_4M_mpi(a[i],
	                                 b[i]);
					
          }
}


/*
     DTLB load MPI.
*/
void
skx_DTLB_load_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_load_mpi(a[i],
	                                b[i]);
					
          }
}


/*
    DTLB large page load MPI.
*/
void
skx_DTLB_2M4M_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_2M_4M_mpi(a[i],
	                                b[i]);
					
          }
}


/*
   DTLB store MPI.
*/
void
skx_DTLB_store_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_store_mpi(a[i],
	                                b[i]);
					
          }
}


/*
    DTLB load miss latency (in core clks)
*/
void
skx_DTLB_load_miss_clks_samples(const double * __restrict  __attribute__((aligned(64))) a,
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
               c[i] = skx_DTLB_load_miss_clks(a[i],
	                                      b[i]);
					
          }
}


/*
     DTLB store miss latency (in core clks).
*/
void
skx_DTLB_store_miss_clks_samples(const double * __restrict  __attribute__((aligned(64))) a,
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
               c[i] = skx_DTLB_store_miss_clks(a[i],
	                                      b[i]);
					
          }
}


/*
      ITLB miss latency (in core clks).
*/
void
skx_ITLB_miss_latency_clks_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                   const double * __restrict   __attribute__((aligned(64))) b,
		                   double * __restrict  __attribute__((aligned(64))) c,
		                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_ITLB_miss_latency_clks(a[i],
	                                      b[i]);
					
          }
}


/*
    NUMA percentage of Reads addressed to local DRAM.
*/
void
skx_numa_reads_local_dram_samples(const double * __restrict   __attribute__((aligned(64))) a,
                                  const double * __restrict   __attribute__((aligned(64))) b,
		                  double * __restrict  __attribute__((aligned(64))) c,
		                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_nume_reads_local_dram(a[i],
	                                        b[i]);
					
          }
}


/*
    NUMA percentage of Reads addressed to remote  DRAM.
*/
void
skx_numa_reads_remote_dram_samples(const double * __restrict   __attribute__((aligned(64))) a,
                                  const double * __restrict   __attribute__((aligned(64))) b,
		                  double * __restrict  __attribute__((aligned(64))) c,
		                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_nume_reads_remote_dram(a[i],
	                                        b[i]);
					
          }
}


/*
     Uncore Frequency Ghz.
*/
void
skx_uncore_frequency_samples(const double * __restrict  __attribute__((aligned(64))) a,
                             const double * __restrict  __attribute__((aligned(64))) b,
			     const double * __restrict  __attribute__((aligned(64))) c,
		             double * __restrict d,
		             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_uncore_frequency_ghz(a[i],
	                                       b[i],
					       c[i]);
					
          }
}


/*
    UPI speed - GT/s.
*/
void
skx_UPI_speed_samples(const double * __restrict __attribute__((aligned(64))) a,
                      const double * __restrict __attribute__((aligned(64))) b,
	              double * __restrict __attribute__((aligned(64))) c,
	              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_UPI_speed(a[i],
	                            b[i]);
					
          }
}


/*
    UPI Data transmit BW (MB/sec) (only data)
*/
void
skx_UPI_data_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                        const double * __restrict __attribute__((aligned(64))) b,
	                double * __restrict __attribute__((aligned(64))) c,
	                const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_UPI_data_bw(a[i],
	                            b[i]);
					
          }
}


/*
     UPI Total transmit BW (MB/sec) (includes control)
*/
void
skx_UPI_tx_total_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                            const double * __restrict __attribute__((aligned(64))) b,
			    const double * __restrict __attribute__((aligned(64))) c,
	                    double * __restrict __attribute__((aligned(64))) d,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(10) if(data_len>=10000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_UPI_transmit_total_bw(a[i],
	                                       b[i],
					       c[i]);
					
          }
}
