
#ifndef __GMS_SPR_METRICS_YMM8R4_HPP__
#define __GMS_SPR_METRICS_YMM8R4_HPP__


#include <immintrin.h>
#include "GMS_config.h"


namespace gms {


/*
        "BriefDescription": "CPU operating frequency (in GHz)",
        "MetricExpr": "( CPU_CLK_UNHALTED.THREAD / CPU_CLK_UNHALTED.REF_TSC * #SYSTEM_TSC_FREQ ) / 1000000000",
        "MetricGroup": "",
        "MetricName": "cpu_operating_frequency",
        "ScaleUnit": "1GHz"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_cpu_op_freq_ymm8r4(const __m256 CPU_CLK_UNHALTED_THREAD,
	                                         const __m256 CPU_CLK_UNHALTED_REF_TSC,
	                                         const __m256 SYSTEM_TSC_FREQ) {
	                  
	                  const __m256 C000000001 = _mm256_set1_ps(1.0e-09f);
	                  register __m256 t0,metric;
	                  t0 = _mm256_div_ps(CPU_CLK_UNHALTED_THREAD,
	                                     CPU_CLK_UNHALTED_REF_TSC);
	                  metric = _mm256_mul_ps(C000000001,
	                                     _mm256_mul_ps(t0,SYSTEM_TSC_FREQ));
	                  return (metric);
	                                             
	          }
	          
	          
/*
        "BriefDescription": "Percentage of time spent in the active CPU power state C0",
        "MetricExpr": "CPU_CLK_UNHALTED.REF_TSC / TSC",
        "MetricGroup": "",
        "MetricName": "cpu_utilization",
        "ScaleUnit": "100%"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_cpu_util_ymm8r4(const __m256 CPU_CLK_UNHALTED_REF_TSC,
	                                      const __m256 TSC_STAMP) {
	                                      
	                  register __m256 metric;
	                  result = _mm256_div_ps(CPU_CLK_UNHALTED_REF_TSC,TSC_STAMP);
	                  return (metric);                             
	        }  
	        
/*
         "BriefDescription": "Cycles per instruction retired; indicating how much time each executed instruction took; in units of cycles.",
        "MetricExpr": "CPU_CLK_UNHALTED.THREAD / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "cpi",
        "ScaleUnit": "1per_instr"
*/      

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_cpi_ymm8r4(const __m256 CPU_CLK_UNHALTED_THREAD,
	                                 const __m256 INST_RETIRED_ANY) {
	                                 
	                  register __m256 metric;
	                  metric = _mm256_div_ps(CPU_CLK_UNHALTED_THREAD,
	                                         INST_RETIRED_ANY);
	                  return (metric);                       
	         }
	         
/*
        "BriefDescription": "The ratio of number of completed memory load instructions to the total number completed instructions",
        "MetricExpr": "MEM_INST_RETIRED.ALL_LOADS / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "loads_per_instr",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_loads_per_instr_ymm8r4(const __m256 MEM_INST_RETIRED_ALL_LOADS,
	                                             const __m256 INST_RETIRED_ANY) {
	                   
	                   register __m256 metric;
	                   metric = _mm256_div_ps(MEM_INST_RETIRED_ALL_LOADS,
	                                          INST_RETIRED_ANY);
	                   return (metric);
	                                                                       
	          }
	          
/*
        "BriefDescription": "The ratio of number of completed memory store instructions to the total number completed instructions",
        "MetricExpr": "MEM_INST_RETIRED.ALL_STORES / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "stores_per_instr",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_stores_per_instr_ymm8r4(const __m256 MEM_INST_RETIRED_ALL_STORES,
	                                              const __m256 INST_RETIRED_ANY) {
	                   
	                   register __m256 metric;
	                   metric = _mm256_div_ps(MEM_INST_RETIRED_ALL_STORES,
	                                          INST_RETIRED_ANY);
	                   return (metric);
	                                                                       
	          }
	          
/*
        "BriefDescription": "Ratio of number of requests missing L1 data cache (includes data+rfo w/ prefetches) to the total number of completed instructions",
        "MetricExpr": "L1D.REPLACEMENT / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "l1d_mpi",
        "ScaleUnit": "1per_instr" 
*/

                   
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1d_mpi_ymm8r4(const __m256 L1D_REPLACEMENT,
	                                     const __m256 INST_RETIRED_ANY) {
	                                     
	                  register __m256 metric;
	                  metric = _mm256_div_ps(L1D_REPLACEMENT,
	                                         INST_RETIRED_ANY);
	                  return (metric);                            
	        }  
	        
/*
        "BriefDescription": "Ratio of number of demand load requests hitting in L1 data cache to the total number of completed instructions ",
        "MetricExpr": "MEM_LOAD_RETIRED.L1_HIT / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "l1d_demand_data_read_hits_per_instr",
        "ScaleUnit": "1per_instr"
*/ 

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1d_ddrh_pi_ymm8r4(const __m256 MEM_LOAD_RETIRED_L1_HIT,
	                                         const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(MEM_LOAD_RETIRED_L1_HIT,
	                                         INST_RETIRED_ANY);
	                  return (metric);                       
	                                                      
	         }
	         
/*
         "BriefDescription": "Ratio of number of code read requests missing in L1 instruction cache (includes prefetches) to the total number of completed instructions",
        "MetricExpr": "L2_RQSTS.ALL_CODE_RD / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "l1_i_code_read_misses_with_prefetches_per_instr",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1i_crmp_pi_ymm8r4(const __m256 L2_RQSTS_ALL_CODE_RD,
	                                         const __m256 INST_RETIRED_ANY) {
	                                         
	                   register __m256 metric;
	                   metric = _mm256_div_ps(L2_RQSTS_ALL_CODE_RD,
	                                          INST_RETIRED_ANY);
	                   return (metric);                               
	         }
	         
/*
         "BriefDescription": "Ratio of number of completed demand load requests hitting in L2 cache to the total number of completed instructions ",
        "MetricExpr": "MEM_LOAD_RETIRED.L2_HIT / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "l2_demand_data_read_hits_per_instr",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_ddrh_pi_ymm8r4(const __m256 MEM_LOAD_RETIRED_L2_HIT,
	                                        const __m256 INST_RETIRED_ANY) {
	                   
	                  register __m256 metric;
	                  metric = _mm256_div_ps(MEM_LOAD_RETIRED_L2_HIT,
	                                         INST_RETIRED_ANY);
	                  return (metric);                             
	          }
	          
/*
        "BriefDescription": "Ratio of number of requests missing L2 cache (includes code+data+rfo w/ prefetches) to the total number of completed instructions",
        "MetricExpr": "L2_LINES_IN.ALL / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "l2_mpi",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_mpi_ymm8r4(const __m256 L2_LINES_IN_ALL,
	                                    const __m256 INST_RETIRED_ANY) {
	                                    
	                  register __m256 metric;
	                  metric = _mm256_div_ps(L2_LINES_IN_ALL,
	                                         INST_RETIRED_ANY);
	                  return (metric);                           
	        }
	        
/*
        "BriefDescription": "Ratio of number of completed data read request missing L2 cache to the total number of completed instructions",
        "MetricExpr": "MEM_LOAD_RETIRED.L2_MISS / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "l2_demand_data_read_mpi",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_ddr_mpi_ymm8r4(const __m256 MEM_LOAD_RETIRED_L2_MISS,
	                                        const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(MEM_LOAD_RETIRED_L2_MISS,
	                                         INST_RETIRED_ANY);
	                  return (metric);                               
	         }  
	         
/*
        "BriefDescription": "Ratio of number of code read request missing L2 cache to the total number of completed instructions",
        "MetricExpr": "L2_RQSTS.CODE_RD_MISS / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "l2_demand_code_mpi",
        "ScaleUnit": "1per_instr"
*/     

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_dc_mpi_ymm8r4(const __m256 L2_RQSTS_CODE_RD_MISS,
	                                       const __m256 INST_RETIRED_ANY) {
	                          
	                  register __m256 metric;
	                  metric = _mm256_div_ps(L2_RQSTS_CODE_RD_MISS,
	                                         INST_RETIRED_ANY);
	                  return (metric);             
	          }    
	          
/*
        "BriefDescription": "Ratio of number of data read requests missing last level core cache (includes demand w/ prefetches) to the total number of completed instructions",
        "MetricExpr": "( UNC_CHA_TOR_INSERTS.IA_MISS_LLCPREFDATA + UNC_CHA_TOR_INSERTS.IA_MISS_DRD + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_PREF ) / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "llc_data_read_mpi_demand_plus_prefetch",
        "ScaleUnit": "1per_instr"   
*/    

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_drmpi_dpp_ymm8r4(const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA,
	                                           const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD,
	                                           const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF,
	                                           const __m256 INST_RETIRED_ANY) {
	                                           
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA,
	                                 _mm256_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD,
	                                               UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF));
	                  metric = _mm256_div_ps(t0,INST_RETIRED_ANY);
	                  return (metric);                               
	         } 
	         
/*
        "BriefDescription": "Ratio of number of code read requests missing last level core cache (includes demand w/ prefetches) to the total number of completed instructions",
        "MetricExpr": "UNC_CHA_TOR_INSERTS.IA_MISS_CRD / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "llc_code_read_mpi_demand_plus_prefetch",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_crmpi_dpp_ymm8r4(const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_CRD,
	                                           const __m256 INST_RETIRED_ANY) {
	                   
	                   register __m256 metric;
	                   metric = _mm256_div_ps(UNC_CHA_TOR_INSERTS_IA_MISS_CRD,
	                                          INST_RETIRED_ANY);
	                   return (metric);                                 
	         }
	      
/*
        "BriefDescription": "Average latency of a last level cache (LLC) demand data read miss (read memory access) in nano seconds",
        "MetricExpr": "( 1000000000 * ( UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD / UNC_CHA_TOR_INSERTS.IA_MISS_DRD ) / ( UNC_CHA_CLOCKTICKS / ( source_count(UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD) *   #num_packages ) ) ) * duration_time",
        "MetricGroup": "",
        "MetricName": "llc_demand_data_read_miss_latency",
        "ScaleUnit": "1ns"
*/

          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddrm_lat_ymm8r4(const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD,
                                                  const __m256 UNC_CHA_CLOCKTICKS,
                                                  const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const float duration_time) {
                                                 
                          const __m256 C1E9 = _mm256_set1_ps(1.0e+09f); 
                          register __m256 t0;
                          register __m256 t1;
                          register __m256 vdt;
                          register __m256 metric;
                        
                          t0  = _mm256_mul_ps(C1E9,
                                          _mm256_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD)); 
                          vdt = _mm256_set1_ps(duration_time);
                          t1  = _mm256_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm256_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD));
                          metric  = _mm256_mul_ps(_mm256_div_ps(t0,t1),vdt);
                          return (metric);
                  }    
                  
/*
        "BriefDescription": "Average latency of a last level cache (LLC) demand data read miss (read memory access) addressed to local memory in nano seconds",
        "MetricExpr": "( 1000000000 * ( UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD_LOCAL / UNC_CHA_TOR_INSERTS.IA_MISS_DRD_LOCAL ) / ( UNC_CHA_CLOCKTICKS / ( source_count(UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD_LOCAL) * #num_packages ) ) ) * duration_time",
        "MetricGroup": "",
        "MetricName": "llc_demand_data_read_miss_latency_for_local_requests",
        "ScaleUnit": "1ns"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddr_mllr_ymm8r4(const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                          const __m256 UNC_CHA_CLOCKTICKS,
	                                          const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const float duration_time) {
	                                          
	                  const __m256 C1E9 = _mm256_set1_ps(1.0e+09f); 
                          register __m256 t0;
                          register __m256 t1;
                          register __m256 vdt;
                          register __m256 metric;  
                          t0  = _mm256_mul_ps(C1E9,
                                          _mm256_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL)); 
                          vdt = _mm256_set1_ps(duration_time);
                          t1  = _mm256_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm256_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL));
                          metric  = _mm256_mul_ps(_mm256_div_ps(t0,t1),vdt);
                          return (metric);                            
	          } 
	          
/*
        "BriefDescription": "Average latency of a last level cache (LLC) demand data read miss (read memory access) addressed to remote memory in nano seconds",
        "MetricExpr": "( 1000000000 * ( UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD_REMOTE / UNC_CHA_TOR_INSERTS.IA_MISS_DRD_REMOTE ) / ( UNC_CHA_CLOCKTICKS / ( source_count(UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD_REMOTE) * #num_packages ) ) ) * duration_time",
        "MetricGroup": "",
        "MetricName": "llc_demand_data_read_miss_latency_for_remote_requests",
        "ScaleUnit": "1ns"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddr_mlrr_ymm8r4(const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                          const __m256 UNC_CHA_CLOCKTICKS,
	                                          const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const float duration_time) {
	                                          
	                  const __m256 C1E9 = _mm256_set1_ps(1.0e+09f); 
                          register __m256 t0;
                          register __m256 t1;
                          register __m256 vdt;
                          register __m256 metric;  
                          t0  = _mm256_mul_ps(C1E9,
                                          _mm256_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE)); 
                          vdt = _mm256_set1_ps(duration_time);
                          t1  = _mm256_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm256_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE));
                          metric  = _mm256_mul_ps(_mm256_div_ps(t0,t1),vdt);
                          return (metric);  
                }
                
/*
        "BriefDescription": "Average latency of a last level cache (LLC) demand data read miss (read memory access) addressed to DRAM in nano seconds",
        "MetricExpr": "( 1000000000 * ( UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD_DDR / UNC_CHA_TOR_INSERTS.IA_MISS_DRD_DDR ) / ( UNC_CHA_CLOCKTICKS / ( source_count(UNC_CHA_TOR_OCCUPANCY.IA_MISS_DRD_DDR) * #num_packages ) ) ) * duration_time",
        "MetricGroup": "",
        "MetricName": "llc_demand_data_read_miss_to_dram_latency",
        "ScaleUnit": "1ns"
*/        

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddr_mdlat_ymm8r4(const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR,
	                                           const __m256 UNC_CHA_CLOCKTICKS,
	                                           const __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const float duration_time) {
	                                          
	                  const __m256 C1E9 = _mm256_set1_ps(1.0e+09f); 
                          register __m256 t0;
                          register __m256 t1;
                          register __m256 vdt;
                          register __m256 metric;  
                          t0  = _mm256_mul_ps(C1E9,
                                          _mm256_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR)); 
                          vdt = _mm256_set1_ps(duration_time);
                          t1  = _mm256_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm256_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR));
                          metric  = _mm256_mul_ps(_mm256_div_ps(t0,t1),vdt);
                          return (metric);  
                }  
                
/*
         "BriefDescription": "Ratio of number of completed page walks (for all page sizes) caused by a code fetch to the total number of completed instructions. This implies it missed in the ITLB (Instruction TLB) and further levels of TLB.",
        "MetricExpr": "ITLB_MISSES.WALK_COMPLETED / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "itlb_2nd_level_mpi",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_itlb_2l_mpi_ymm8r4(const __m256 ITLB_MISSES_WALK_COMPLETED,
	                                         const __m256 INST_RETIRED_ANY) {
	                                         
	                  register __m256 metric;
	                  metric = _mm256_div_ps(ITLB_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                              
	          }
	          
/*
        "BriefDescription": "Ratio of number of completed page walks (for 2 megabyte and 4 megabyte page sizes) caused by a code fetch to the total number of completed instructions. This implies it missed in the Instruction Translation Lookaside Buffer (ITLB) and further levels of TLB.",
        "MetricExpr": "ITLB_MISSES.WALK_COMPLETED_2M_4M / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "itlb_2nd_level_large_page_mpi",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_itlb_2l_llp_mpi_ymm8r4(const __m256 ITLB_MISSES_WALK_COMPLETED_2M_4M,
	                                             const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(ITLB_MISSES_WALK_COMPLETED_2M_4M,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }
	          
/*
        "BriefDescription": "Ratio of number of completed page walks (for all page sizes) caused by demand data loads to the total number of completed instructions. This implies it missed in the DTLB and further levels of TLB.",
        "MetricExpr": "DTLB_LOAD_MISSES.WALK_COMPLETED / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "dtlb_2nd_level_load_mpi",
        "ScaleUnit": "1per_instr"
*/

                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_llmpi_ymm8r4(const __m256 DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }
	          
/*
         "BriefDescription": "Ratio of number of completed page walks (for 2 megabyte page sizes) caused by demand data loads to the total number of completed instructions. This implies it missed in the Data Translation Lookaside Buffer (DTLB) and further levels of TLB.",
        "MetricExpr": "DTLB_LOAD_MISSES.WALK_COMPLETED_2M_4M / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "dtlb_2nd_level_2mb_large_page_load_mpi",
        "ScaleUnit": "1per_instr"
*/

             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_2mblp_mpi_ymm8r4(const __m256 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                               const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }
	          
/*
    "BriefDescription": "Ratio of number of completed page walks (for all page sizes) caused by demand data loads to the total number of completed instructions. This implies it missed in the DTLB and further levels of TLB.",
        "MetricExpr": "DTLB_LOAD_MISSES.WALK_COMPLETED / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "dtlb_2nd_level_load_mpi",
        "ScaleUnit": "1per_instr"    
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_llmpi_ymm8r4(const __m256 DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }  
	          
/*
      "BriefDescription": "Ratio of number of completed page walks (for 2 megabyte page sizes) caused by demand data loads to the total number of completed instructions. This implies it missed in the Data Translation Lookaside Buffer (DTLB) and further levels of TLB.",
        "MetricExpr": "DTLB_LOAD_MISSES.WALK_COMPLETED_2M_4M / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "dtlb_2nd_level_2mb_large_page_load_mpi",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_2mllmpi_ymm8r4(const __m256 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                             const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }  
	          
/*
       "BriefDescription": "Ratio of number of completed page walks (for all page sizes) caused by demand data stores to the total number of completed instructions. This implies it missed in the DTLB and further levels of TLB.",
        "MetricExpr": "DTLB_STORE_MISSES.WALK_COMPLETED / INST_RETIRED.ANY",
        "MetricGroup": "",
        "MetricName": "dtlb_2nd_level_store_mpi",
        "ScaleUnit": "1per_instr"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_store_mpi_ymm8r4(const __m256 DTLB_STORE_MISSES_WALK_COMPLETED,
	                                               const __m256 INST_RETIRED_ANY) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(DTLB_STORE_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }  
	          
/*
        "BriefDescription": "Memory read that miss the last level cache (LLC) addressed to local DRAM as a percentage of total memory read accesses, does not include LLC prefetches.",
        "MetricExpr": "( UNC_CHA_TOR_INSERTS.IA_MISS_DRD_LOCAL + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_PREF_LOCAL ) / ( UNC_CHA_TOR_INSERTS.IA_MISS_DRD_LOCAL + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_PREF_LOCAL + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_REMOTE + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_PREF_REMOTE )",
        "MetricGroup": "",
        "MetricName": "numa_reads_addressed_to_local_dram",
        "ScaleUnit": "100%"
*/

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_numa_reads_ald_ymm8r4(const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL,
	                                            const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE) {
	               
	                  register __m256 t0,t1,t2;
	                  register __m256 metric;
	                  t0 = _mm256_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL);
	                  t1 = _mm256_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE);
	                  t2 = _mm256_add_ps(t0,t1);
	                  metric = _mm256_div_ps(t0,t2);
	                  return (metric);                             
	         }
	         
/*
          "BriefDescription": "Memory reads that miss the last level cache (LLC) addressed to remote DRAM as a percentage of total memory read accesses, does not include LLC prefetches.",
        "MetricExpr": "( UNC_CHA_TOR_INSERTS.IA_MISS_DRD_REMOTE + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_PREF_REMOTE ) / ( UNC_CHA_TOR_INSERTS.IA_MISS_DRD_LOCAL + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_PREF_LOCAL + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_REMOTE + UNC_CHA_TOR_INSERTS.IA_MISS_DRD_PREF_REMOTE )",
        "MetricGroup": "",
        "MetricName": "numa_reads_addressed_to_remote_dram",
        "ScaleUnit": "100%"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_numa_reads_ard_ymm8r4(const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE,
	                                            const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL) {
	               
	                  register __m256 t0,t1,t2;
	                  register __m256 metric;
	                  t0 = _mm256_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE);
	                  t1 = _mm256_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL);
	                  t2 = _mm256_add_ps(t0,t1);
	                  metric = _mm256_div_ps(t0,t2);
	                  return (metric);                             
	         }    
	         
/*
        "BriefDescription": "Uncore operating frequency in GHz",
        "MetricExpr": "( UNC_CHA_CLOCKTICKS / ( source_count(UNC_CHA_CLOCKTICKS) * #num_packages ) / 1000000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "uncore_frequency",
        "ScaleUnit": "1GHz"
*/ 

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uncore_freq_ymm8r4(const __m256 UNC_CHA_CLOCKTICKS,
	                                     const float n_scks,
	                                     const float idt) { // shall be inverse value passed.
	                 
	                 const __m256 C000000001 = _mm256_set1_ps(1.0e-09f);
	                 register  __m256 vscks;
	                 register  __m256 vidt;
	                 register  __m256 t0;
	                 register  __m256 t1;
	                 register  __m256 metric;
	                 vscks = _mm256_set1_ps(n_scks);
	                 vdt   = _mm256_set1_ps(dt);
	                 t0    = _mm256_mul_ps(C000000001,
	                                   _mm256_mul_ps(vscks,UNC_CHA_CLOCKTICKS));
	                 t1    = _mm256_div_ps(UNC_CHA_CLOCKTICKS,t0);
	                 metric= _mm256_mul_ps(t1,vidt);
	                 return (metric); 
	         }
	         
/*
         "BriefDescription": "Intel(R) Ultra Path Interconnect (UPI) data transmit bandwidth (MB/sec)",
        "MetricExpr": "( UNC_UPI_TxL_FLITS.ALL_DATA * (64 / 9.0) / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "upi_data_transmit_bw",
        "ScaleUnit": "1MB/s" 
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_upi_data_tx_bw_ymm8r4(const __m256 UNC_UPI_TxL_FLITS_ALL_DATA,
	                                        const float idt) {
	             
	                  const __m256 C000001                     =
	                                    _mm256_set1_ps(1.0e-06f);
	                  const __m256 C711111111111111111111111   =
	                                    _mm256_set1_ps(7.11111111111111111111111f);
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(UNC_UPI_TxL_FLITS_ALL_DATA,
	                                   _mm256_mul_ps(C711111111111111111111111,
	                                                                        C000001));
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);                      
	          }
	          
/*
        "BriefDescription": "DDR memory read bandwidth (MB/sec)",
        "MetricExpr": "( UNC_M_CAS_COUNT.RD * 64 / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "memory_bandwidth_read",
        "ScaleUnit": "1MB/s"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mem_bw_rd_ymm8r4(const __m256 UNC_M_CAS_COUNT_RD,
	                                       const float idt) {
	                
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                            UNC_M_CAS_COUNT_RD,C640),
	                                                                 C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          


} // gms








#endif /*__GMS_SPR_METRICS_YMM8R4*/
