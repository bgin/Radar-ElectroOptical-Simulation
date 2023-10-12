
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_cpu_op_freq_ymm8r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                         const float * __restrict pCPU_CLK_UNHALTED_REF_TSC,
	                                         const float * __restrict pSYSTEM_TSC_FREQ) {
	                  
	                  register __m256 CPU_CLK_UNHALTED_THREAD   =  
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  register __m256 CPU_CLK_UNHALTED_REF_TSC  =  
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_REF_TSC[0]);
	                  register __m256 SYSTEM_TSC_FREQ           =  
	                                        _mm256_loadu_ps(&pSYSTEM_TSC_FREQ[0]);
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
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_cpu_util_ymm8r4(const float * __restrict pCPU_CLK_UNHALTED_REF_TSC,
	                                      const float * __restrict pTSC_STAMP) {
	                                      
	                  register __m256 CPU_CLK_UNHALTED_REF_TSC = 
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_REF_TSC[0]);
	                  register __m256 TSC_STAMP                =
	                                        _mm256_loadu_ps(&pTSC_STAMP[0]);   
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
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_cpi_ymm8r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                 const float * __restrict pINST_RETIRED_ANY) {
	                       
	                  register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                         _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  register __m256 INST_RETIRED_ANY        = 
	                                         _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);          
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_loads_per_instr_ymm8r4(const float * __restrict pMEM_INST_RETIRED_ALL_LOADS,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                   
	                   register __m256 MEM_INST_RETIRED_ALL_LOADS = 
	                                           _mm256_loadu_ps(&pMEM_INST_RETIRED_ALL_LOADS[0]);
	                   register __m256 INST_RETIRED_ANY           =
	                                           _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_stores_per_instr_ymm8r4(const float * __restrict pMEM_INST_RETIRED_ALL_STORES,
	                                              const float * __restrict pINST_RETIRED_ANY) {
	                   
	                   register __m256 MEM_INST_RETIRED_ALL_STORES = 
	                                           _mm256_loadu_ps(&pMEM_INST_RETIRED_ALL_STORES[0]);
	                   register __m256 INST_RETIRED_ANY            =
	                                           _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	        
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1d_mpi_ymm8r4(const float * __restrict pL1D_REPLACEMENT,
	                                     const float * __restrict pINST_RETIRED_ANY) {
	                                     
	                  register __m256 L1D_REPLACEMENT  = 
	                                     _mm256_loadu_ps(&pL1D_REPLACEMENT[0]);
	                  register __m256 INST_RETIRED_ANY =
	                                     _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1d_ddrh_pi_ymm8r4(const float * __restrict pMEM_LOAD_RETIRED_L1_HIT,
	                                         const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 MEM_LOAD_RETIRED_L1_HIT = 
	                                        _mm256_loadu_ps(&pMEM_LOAD_RETIRED_L1_HIT[0]);
	                  register __m256 INST_RETIRED_ANY        =
	                                        _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1i_crmp_pi_ymm8r4(const float * __restrict pL2_RQSTS_ALL_CODE_RD,
	                                         const float * __restrict pINST_RETIRED_ANY) {
	                                 
	                   register __m256 L2_RQSTS_ALL_CODE_RD = 
	                                           _mm256_loadu_ps(&pL2_RQSTS_ALL_CODE_RD[0]);
	                   register __m256 INST_RETIRED_ANY     =
	                                           _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);        
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_ddrh_pi_ymm8r4(const float * __restrict pMEM_LOAD_RETIRED_L2_HIT,
	                                        const float * __restrict pINST_RETIRED_ANY) {
	                   
	                  register __m256 MEM_LOAD_RETIRED_L2_HIT = 
	                                          _mm256_loadu_ps(&pMEM_LOAD_RETIRED_L2_HIT[0]);
	                  register __m256 INST_RETIRED_ANY        = 
	                                          _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_mpi_ymm8r4(const float * __restrict pL2_LINES_IN_ALL,
	                                    const float * __restrict pINST_RETIRED_ANY) {
	                     
	                  register __m256 L2_LINES_IN_ALL  = 
	                                      _mm256_loadu_ps(&pL2_LINES_IN_ALL[0]);
	                  register __m256 INST_RETIRED_ANY =
	                                      _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);               
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
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_ddr_mpi_ymm8r4(const float * __restrict pMEM_LOAD_RETIRED_L2_MISS,
	                                        const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 MEM_LOAD_RETIRED_L2_MISS = 
	                                          _mm256_loadu_ps(&pMEM_LOAD_RETIRED_L2_MISS[0]);
	                  register __m256 INST_RETIRED_ANY         =
	                                          _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l2_dc_mpi_ymm8r4(const float * __restrict pL2_RQSTS_CODE_RD_MISS,
	                                       const float * __restrict pINST_RETIRED_ANY) {
	                          
	                  register __m256 L2_RQSTS_CODE_RD_MISS = 
	                                          _mm256_loadu_ps(&pL2_RQSTS_CODE_RD_MISS[0]);
	                  register __m256 INST_RETIRED_ANY      =
	                                          _mm256_loadu_ps(pINST_RETIRED_ANY[0]);
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
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_drmpi_dpp_ymm8r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA,
	                                           const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD,
	                                           const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                            
	                  restrict __m256 UNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA = 
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA[0]);
	                  restrict __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD         =
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF    =
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF[0]);
	                  register __m256 INST_RETIRED_ANY                        =
	                                             _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);              
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
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_crmpi_dpp_ymm8r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_CRD,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                   
	                   register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_CRD = 
	                                               _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_CRD[0]);
	                   register __m256 INST_RETIRED_ANY                = 
	                                               _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddrm_lat_ymm8r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD,
                                                  const float * __restrict pUNC_CHA_CLOCKTICKS,
                                                  const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const float duration_time) {
                          
                          register __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD = 
                                                     _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD[0]);
                          register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD   =
                                                     _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD[0]);
                          register __m256 UNC_CHA_CLOCKTICKS                =
                                                     _mm256_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
                          register __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD =
                                                     _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD[0]);                      
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddr_mllr_ymm8r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                          const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                          const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const float duration_time) {
	                      
	                  register __m256  UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL = 
	                                              _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL[0]);
	                  register __m256  UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL   =
	                                              _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL[0]);
	                  register __m256  UNC_CHA_CLOCKTICKS                      =
	                                              _mm256_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
	                  register __m256  UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL =
	                                              _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL[0]);                   
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
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddr_mlrr_ymm8r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                          const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                          const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const float duration_time) {
	                     
	                  register __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE = 
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE   =
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE[0]);
	                  register __m256 UNC_CHA_CLOCKTICKS                       =
	                                             _mm256_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
	                  register __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE =
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE[0]);
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
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_ddr_mdlat_ymm8r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR,
	                                           const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                           const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const float duration_time) {
	                          
	                  register __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR = 
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR   =
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR[0]);
	                  register __m256 UNC_CHA_CLOCKTICKS                    =
	                                             _mm256_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
	                  register __m256 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR =
	                                             _mm256_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR[0]);                
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_itlb_2l_mpi_ymm8r4(const float * __restrict pITLB_MISSES_WALK_COMPLETED,
	                                         const float * __restrict pINST_RETIRED_ANY) {
	                          
	                  register __m256 ITLB_MISSES_WALK_COMPLETED = 
	                                             _mm256_loadu_ps(&pITLB_MISSES_WALK_COMPLETED[0]);
	                  register __m256 INST_RETIRED_ANY           =
	                                             _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);               
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
	          
	            __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_itlb_2l_llp_mpi_ymm8r4(const float * __restrict pITLB_MISSES_WALK_COMPLETED_2M_4M,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 ITLB_MISSES_WALK_COMPLETED_2M_4M = 
	                                             _mm256_loadu_ps(&pITLB_MISSES_WALK_COMPLETED_2M_4M[0]);
	                  register __m256 INST_RETIRED_ANY                 =
	                                             _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_llmpi_ymm8r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 DTLB_LOAD_MISSES_WALK_COMPLETED = 
	                                           _mm256_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED[0]);
	                  register __m256 INST_RETIRED_ANY                =
	                                           _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	          
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_2mblp_mpi_ymm8r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                               const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M = 
	                                           _mm256_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M[0]);
	                  register __m256 INST_RETIRED_ANY                      =
	                                           _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	          
	            __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_llmpi_ymm8r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 DTLB_LOAD_MISSES_WALK_COMPLETED = 
	                                           _mm256_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED[0]);
	                  register __m256 INST_RETIRED_ANY                =
	                                           _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_2mllmpi_ymm8r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M = 
	                                           _mm256_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M[0]);
	                  register __m256 INST_RETIRED_ANY                      =
	                                           _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_2l_store_mpi_ymm8r4(const float * __restrict pDTLB_STORE_MISSES_WALK_COMPLETED,
	                                               const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m256 DTLB_STORE_MISSES_WALK_COMPLETED = 
	                                            _mm256_loadu_ps(&pDTLB_STORE_MISSES_WALK_COMPLETED[0]);
	                  register __m256 INST_RETIRED_ANY                 = 
	                                            _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
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
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_numa_reads_ald_ymm8r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE) {
	               
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL      = 
	                                         _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL = 
	                                         _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE     =
	                                         _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE=
	                                         _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE[0]);
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
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_numa_reads_ard_ymm8r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL) {
	               
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE      = 
	                                                     _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE =
	                                                     _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL       =
	                                                     _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL  = 
	                                                     _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL[0]);
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
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uncore_freq_ymm8r4(const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                         const float n_scks,
	                                         const float idt) { // shall be inverse value passed.
	                 
	                 register __m256 UNC_CHA_CLOCKTICKS = _mm256_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
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
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_upi_data_tx_bw_ymm8r4(const float * __restrict pUNC_UPI_TxL_FLITS_ALL_DATA,
	                                            const float idt) {
	             
	                  register __m256 UNC_UPI_TxL_FLITS_ALL_DATA = 
	                                         _mm256_loadu_ps(&pUNC_UPI_TxL_FLITS_ALL_DATA[0]);
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
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                            UNC_M_CAS_COUNT_RD,C640),
	                                                                 C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mem_bw_rd_ymm8r4(const float * __restrict pUNC_M_CAS_COUNT_RD,
	                                       const float idt) {
	                
	                  register __m256 UNC_M_CAS_COUNT_RD = 
	                                       _mm256_loadu_ps(&pUNC_M_CAS_COUNT_RD[0]);
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                            UNC_M_CAS_COUNT_RD,C640),
	                                                                 C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
/*
        "BriefDescription": "DDR memory write bandwidth (MB/sec)",
        "MetricExpr": "( UNC_M_CAS_COUNT.WR * 64 / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "memory_bandwidth_write",
        "ScaleUnit": "1MB/s"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mem_bw_wr_ymm8r4(const __m256 UNC_M_CAS_COUNT_WR,
	                                       const float idt) {
	                
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                            UNC_M_CAS_COUNT_WR,C640),
	                                                                 C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mem_bw_wr_ymm8r4(const float * __restrict pUNC_M_CAS_COUNT_WR,
	                                       const float idt) {
	                
	                  register __m256 UNC_M_CAS_COUNT_WR = 
	                                       _mm256_loadu_ps(&pUNC_M_CAS_COUNT_WR[0]);
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                            UNC_M_CAS_COUNT_WR,C640),
	                                                                 C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
/*
         "BriefDescription": "DDR memory bandwidth (MB/sec)",
        "MetricExpr": "(( UNC_M_CAS_COUNT.RD + UNC_M_CAS_COUNT.WR ) * 64 / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "memory_bandwidth_total",
        "ScaleUnit": "1MB/s"  
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mem_bw_total_ymm8r4(const __m256 UNC_M_CAS_COUNT_RD,
	                                          const __m256 UNC_M_CAS_COUNT_WR,
	                                          const float idt) {
	                  
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;  
	                  register __m256 t1;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0 =  _mm256_add_ps(UNC_M_CAS_COUNT_RD,
	                                      UNC_M_CAS_COUNT_WR);
	                  t1 =  _mm256_mul_ps(
	                                _mm256_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm256_mul_ps(t1,vidt);
	                  return (metric);                                                
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mem_bw_total_ymm8r4(const float * __restrict pUNC_M_CAS_COUNT_RD,
	                                          const float * __restrict pUNC_M_CAS_COUNT_WR,
	                                          const float idt) {
	                  
	                  register __m256 UNC_M_CAS_COUNT_RD = 
	                                       _mm256_loadu_ps(&pUNC_M_CAS_COUNT_RD[0]);
	                  register __m256 UNC_M_CAS_COUNT_WR =
	                                       _mm256_loadu_ps(&pUNC_M_CAS_COUNT_WR[0]);
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;  
	                  register __m256 t1;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0 =  _mm256_add_ps(UNC_M_CAS_COUNT_RD,
	                                      UNC_M_CAS_COUNT_WR);
	                  t1 =  _mm256_mul_ps(
	                                _mm256_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm256_mul_ps(t1,vidt);
	                  return (metric);                                                
	        }
	        
/*
         "BriefDescription": "Bandwidth of IO reads that are initiated by end device controllers that are requesting memory from the CPU.",
        "MetricExpr": "( UNC_CHA_TOR_INSERTS.IO_PCIRDCUR * 64 / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "io_bandwidth_read",
        "ScaleUnit": "1MB/s"
*/

           
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_io_bw_read_ymm8r4(const __m256 UNC_CHA_TOR_INSERTS_IO_PCIRDCUR,
	                                       const float idt) {
	                
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                         UNC_CHA_TOR_INSERTS_IO_PCIRDCUR,C640),
	                                                                        C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_io_bw_read_ymm8r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IO_PCIRDCUR,
	                                        const float idt) {
	                  
	                  register __m256 UNC_CHA_TOR_INSERTS_IO_PCIRDCUR = 
	                                    _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IO_PCIRDCUR[0]);
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                         UNC_CHA_TOR_INSERTS_IO_PCIRDCUR,C640),
	                                                                        C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
/* 
        "BriefDescription": "Bandwidth of IO writes that are initiated by end device controllers that are writing memory to the CPU.",
        "MetricExpr": "(( UNC_CHA_TOR_INSERTS.IO_ITOM + UNC_CHA_TOR_INSERTS.IO_ITOMCACHENEAR ) * 64 / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "io_bandwidth_write",
        "ScaleUnit": "1MB/s"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_io_bw_read_ymm8r4(  const __m256 UNC_CHA_TOR_INSERTS_IO_ITOM,
	                                          const __m256 UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR ,
	                                          const float idt) {
	                  
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;  
	                  register __m256 t1;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0 =  _mm256_add_ps(UNC_CHA_TOR_INSERTS_IO_ITOM,
	                                      UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR);
	                  t1 =  _mm256_mul_ps(
	                                _mm256_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm256_mul_ps(t1,vidt);
	                  return (metric);                                                
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_io_bw_read_ymm8r4(  const float * __restrict pUNC_CHA_TOR_INSERTS_IO_ITOM,
	                                          const float * __restrict pUNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR ,
	                                          const float idt) {
	                  
	                  register __m256 UNC_CHA_TOR_INSERTS_IO_ITOM          = 
	                                         _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IO_ITOM[0]);
	                  register __m256 UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR =
	                                         _mm256_loadu_ps(&pUNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR[0]);
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;  
	                  register __m256 t1;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0 =  _mm256_add_ps(UNC_CHA_TOR_INSERTS_IO_ITOM,
	                                      UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR);
	                  t1 =  _mm256_mul_ps(
	                                _mm256_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm256_mul_ps(t1,vidt);
	                  return (metric);                                                
	        }
	        
/*
        "BriefDescription": "Uops delivered from decoded instruction cache (decoded stream buffer or DSB) as a percent of total uops delivered to Instruction Decode Queue",
        "MetricExpr": "( IDQ.DSB_UOPS / ( IDQ.DSB_UOPS + IDQ.MITE_UOPS + IDQ.MS_UOPS + LSD.UOPS ) )",
        "MetricGroup": "",
        "MetricName": "percent_uops_delivered_from_decoded_icache",
        "ScaleUnit": "100%"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uops_ddic_ymm8r4(const __m256  IDQ_DSB_UOPS,
	                                       const __m256  IDQ_MITE_UOPS,
	                                       const __m256  IDQ_MS_UOPS,
	                                       const __m256  LSD_UOPS) {
	             
	                  register __m256 t0;
	               	  register __m256 metric;
	                  t0 = _mm256_add_ps(_mm256_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm256_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm256_div_ps(IDQ_DSB_UOPS,t0);
	                  return (metric);                             
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uops_ddic_ymm8r4(const float * __restrict  pIDQ_DSB_UOPS,
	                                       const float * __restrict  pIDQ_MITE_UOPS,
	                                       const float * __restrict  pIDQ_MS_UOPS,
	                                       const float * __restrict  pLSD_UOPS) {
	                  
	                  register __m256 IDQ_DSB_UOPS  = 
	                                     _mm256_loadu_ps(&pIDQ_DSB_UOPS[0]);
	                  register __m256 IDQ_MITE_UOPS =
	                                     _mm256_loadu_ps(&IDQ_MITE_UOPS[0]);
	                  register __m256 IDQ_MS_UOPS   =
	                                     _mm256_loadu_ps(&IDQ_MS_UOPS[0]);
	                  register __m256 LSD_UOPS      = 
	                                     _mm256_loadu_ps(&LSD_UOPS[0]);
	                  register __m256 t0;
	               	  register __m256 metric;
	                  t0 = _mm256_add_ps(_mm256_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm256_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm256_div_ps(IDQ_DSB_UOPS,t0);
	                  return (metric);                             
	         } 
	         
/*
         "BriefDescription": "Uops delivered from legacy decode pipeline (Micro-instruction Translation Engine or MITE) as a percent of total uops delivered to Instruction Decode Queue",
        "MetricExpr": "( IDQ.MITE_UOPS / ( IDQ.DSB_UOPS + IDQ.MITE_UOPS + IDQ.MS_UOPS + LSD.UOPS ) )",
        "MetricGroup": "",
        "MetricName": "percent_uops_delivered_from_legacy_decode_pipeline",
        "ScaleUnit": "100%"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uops_dldp_ymm8r4(const __m256  IDQ_DSB_UOPS,
	                                       const __m256  IDQ_MITE_UOPS,
	                                       const __m256  IDQ_MS_UOPS,
	                                       const __m256  LSD_UOPS) {
	             
	                  register __m256 t0;
	               	  register __m256 metric;
	                  t0 = _mm256_add_ps(_mm256_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm256_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm256_div_ps(IDQ_MITE_UOPS,t0);
	                  return (metric);                             
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uops_dldp_ymm8r4(const float * __restrict  pIDQ_DSB_UOPS,
	                                       const float * __restrict  pIDQ_MITE_UOPS,
	                                       const float * __restrict  pIDQ_MS_UOPS,
	                                       const float * __restrict  pLSD_UOPS) {
	             
	                  register __m256 IDQ_DSB_UOPS  = 
	                                     _mm256_loadu_ps(&pIDQ_DSB_UOPS[0]);
	                  register __m256 IDQ_MITE_UOPS =
	                                     _mm256_loadu_ps(&IDQ_MITE_UOPS[0]);
	                  register __m256 IDQ_MS_UOPS   =
	                                     _mm256_loadu_ps(&IDQ_MS_UOPS[0]);
	                  register __m256 LSD_UOPS      = 
	                                     _mm256_loadu_ps(&LSD_UOPS[0]);
	                  register __m256 t0;
	               	  register __m256 metric;
	                  t0 = _mm256_add_ps(_mm256_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm256_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm256_div_ps(IDQ_MITE_UOPS,t0);
	                  return (metric);                             
	         } 
	         
/*
        "BriefDescription": "Uops delivered from microcode sequencer (MS) as a percent of total uops delivered to Instruction Decode Queue",
        "MetricExpr": "( IDQ.MS_UOPS / ( IDQ.DSB_UOPS + IDQ.MITE_UOPS + IDQ.MS_UOPS + LSD.UOPS ) )",
        "MetricGroup": "",
        "MetricName": "percent_uops_delivered_from_microcode_sequencer",
        "ScaleUnit": "100%"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uops_dmseq_ymm8r4(const __m256  IDQ_DSB_UOPS,
	                                       const __m256  IDQ_MITE_UOPS,
	                                       const __m256  IDQ_MS_UOPS,
	                                       const __m256  LSD_UOPS) {
	             
	                  register __m256 t0;
	               	  register __m256 metric;
	                  t0 = _mm256_add_ps(_mm256_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm256_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm256_div_ps(IDQ_MS_UOPS,t0);
	                  return (metric);                             
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_uops_dmseq_ymm8r4(const float * __restrict  pIDQ_DSB_UOPS,
	                                       const  float * __restrict  pIDQ_MITE_UOPS,
	                                       const  float * __restrict  pIDQ_MS_UOPS,
	                                       const  float * __restrict  pLSD_UOPS) {
	             
	                  register __m256 IDQ_DSB_UOPS  = 
	                                     _mm256_loadu_ps(&pIDQ_DSB_UOPS[0]);
	                  register __m256 IDQ_MITE_UOPS =
	                                     _mm256_loadu_ps(&IDQ_MITE_UOPS[0]);
	                  register __m256 IDQ_MS_UOPS   =
	                                     _mm256_loadu_ps(&IDQ_MS_UOPS[0]);
	                  register __m256 LSD_UOPS      = 
	                                     _mm256_loadu_ps(&LSD_UOPS[0]);
	                  register __m256 t0;
	               	  register __m256 metric;
	                  t0 = _mm256_add_ps(_mm256_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm256_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm256_div_ps(IDQ_MS_UOPS,t0);
	                  return (metric);                             
	         } 
	         
/*
        "BriefDescription": "Bandwidth (MB/sec) of read requests that miss the last level cache (LLC) and go to local memory.",
        "MetricExpr": "( UNC_CHA_REQUESTS.READS_LOCAL * 64 / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "llc_miss_local_memory_bandwidth_read",
        "ScaleUnit": "1MB/s"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_mlmbw_rd_ymm8r4(const __m256 UNC_CHA_REQUESTS_READS_LOCAL,
	                                          const float idt) {
	                
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_LOCAL,C640),
	                                                                        C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_mlmbw_rd_ymm8r4(const float * __restrict pUNC_CHA_REQUESTS_READS_LOCAL,
	                                          const float idt) {
	                
	                  register __m256 UNC_CHA_REQUESTS_READS_LOCAL = 
	                                         _mm256_loadu_ps(&pUNC_CHA_REQUESTS_READS_LOCAL[0]);
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_LOCAL,C640),
	                                                                        C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
/*
        "BriefDescription": "Bandwidth (MB/sec) of read requests that miss the last level cache (LLC) and go to remote memory.",
        "MetricExpr": "( UNC_CHA_REQUESTS.READS_REMOTE * 64 / 1000000) / duration_time",
        "MetricGroup": "",
        "MetricName": "llc_miss_remote_memory_bandwidth_read",
        "ScaleUnit": "1MB/s"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_mrmbw_rd_ymm8r4(const __m256 UNC_CHA_REQUESTS_READS_REMOTE,
	                                          const float idt) {
	                
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_REMOTE,C640),
	                                                                        C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	            __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_llc_mrmbw_rd_ymm8r4(const float * __restrict pUNC_CHA_REQUESTS_READS_REMOTE,
	                                          const float idt) {
	                
	                  register __m256 UNC_CHA_REQUESTS_READS_REMOTE = 
	                                         _mm256_loadu_ps(&pUNC_CHA_REQUESTS_READS_REMOTE[0]);
	                  const __m256 C640    = _mm256_mul_ps(64.0f);
	                  const __m256 C000001 = _mm256_set1_ps(1.0e-06f);              
	                  register __m256 vidt;
	                  register __m256 t0;
	                  register __m256 metric;
	                  vidt = _mm256_set1_ps(idt);
	                  t0   = _mm256_mul_ps(_mm256_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_REMOTE,C640),
	                                                                        C000001);
	                  metric = _mm256_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
/*
     "BriefDescription": "This metric represents fraction of slots the CPU was stalled due to Frontend latency issues.  For example; instruction-cache misses; iTLB misses or fetch stalls after a branch misprediction are categorized under Frontend Latency. In such cases; the Frontend eventually delivers no uops for some period.",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           float spr_fetch_latency_r4(const float PERF_METRICS_FETCH_LATENCY,
	                                  const float PERF_METRICS_FRONTEND_BOUND,
	                                  const float PERF_METRICS_BAD_SPECULATION,
	                                  const float PERF_METRICS_RETIRING,
	                                  const float PERF_METRICS_BACKEND_BOUND,
	                                  const float INT_MISC_UOP_DROPPING,
	                                  const float TOPDOWN_SLOTS_perf_metrics) {
	                 
	                  constexpr C100 = 100.0f;
	                  float t0,t1,t2;
	                  float metric = 0.0f;
	                  t0 =  INT_MISC_UOP_DROPPING/TOPDOWN_SLOTS_perf_metrics;
	                  t1 =  PERF_METRICS_FRONTEND_BOUND+
	                        PERF_METRICS_BAD_SPECULATION+
	                        PERF_METRICS_RETIRING+
	                        PERF_METRICS_BACKEND_BOUND;
	                  t2 = t1-t0;
	                  metric = C100*t2;
	                  return (metric);                     
	         }
	         
/*
     "BriefDescription": "This category represents fraction of slots where the processor's Frontend undersupplies its Backend. Frontend denotes the first part of the processor core responsible to fetch operations that are executed later on by the Backend part. Within the Frontend; a branch predictor predicts the next address to fetch; cache-lines are fetched from the memory subsystem; parsed into instructions; and lastly decoded into micro-operations (uops). Ideally the Frontend can issue Machine_Width uops every cycle to the Backend. Frontend Bound denotes unutilized issue-slots when there is no Backend stall; i.e. bubbles where Frontend delivered no uops while Backend could have accepted them. For example; stalls due to instruction-cache misses would be categorized under Frontend Bound.",
      "UnitOfMeasure": "percent",
*/	

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           float spr_frontend_bound_r4(
	                                   const float PERF_METRICS_FRONTEND_BOUND,
	                                   const float PERF_METRICS_BAD_SPECULATION,
	                                   const float PERF_METRICS_RETIRING,
	                                   const float PERF_METRICS_BACKEND_BOUND,
	                                   const float INT_MISC_UOP_DROPPING,
	                                   const float TOPDOWN_SLOTS_perf_metrics) {
	                 
	                  constexpr C100 = 100.0f;
	                  float t0,t1,t2;
	                  float metric = 0.0f;
	                  t0 =  INT_MISC_UOP_DROPPING/TOPDOWN_SLOTS_perf_metrics;
	                  t1 =  PERF_METRICS_FRONTEND_BOUND+
	                        PERF_METRICS_BAD_SPECULATION+
	                        PERF_METRICS_RETIRING+
	                        PERF_METRICS_BACKEND_BOUND;
	                  t2 = t1-t0;
	                  metric = C100*t2;
	                  return (metric);                     
	         }   
	         
/*
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to instruction cache misses.",
      "UnitOfMeasure": "percent",
*/     

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_icache_misses_ymm8r4(const __m256 ICACHE_DATA_STALLS,
	                                           const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(ICACHE_DATA_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                                 
	         }   
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_icache_misses_ymm8r4(const float * __restrict pICACHE_DATA_STALLS,
	                                           const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m256 ICACHE_DATA_STALLS = 
	                                        _mm256_loadu_ps(&pICACHE_DATA_STALLS[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(ICACHE_DATA_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                                 
	         }   
	         
/*
      "MetricName": "ITLB_Misses",
      "LegacyName": "metric_TMA_....ITLB_Misses(%)",
      "ParentCategory": "Fetch_Latency",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to Instruction TLB (ITLB) misses.",
      "UnitOfMeasure": "percent",
*/   

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_itlb_misses_ymm8r4(const __m256 ICACHE_TAG_STALLS,
	                                         const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(ICACHE_TAG_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                                 
	         }   
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_itlb_misses_ymm8r4(const float * __restrict pICACHE_TAG_STALLS,
	                                         const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m256 ICACHE_TAG_STALLS = 
	                                        _mm256_loadu_ps(&pICACHE_TAG_STALLS[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(ICACHE_TAG_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                                 
	         }   
	         
/*
      "MetricName": "Branch_Resteers",
      "LegacyName": "metric_TMA_....Branch_Resteers(%)",
      "ParentCategory": "Fetch_Latency",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to Branch Resteers. Branch Resteers estimates the Frontend delay in fetching operations from corrected path; following all sorts of miss-predicted branches. For example; branchy code with lots of miss-predictions might get categorized under Branch Resteers. Note the value of this node may overlap with its siblings.",
      "UnitOfMeasure": "percent",
*/    

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_branch_resteers_ymm8r4(const __m256 INT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const __m256 CPU_CLK_UNHALTED_THREAD,
	                                             const __m256 INT_MISC_UNKNOWN_BRANCH_CYCLES) {
	                  
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 metric;  
	                  t0  = _mm256_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1  = _mm256_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,
	                                     _mm256_add_ps(t0,t1));
	                  return (metric);                                
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_branch_resteers_ymm8r4(const float * __restrict pINT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                             const float * __restrict pINT_MISC_UNKNOWN_BRANCH_CYCLES) {
	                  
	                  register __m256 INT_MISC_CLEAR_RESTEER_CYCLES = 
	                                          _mm256_loadu_ps(&pINT_MISC_CLEAR_RESTEER_CYCLES[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  register __m256 INT_MISC_UNKNOWN_BRANCH_CYCLES=
	                                          _mm256_loadu_ps(&pINT_MISC_UNKNOWN_BRANCH_CYCLES[0]);
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 metric;  
	                  t0  = _mm256_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1  = _mm256_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,
	                                     _mm256_add_ps(t0,t1));
	                  return (metric);                                
	         }
	         
/*
        "MetricName": "Mispredicts_Resteers",
      "LegacyName": "metric_TMA_......Mispredicts_Resteers(%)",
      "ParentCategory": "Branch_Resteers",
      "Level": 4,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to Branch Resteers as a result of Branch Misprediction at execution stage. ",
      "UnitOfMeasure": "percent",
*/

#include <algorithm>

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           float spr_mispredict_resteers_r4(const float PERF_METRICS_BRANCH_MISPREDICTS,
	                                            const float PERF_METRICS_FRONTEND_BOUND,
	                                            const float PERF_METRICS_BAD_SPECULATION,
	                                            const float PERF_METRICS_RETIRING,
	                                            const float PERF_METRICS_BACKEND_BOUND,
	                                            const float INT_MISC_UOP_DROPPING,
	                                            const float TOPDOWN_SLOTS_perf_metrics,
	                                            const float INT_MISC_CLEAR_RESTEER_CYCLES,
	                                            const float CPU_CLK_UNHALTED_THREAD) {
	                 
	                 constexpr float C100 = 100.0f;
	                 float t0,t1,t2,t3,t4,t5,t6,t7;
	                 float metric;
	                 t0 = PERF_METRICS_FRONTEND_BOUND+
	                      PERF_METRICS_BAD_SPECULATION+
	                      PERF_METRICS_RETIRING+
	                      PERF_METRICS_BACKEND_BOUND;
	                 t1 = INT_MISC_UOP_DROPPING/ 
	                      TOPDOWN_SLOTS_perf_metrics;
	                 t2 = PERF_METRICS_BACKEND_BOUND/t0;
	                 t3 = PERF_METRICS_RETIRING/t0;
	                 t4 = PERF_METRICS_BRANCH_MISPREDICTS/t0;
	                 t5 = PERF_METRICS_FRONTEND_BOUND/t0;
	                 t6 = t5-t1+t2+t3;
	                 t7 = 1.0f-t6;
	                 t2 = std::max(t7,0.0f);
	                 t3 = INT_MISC_CLEAR_RESTEER_CYCLES/
	                      CPU_CLK_UNHALTED_THREAD;
	                 metric = C100*(t4/t2)*t3;
	                 return (metric);                    
	         }
	         
/*
      "MetricName": "Clears_Resteers",
      "LegacyName": "metric_TMA_......Clears_Resteers(%)",
      "ParentCategory": "Branch_Resteers",
      "Level": 4,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to Branch Resteers as a result of Machine Clears. ",
      "UnitOfMeasure": "percent",
*/
	         
 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           float spr_clears_resteers_r4(    const float PERF_METRICS_BRANCH_MISPREDICTS,
	                                            const float PERF_METRICS_FRONTEND_BOUND,
	                                            const float PERF_METRICS_BAD_SPECULATION,
	                                            const float PERF_METRICS_RETIRING,
	                                            const float PERF_METRICS_BACKEND_BOUND,
	                                            const float INT_MISC_UOP_DROPPING,
	                                            const float TOPDOWN_SLOTS_perf_metrics,
	                                            const float INT_MISC_CLEAR_RESTEER_CYCLES,
	                                            const float CPU_CLK_UNHALTED_THREAD) {
	                 
	                 constexpr float C100 = 100.0f;
	                 float t0,t1,t2,t3,t4,t5,t6,t7;
	                 float metric;
	                 t0 = PERF_METRICS_FRONTEND_BOUND+
	                      PERF_METRICS_BAD_SPECULATION+
	                      PERF_METRICS_RETIRING+
	                      PERF_METRICS_BACKEND_BOUND;
	                 t1 = INT_MISC_UOP_DROPPING/ 
	                      TOPDOWN_SLOTS_perf_metrics;
	                 t2 = PERF_METRICS_BACKEND_BOUND/t0;
	                 t3 = PERF_METRICS_RETIRING/t0;
	                 t4 = PERF_METRICS_BRANCH_MISPREDICTS/t0;
	                 t5 = PERF_METRICS_FRONTEND_BOUND/t0;
	                 t6 = t5-t1+t2+t3;
	                 t7 = 1.0f-t6;
	                 t2 = std::max(t7,0.0f);
	                 t3 = INT_MISC_CLEAR_RESTEER_CYCLES/
	                      CPU_CLK_UNHALTED_THREAD;
	                 metric = C100*(1.0f-(t4/t2))*t3;
	                 return (metric);                    
	         }
	         
/*
      "MetricName": "Unknown_Branches",
      "LegacyName": "metric_TMA_......Unknown_Branches(%)",
      "ParentCategory": "Branch_Resteers",
      "Level": 4,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to new branch address clears. These are fetched branches the Branch Prediction Unit was unable to recognize (First fetch or hitting BPU capacity limit).",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_unknown_branches_ymm8r4(const __m256 INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                              const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);
	                                                      
	         }   
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_unknown_branches_ymm8r4(const float * __restrict pINT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                              const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m256 INT_MISC_UNKNOWN_BRANCH_CYCLES = 
	                                       _mm256_loadu_ps(&pINT_MISC_UNKNOWN_BRANCH_CYCLES[0]);
	                   register __m256 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);
	                                                      
	         }   
	         
/*
     "MetricName": "DSB_Switches",
      "LegacyName": "metric_TMA_....DSB_Switches(%)",
      "ParentCategory": "Fetch_Latency",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to switches from DSB to MITE pipelines. The DSB (decoded i-cache) is a Uop Cache where the front-end directly delivers Uops (micro operations) avoiding heavy x86 decoding. The DSB pipeline has shorter latency and delivered higher bandwidth than the MITE (legacy instruction decode pipeline). Switching between the two pipelines can cause penalties hence this metric measures the exposed penalty.",
      "UnitOfMeasure": "percent",
*/  

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dsb_switches_ymm8r4(const __m256 DSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                          const __m256 CPU_CLK_UNHALTED_THREAD) {
	                 
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;  
	                  t0 = _mm256_div_ps(DSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                              
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dsb_switches_ymm8r4(const float * __restrict pDSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                          const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                  register __m256 DSB2MITE_SWITCHES_PENALTY_CYCLES = 
	                                          _mm256_loadu_ps(&pDSB2MITE_SWITCHES_PENALTY_CYCLES[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;  
	                  t0 = _mm256_div_ps(DSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                              
	         }
	         
/*
        "MetricName": "LCP",
      "LegacyName": "metric_TMA_....LCP(%)",
      "ParentCategory": "Fetch_Latency",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of cycles CPU was stalled due to Length Changing Prefixes (LCPs). Using proper compiler flags or Intel Compiler by default will certainly avoid this. #Link: Optimization Guide about LCP BKMs.",
      "UnitOfMeasure": "percent",
*/

              	   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_lcp_ymm8r4(const __m256 DECODE_LCP,
	                                 const __m256 CPU_CLK_UNHALTED_THREAD) {
	                 
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;  
	                  t0 = _mm256_div_ps(DECODE_LCP,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                            
	         }    
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_lcp_ymm8r4(const float * __restrict pDECODE_LCP,
	                                 const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                  register __m256 DECODE_LCP              =
	                                  _mm256_loadu_ps(&pDECODE_LCP[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD =
	                                          _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C1000 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;  
	                  t0 = _mm256_div_ps(DECODE_LCP,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C1000,t0);
	                  return (metric);                            
	         }    
	         
/*
      "MetricName": "Info_Memory_L2MPKI",
      "LegacyName": "metric_TMA_Info_Memory_L2MPKI",
      "Level": 1,
      "BriefDescription": "L2 cache true misses per kilo instruction for retired demand loads",
      "UnitOfMeasure": ""
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_info_mem_l2mpki_ymm8r4(const __m256 MEM_LOAD_RETIRED_L2_MISS,
	                                             const __m256 INST_RETIRED_ANY) {
	                                             
	                  const __m256 C10000 = _mm256_set1_ps(1000.0f);
	                  register __m256 t0;
	                  register __m256 metric; 
	                  t0 = _mm256_div_ps(MEM_LOAD_RETIRED_L2_MISS,
	                                     INST_RETIRED_ANY);
	                  metric = _mm256_mul_ps(C10000,t0);
	                  return (metric);                                                       
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_info_mem_l2mpki_ymm8r4(const float * __restrict pMEM_LOAD_RETIRED_L2_MISS,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                                  
	                  register __m256 MEM_LOAD_RETIRED_L2_MISS = 
	                                          _mm256_loadu_ps(&pMEM_LOAD_RETIRED_L2_MISS[0]);
	                  register __m256 INST_RETIRED_ANY         =
	                                          _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);          
	                  const __m256 C10000 = _mm256_set1_ps(1000.0f);
	                  register __m256 t0;
	                  register __m256 metric; 
	                  t0 = _mm256_div_ps(MEM_LOAD_RETIRED_L2_MISS,
	                                     INST_RETIRED_ANY);
	                  metric = _mm256_mul_ps(C10000,t0);
	                  return (metric);                                                       
	         }
	         
/*
      "MetricName": "Info_Bad_Spec_IpMispredict",
      "LegacyName": "metric_TMA_Info_Bad_Spec_IpMispredict",
      "Level": 1,
      "BriefDescription": "Number of Instructions per non-speculative Branch Misprediction (JEClear) (lower number means higher occurrence rate)",
      "UnitOfMeasure": ""
*/	

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_bad_spec_ipm_ymm8r4(const __m256 INST_RETIRED_ANY,
	                                          const __m256 BR_MISP_RETIRED_ALL_BRANCHES) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(INST_RETIRED_ANY,
	                                     BR_MISP_RETIRED_ALL_BRANCHES);
	                  return (metric);                          
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_bad_spec_ipm_ymm8r4(const float * __restrict pINST_RETIRED_ANY,
	                                          const float * __restrict pBR_MISP_RETIRED_ALL_BRANCHES) {
	                  
	                  register __m256 INST_RETIRED_ANY = 
	                                      _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m256  BR_MISP_RETIRED_ALL_BRANCHES = 
	                                      _mm256_loadu_ps(&pBR_MISP_RETIRED_ALL_BRANCHES[0]);
	                  register __m256 metric;
	                  metric = _mm256_div_ps(INST_RETIRED_ANY,
	                                     BR_MISP_RETIRED_ALL_BRANCHES);
	                  return (metric);                          
	         }   
	         
/*
       "MetricName": "Info_Inst_Mix_IpTB",
      "LegacyName": "metric_TMA_Info_Inst_Mix_IpTB",
      "Level": 1,
      "BriefDescription": "Instruction per taken branch",
      "UnitOfMeasure": "",
*/    

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_inst_mix_iptb_ymm8r4(const __m256 INST_RETIRED_ANY,
	                                           const __m256 BR_INST_RETIRED_NEAR_TAKEN) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(INST_RETIRED_ANY,
	                                         BR_INST_RETIRED_NEAR_TAKEN);
	                  return (metric);                          
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_inst_mix_iptb_ymm8r4(const float * __restrict pINST_RETIRED_ANY,
	                                           const float * __restrict pBR_INST_RETIRED_NEAR_TAKEN) {
	                  
	                  register __m256 INST_RETIRED_ANY = 
	                                      _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m256  BR_INST_RETIRED_NEAR_TAKEN = 
	                                      _mm256_loadu_ps(&pBR_INST_RETIRED_NEAR_TAKEN[0]);
	                  register __m256 metric;
	                  metric = _mm256_div_ps(INST_RETIRED_ANY,
	                                         BR_INST_RETIRED_NEAR_TAKEN);
	                  return (metric);                          
	         }   
	         
/*
      "MetricName": "Info_Core_CoreIPC",
      "LegacyName": "metric_TMA_Info_Core_CoreIPC",
      "Level": 1,
      "BriefDescription": "Instructions Per Cycle across hyper-threads (per physical core)",
      "UnitOfMeasure": "",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_core_coreipc_ymm8r4( const __m256 INST_RETIRED_ANY,
	                                           const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                  register __m256 metric;
	                  metric = _mm256_div_ps(INST_RETIRED_ANY,
	                                         CPU_CLK_UNHALTED_DISTRIBUTED);
	                  return (metric);                          
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_core_coreipc_ymm8r4( const float * __restrict pINST_RETIRED_ANY,
	                                           const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                  register __m256 INST_RETIRED_ANY = 
	                                      _mm256_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m256  CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                      _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                  register __m256 metric;
	                  metric = _mm256_div_ps(INST_RETIRED_ANY,
	                                         CPU_CLK_UNHALTED_DISTRIBUTED);
	                  return (metric);                          
	         }   
	         
/*
       "MetricName": "CISC",
      "LegacyName": "metric_TMA_......CISC(%)",
      "ParentCategory": "Microcode_Sequencer",
      "Level": 4,
      "BriefDescription": "This metric estimates fraction of cycles the CPU retired uops originated from CISC (complex instruction set computer) instruction. A CISC instruction has multiple  uops that are required to perform the instruction's functionality as in the case of read-modify-write as an example. Since these instructions require multiple uops they may or may not imply sub-optimal use of machine resources.",
      "UnitOfMeasure": "percent"
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline 
                   float spr_cisc_frac_r4(const float UOPS_RETIRED_MS,
                                          const float TOPDOWN_SLOTS_perf_metrics,
                                          const float ASSISTS_ANY_u0x1B) {
                         
                         constexpr float C100 = 100.0f;
                         float ab,cb;
                         float t0,t1;
                         float metric;
                         ab = UOPS_RETIRED_MS/TOPDOWN_SLOTS_perf_metrics;
                         cb = ASSISTS_ANY_u0x1B/TOPDOWN_SLOTS_perf_metrics;  
                         t0 = std::min(C100*cb,1.0f);   
                         t1 = std::max(0.0f,ab-t0);
                         metric = C100*t1;
                         return (metric);                   
                 } 
                 
/*
      "MetricName": "AVX_Assists",
      "LegacyName": "metric_TMA_........AVX_Assists(%)",
      "ParentCategory": "Assists",
      "Level": 5,
      "BriefDescription": "This metric estimates fraction of slots the CPU retired uops as a result of handing SSE to AVX* or AVX* to SSE transition Assists.",
      "UnitOfMeasure": "percent",
*/ 

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_avx_assists_ymm8r4(const __m256 ASSISTS_SSE_AVX_MIX,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C63  = _mm256_set1_ps(63.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_mul_ps(C63,
	                             _mm256_div_ps(ASSISTS_SSE_AVX_MIX,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_avx_assists_ymm8r4(const float * __restrict pASSISTS_SSE_AVX_MIX,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m256 ASSISTS_SSE_AVX_MIX = 
	                                         _mm256_loadu_ps(&pASSISTS_SSE_AVX_MIX[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C63  = _mm256_set1_ps(63.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_mul_ps(C63,
	                             _mm256_div_ps(ASSISTS_SSE_AVX_MIX,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
/*
      "MetricName": "FP_Assists",
      "LegacyName": "metric_TMA_........FP_Assists(%)",
      "ParentCategory": "Assists",
      "Level": 5,
      "BriefDescription": "This metric roughly estimates fraction of slots the CPU retired uops as a result of handing Floating Point (FP) Assists. FP Assist may apply when working with very small floating point values (so-called Denormals).",
      "UnitOfMeasure": "percent",         
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_fp_assists_ymm8r4(const __m256 ASSISTS_FP,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C30  = _mm256_set1_ps(30.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_mul_ps(C30,
	                             _mm256_div_ps(ASSISTS_FP,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_fp_assists_ymm8r4(const float * __restrict pASSISTS_FP,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m256 ASSISTS_FP = 
	                                         _mm256_loadu_ps(&pASSISTS_FP[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C30  = _mm256_set1_ps(30.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_mul_ps(C30,
	                             _mm256_div_ps(ASSISTS_FP,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
/*
      "MetricName": "Page_Faults",
      "LegacyName": "metric_TMA_........Page_Faults(%)",
      "ParentCategory": "Assists",
      "Level": 5,
      "BriefDescription": "This metric roughly estimates fraction of slots the CPU retired uops as a result of handing Page Faults. A Page Fault may apply on first application access to a memory page. Note operating system handling of page faults accounts for the majority of its cost.",
      "UnitOfMeasure": "percent",
*/

                  __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_page_faults_ymm8r4(const __m256 ASSISTS_PAGE_FAULTS,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C99  = _mm256_set1_ps(99.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_mul_ps(C99,
	                             _mm256_div_ps(ASSISTS_PAGE_FAULTS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_page_faults_ymm8r4(const float * __restrict pASSISTS_PAGE_FAULTS,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m256 ASSISTS_PAGE_FAULTS = 
	                                         _mm256_loadu_ps(&pASSISTS_PAGE_FAULTS[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C99  = _mm256_set1_ps(99.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_mul_ps(C99,
	                             _mm256_div_ps(ASSISTS_PAGE_FAULTS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
/*
        "MetricName": "Assists",
      "LegacyName": "metric_TMA_......Assists(%)",
      "ParentCategory": "Microcode_Sequencer",
      "Level": 4,
      "BriefDescription": "This metric estimates fraction of cycles the CPU retired uops delivered by the Microcode_Sequencer as a result of Assists. Assists are long sequences of uops that are required in certain corner-cases for operations that cannot be handled natively by the execution pipeline. For example; when working with very small floating point values (so-called Denormals); the FP units are not set up to perform these operations natively. Instead; a sequence of instructions to perform the computation on the Denormals is injected into the pipeline. Since these microcode sequences might be hundreds of uops long; Assists can be extremely deleterious to performance and they can be avoided in many cases. #Link: article for compiler flags of DAZ and FTZ",
      "UnitOfMeasure": "percent"  
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline 
                   float spr_assists_r4(const float ASSISTS_ANY_u0x1B,
                                        const float TOPDOWN_SLOTS_perf_metrics) {
                                        
                         constexpr float C100 = 100.0f;
                         float ab,t0;
                         float metric;
                         ab = ASSISTS_ANY_u0x1B/TOPDOWN_SLOTS_perf_metrics;
                         t0 = std::min(C100*ab,1.0f);
                         metric = C100*t0;
                         return (metric);                    
                 }
                 
/*
        "MetricName": "Microcode_Sequencer",
      "LegacyName": "metric_TMA_....Microcode_Sequencer(%)",
      "ParentCategory": "Heavy_Operations",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of slots the CPU was retiring uops fetched by the Microcode Sequencer (MS) unit.  The MS is used for CISC instructions not supported by the default decoders (like repeat move strings; or CPUID); or by microcode assists used to address some operation modes (like in Floating Point assists). These cases can often be avoided.",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_micrcode_seq_ymm8r4(const __m256 UOPS_RETIRED_MS,
	                                          const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_div_ps(UOPS_RETIRED_MS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_micrcode_seq_ymm8r4(const float * __restrict pUOPS_RETIRED_MS,
	                                          const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m256 UOPS_RETIRED_MS = 
	                                      _mm256_loadu_ps(&pUOPS_RETIRED_MS[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 vx;
	                  register __m256 t0;
	                  register __metric;
	                  vx = _mm256_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm256_div_ps(UOPS_RETIRED_MS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
/*
     "MetricName": "Few_Uops_Instructions",
      "LegacyName": "metric_TMA_....Few_Uops_Instructions(%)",
      "ParentCategory": "Heavy_Operations",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of slots where the CPU was retiring instructions that that are decoder into two or up to ([SNB+] four; [ADL+] five) uops. This highly-correlates with the number of uops in such instructions.",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline 
                   float  spr_few_uops_instr_r4(const float PERF_METRICS_HEAVY_OPERATIONS,
                                                const float PERF_METRICS_FRONTEND_BOUND,
                                                const float PERF_METRICS_BAD_SPECULATION,
                                                const float PERF_METRICS_RETIRING,
                                                const float PERF_METRICS_BACKEND_BOUND,
                                                const float UOPS_RETIRED_MS,
                                                const float TOPDOWN_SLOTS_perf_metrics) {
                       
                          constexpr C100 = 100.0f;
                          float bcde,fg;
                          float t0,t1,metric;
                          bcde = PERF_METRICS_FRONTEND_BOUND+
                                 PERF_METRICS_BAD_SPECULATION+
                                 PERF_METRICS_RETIRING+
                                 PERF_METRICS_BACKEND_BOUND;
                         fg    = UOPS_RETIRED_MS/ 
                                 TOPDOWN_SLOTS_perf_metrics;
                         t0    = PERF_METRICS_HEAVY_OPERATIONS/bcde;   
                         t1    = std::max(0.0f,t0-fg);
                         metric = C100*t1;
                         return (metric);                         
                 }
                 
/*
    "MetricName": "Heavy_Operations",
      "LegacyName": "metric_TMA_..Heavy_Operations(%)",
      "ParentCategory": "Retiring",
      "Level": 2,
      "BriefDescription": "This metric represents fraction of slots where the CPU was retiring heavy-weight operations -- instructions that require two or more uops. This highly-correlates with the uop length of these instructions/flows.",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline 
                   float spr_heavy_ops_r4(const float PERF_METRICS_HEAVY_OPERATIONS,
                                          const float PERF_METRICS_FRONTEND_BOUND,
                                          const float PERF_METRICS_BAD_SPECULATION,
                                          const float PERF_METRICS_RETIRING,
                                          const float PERF_METRICS_BACKEND_BOUND) {
                     
                         constexpr C100 = 100.0f;
                         float bcde,t0;
                         float metric;
                         bcde = PERF_METRICS_FRONTEND_BOUND+
                                PERF_METRICS_BAD_SPECULATION+
                                PERF_METRICS_RETIRING+
                                PERF_METRICS_BACKEND_BOUND;
                         t0   = PERF_METRICS_HEAVY_OPERATIONS/bcde;
                         metric = C100*t0;
                         return (metric);                     
                  }
                  
/*
       "MetricName": "MS_Switches",
      "LegacyName": "metric_TMA_....MS_Switches(%)",
      "ParentCategory": "Fetch_Latency",
      "Level": 3,
      "BriefDescription": "This metric estimates the fraction of cycles when the CPU was stalled due to switches of uop delivery to the Microcode Sequencer (MS). Commonly used instructions are optimized for delivery by the DSB (decoded i-cache) or MITE (legacy instruction decode) pipelines. Certain operations cannot be handled natively by the execution pipeline; and must be performed by microcode (small programs injected into the execution stream). Switching to the MS too often can negatively impact performance. The MS is designated to deliver long uop flows required by CISC instructions like CPUID; or uncommon conditions like Floating Point Assists when dealing with Denormals.",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline 
                   float spr_ms_switches_r4(const float UOPS_RETIRED_MS_c1_e1,
                                            const float PERF_METRICS_RETIRING,
                                            const float PERF_METRICS_FRONTEND_BOUND,
                                            const float PERF_METRICS_BAD_SPECULATION,
                                            const float PERF_METRICS_BACKEND_BOUND,
                                            const float TOPDOWN_SLOTS_perf_metrics,
                                            const float UOPS_ISSUED_ANY,
                                            const float CPU_CLK_UNHALTED_THREAD) {
                                            
                         constexpr C100 = 100.0f;
                         constexpr C3   = 3.0f;
                         float cdbe,t0;
                         float t1,t2,
                         float metric;
                         cdbe = PERF_METRICS_FRONTEND_BOUND+
                                PERF_METRICS_BAD_SPECULATION+
                                PERF_METRICS_RETIRING+
                                PERF_METRICS_BACKEND_BOUND;
                         t0   = PERF_METRICS_RETIRING/cdbe*
                                TOPDOWN_SLOTS_perf_metrics;
                         t1   = UOPS_RETIRED_MS_c1_e1/t0/
                                UOPS_ISSUED_ANY;
                         t2   = (C3*t1)/CPU_CLK_UNHALTED_THREAD;
                         metric = C100*std::min(t2,1.0f);
                         return (metric);
                                                       
                }
                
/*
        "MetricName": "Branch_Resteers",
      "LegacyName": "metric_TMA_....Branch_Resteers(%)",
      "ParentCategory": "Fetch_Latency",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to Branch Resteers. Branch Resteers estimates the Frontend delay in fetching operations from corrected path; following all sorts of miss-predicted branches. For example; branchy code with lots of miss-predictions might get categorized under Branch Resteers. Note the value of this node may overlap with its siblings.",
      "UnitOfMeasure": "percent"  
*/  

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_branch_resteers_ymm8r4(const __m256 INT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const __m256 CPU_CLK_UNHALTED_THREAD,
	                                             const __m256 INT_MISC_UNKNOWN_BRANCH_CYCLES) {
	           
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0,t1;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,
	                                     _mm256_add_ps(t0,t1));
	                  return (metric);                                  
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_branch_resteers_ymm8r4(const float * __restrict pINT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                             const float * __restrict pINT_MISC_UNKNOWN_BRANCH_CYCLES) {
	            
	                  register __m256 INT_MISC_CLEAR_RESTEER_CYCLES = 
	                                          _mm256_loadu_ps(&pINT_MISC_CLEAR_RESTEER_CYCLES[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm256_loadu_ps(&pINT_MISC_CLEAR_RESTEER_CYCLES[0]);
	                  register __m256 INT_MISC_UNKNOWN_BRANCH_CYCLES=
	                                          _mm256_loadu_ps(&pINT_MISC_UNKNOWN_BRANCH_CYCLES[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0,t1;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,
	                                     _mm256_add_ps(t0,t1));
	                  return (metric);                                  
	          }     
	          
/*
     "MetricName": "MITE",
      "LegacyName": "metric_TMA_....MITE(%)",
      "ParentCategory": "Fetch_Bandwidth",
      "Level": 3,
      "BriefDescription": "This metric represents Core fraction of cycles in which CPU was likely limited due to the MITE pipeline (the legacy decode pipeline). This pipeline is used for code that was not pre-cached in the DSB or LSD. For example; inefficiencies due to asymmetric decoders; use of long immediate or LCP can manifest as MITE fetch bandwidth bottleneck.",
      "UnitOfMeasure": "percent",
*/	

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mite_ymm8r4(const __m256 IDQ_MITE_CYCLES_ANY,
	                                  const __m256 IDQ_MITE_CYCLES_OK,
	                                  const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                    const __m256 C100 = _mm256_set1_ps(100.0f);
	                    const __m256 C05  = _mm256_set1_ps(0.5f);
	                    register __m256 t0,t1;
	                    register __m256 metric;
	                    t0 = _mm256_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm256_sub_ps(IDQ_MITE_CYCLES_ANY,
	                                       IDQ_MITE_CYCLES_OK);
	                    metric = _mm256_mul_ps(C100,
	                                   _mm256_div_ps(t1,t0));
	                    return (metric);                        
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_mite_ymm8r4(const float * __restrict pIDQ_MITE_CYCLES_ANY,
	                                  const float * __restrict pIDQ_MITE_CYCLES_OK,
	                                  const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                    register __m256 IDQ_MITE_CYCLES_ANY  =
	                                            _mm256_loadu_ps(&pIDQ_MITE_CYCLES_ANY[0]);
	                    register __m256 IDQ_MITE_CYCLES_OK   =
	                                            _mm256_loadu_ps(&pIDQ_MITE_CYCLES_OK[0]);
	                    register __m256 CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                            _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                    const __m256 C100 = _mm256_set1_ps(100.0f);
	                    const __m256 C05  = _mm256_set1_ps(0.5f);
	                    register __m256 t0,t1;
	                    register __m256 metric;
	                    t0 = _mm256_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm256_sub_ps(IDQ_MITE_CYCLES_ANY,
	                                       IDQ_MITE_CYCLES_OK);
	                    metric = _mm256_mul_ps(C100,
	                                   _mm256_div_ps(t1,t0));
	                    return (metric);                        
	         }   
	         
/*
    "MetricName": "Decoder0_Alone",
      "LegacyName": "metric_TMA_......Decoder0_Alone(%)",
      "ParentCategory": "MITE",
      "Level": 4,
      "BriefDescription": "This metric represents fraction of cycles where decoder-0 was the only active decoder",
      "UnitOfMeasure": "percent",
*/   

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_decoder0_only_ymm8r4(const __m256 INST_DECODED_DECODERS_c1,
	                                           const __m256 INST_DECODED_DECODERS_c2,
	                                           const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    const __m256 C100 = _mm256_set1_ps(100.0f);
	                    const __m256 C05  = _mm256_set1_ps(0.5f);
	                    register __m256 t0,t1;
	                    register __m256 metric;
	                    t0 = _mm256_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm256_sub_ps(INST_DECODED_DECODERS_c1,
	                                       INST_DECODED_DECODERS_c2);
	                    metric = _mm256_mul_ps(C100,
	                                   _mm256_div_ps(t1,t0));
	                    return (metric);                                        
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_decoder0_only_ymm8r4(const float * __restrict pINST_DECODED_DECODERS_c1,
	                                           const float * __restrict pINST_DECODED_DECODERS_c2,
	                                           const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    register __m256 INST_DECODED_DECODERS_c1 = 
	                                        _mm256_loadu_ps(&pINST_DECODED_DECODERS_c1[0]);
	                    register __m256 INST_DECODED_DECODERS_c2 = 
	                                        _mm256_loadu_ps(&pINST_DECODED_DECODERS_c2[0]);
	                    register __m256 CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                    const __m256 C100 = _mm256_set1_ps(100.0f);
	                    const __m256 C05  = _mm256_set1_ps(0.5f);
	                    register __m256 t0,t1;
	                    register __m256 metric;
	                    t0 = _mm256_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm256_sub_ps(INST_DECODED_DECODERS_c1,
	                                       INST_DECODED_DECODERS_c2);
	                    metric = _mm256_mul_ps(C100,
	                                   _mm256_div_ps(t1,t0));
	                    return (metric);                                        
	        }
	        
/*
      "MetricName": "DSB",
      "LegacyName": "metric_TMA_....DSB(%)",
      "ParentCategory": "Fetch_Bandwidth",
      "Level": 3,
      "BriefDescription": "This metric represents Core fraction of cycles in which CPU was likely limited due to DSB (decoded uop cache) fetch pipeline.  For example; inefficient utilization of the DSB cache structure or bank conflict when reading from it; are categorized here.",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dsb_ymm8r4(          const __m256 IDQ_DSB_CYCLES_ANY,
	                                           const __m256 IDQ_DSB_CYCLES_OK,
	                                           const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    const __m256 C100 = _mm256_set1_ps(100.0f);
	                    const __m256 C05  = _mm256_set1_ps(0.5f);
	                    register __m256 t0,t1;
	                    register __m256 metric;
	                    t0 = _mm256_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm256_sub_ps(IDQ_DSB_CYCLES_ANY,
	                                       IDQ_DSB_CYCLES_OK);
	                    metric = _mm256_mul_ps(C100,
	                                   _mm256_div_ps(t1,t0));
	                    return (metric);                                        
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dsb_only_ymm8r4(     const float * __restrict pIDQ.DSB_CYCLES_ANY,
	                                           const float * __restrict pIDQ.DSB_CYCLES_OK,
	                                           const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    register __m256 IDQ_DSB_CYCLES_ANY = 
	                                        _mm256_loadu_ps(&pIDQ_DSB_CYCLES_ANY[0]);
	                    register __m256 IDQ_DSB_CYCLES_OK = 
	                                        _mm256_loadu_ps(&pIDQ_DSB_CYCLES_OK[0]);
	                    register __m256 CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                    const __m256 C100 = _mm256_set1_ps(100.0f);
	                    const __m256 C05  = _mm256_set1_ps(0.5f);
	                    register __m256 t0,t1;
	                    register __m256 metric;
	                    t0 = _mm256_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm256_sub_ps(IDQ.DSB_CYCLES_ANY,
	                                       IDQ.DSB_CYCLES_OK);
	                    metric = _mm256_mul_ps(C100,
	                                   _mm256_div_ps(t1,t0));
	                    return (metric);                                        
	        }   
	        
/*
      "MetricName": "Backend_Bound",
      "LegacyName": "metric_TMA_Backend_Bound(%)",
      "Level": 1,
      "BriefDescription": "This category represents fraction of slots where no uops are being delivered due to a lack of required resources for accepting new uops in the Backend. Backend is the portion of the processor core where the out-of-order scheduler dispatches ready uops into their respective execution units; and once completed these uops get retired according to program order. For example; stalls due to data-cache misses or stalls due to the divider unit being overloaded are both categorized under Backend Bound. Backend Bound is further divided into two main categories: Memory Bound and Core Bound.",
      "UnitOfMeasure": "percent",    
*/ 

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline 
                   float spr_backend_bound_r4(const float PERF_METRICS_BACKEND_BOUND,
                                          const float PERF_METRICS_FRONTEND_BOUND,
                                          const float PERF_METRICS_BAD_SPECULATION,
                                          const float PERF_METRICS_RETIRING) {
                                          
                           constexpr C100 = 100.0f;
                           float bcda,t0;
                           float metric;
                           bcda = PERF_METRICS_FRONTEND_BOUND+
                                  PERF_METRICS_BAD_SPECULATION+
                                  PERF_METRICS_RETIRING+
                                  PERF_METRICS_BACKEND_BOUND;
                           t0   = PERF_METRICS_BACKEND_BOUND/bcda;
                           metric = C100*t0;
                           return (metric);                
                 }
                 
/*
       "MetricName": "L1_Bound",
      "LegacyName": "metric_TMA_....L1_Bound(%)",
      "ParentCategory": "Memory_Bound",
      "Level": 3,
      "BriefDescription": "This metric estimates how often the CPU was stalled without loads missing the L1 data cache.  The L1 data cache typically has the shortest latency.  However; in certain cases like loads blocked on older stores; a load might suffer due to high latency even though it is being satisfied by the L1. Another example is loads who miss in the TLB. These cases are characterized by execution unit stalls; while some non-completed demand load lives in the machine without having that demand load missing the L1 cache.",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1_bound_ymm8r4(const __m256 EXE_ACTIVITY_BOUND_ON_LOADS,
	                                      const __m256 MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                      const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0,t1;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(_mm256_sub_ps(EXE_ACTIVITY_BOUND_ON_LOADS,
	                                                   MEMORY_ACTIVITY_STALLS_L1D_MISS),
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_max_ps(t0,_mm256_setzero_ps());
	                  metric = _mm256_mul_ps(C100,t1);
	                  return (metric);                           
	          }
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_l1_bound_ymm8r4(const float * __restrict pEXE_ACTIVITY_BOUND_ON_LOADS,
	                                      const float * __restrict pMEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                      const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m256 EXE_ACTIVITY_BOUND_ON_LOADS = 
	                                              _mm256_loadu_ps(&pEXE_ACTIVITY_BOUND_ON_LOADS[0]);
	                  register __m256 MEMORY_ACTIVITY_STALLS_L1D_MISS =
	                                              _mm256_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L1D_MISS[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD     = 
	                                              _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0,t1;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(_mm256_sub_ps(EXE_ACTIVITY_BOUND_ON_LOADS,
	                                                   MEMORY_ACTIVITY_STALLS_L1D_MISS),
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_max_ps(t0,_mm256_setzero_ps());
	                  metric = _mm256_mul_ps(C100,t1);
	                  return (metric);                           
	          }
	          
/*
   "MetricName": "DTLB_Load",
      "LegacyName": "metric_TMA_......DTLB_Load(%)",
      "ParentCategory": "L1_Bound",
      "Level": 4,
      "BriefDescription": "This metric roughly estimates the fraction of cycles where the Data TLB (DTLB) was missed by load accesses. TLBs (Translation Look-aside Buffers) are processor caches for recently used entries out of the Page Tables that are used to map virtual- to physical-addresses by the operating system. This metric approximates the potential delay of demand loads missing the first-level data TLB (assuming worst case scenario with back to back misses to different pages). This includes hitting in the second-level TLB (STLB) as well as performing a hardware page walk on an STLB miss.",
      "UnitOfMeasure": "percent",   
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_load_ymm8r4(const __m256 DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                       const __m256 DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                       const __m256 CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                       const __m256 MEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                       const __m256 CPU_CLK_UNHALTED_THREAD) {
	                            
	                  const __m256 C7 = _mm256_set1_ps(7.0f);
	                  const __m256 C0 = _mm256_setzero_ps();
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 t2;
	                  register __m256 metric;
	                  t0 = _mm256_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t1 = _mm256_max_ps(_mm256_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm256_min_ps(t0,t1);
	                  metric = _mm256_mul_ps(C100,
	                                     _mm256_div_ps(t2,CPU_CLK_UNHALTED_THREAD));
	                  return (metric);      
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_dtlb_load_ymm8r4(const float * __restrict pDTLB_LOAD_MISSES_STLB_HIT_c1,
	                                       const float * __restrict pDTLB_LOAD_MISSES_WALK_ACTIVE,
	                                       const float * __restrict pCYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                       const float * __restrict pMEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                       const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                            
	                  const __m256 DTLB_LOAD_MISSES_STLB_HIT_c1 = 
	                                        _mm256_loadu_ps(&pDTLB_LOAD_MISSES_STLB_HIT_c1[0]);
	                  const __m256 DTLB_LOAD_MISSES_WALK_ACTIVE =
	                                        _mm256_loadu_ps(&pDTLB_LOAD_MISSES_WALK_ACTIVE[0]);
	                  const __m256 CYCLE_ACTIVITY_CYCLES_MEM_ANY =
	                                        _mm256_loadu_ps(&pCYCLE_ACTIVITY_CYCLES_MEM_ANY[0]);
	                  const __m256 MEMORY_ACTIVITY_CYCLES_L1D_MISS=
	                                        _mm256_loadu_ps(&pMEMORY_ACTIVITY_CYCLES_L1D_MISS[0]);
	                  const __m256 CPU_CLK_UNHALTED_THREAD        =
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C7 = _mm256_set1_ps(7.0f);
	                  const __m256 C0 = _mm256_setzero_ps();
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 t2;
	                  register __m256 metric;
	                  t0 = _mm256_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t1 = _mm256_max_ps(_mm256_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm256_min_ps(t0,t1);
	                  metric = _mm256_mul_ps(C100,
	                                     _mm256_div_ps(t2,CPU_CLK_UNHALTED_THREAD));
	                  return (metric);      
	         }
	          
/*
       "MetricName": "Load_STLB_Hit",
      "LegacyName": "metric_TMA_........Load_STLB_Hit(%)",
      "ParentCategory": "DTLB_Load",
      "Level": 5,
      "BriefDescription": "This metric roughly estimates the fraction of cycles where the (first level) DTLB was missed by load accesses, that later on hit in second-level TLB (STLB)",
      "UnitOfMeasure": "percent",
*/       

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_load_stlb_hit_ymm8r4(const __m256 DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                          const __m256 DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                          const __m256 CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                          const __m256 MEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                          const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 C7 = _mm256_set1_ps(7.0f);
	                  const __m256 C0 = _mm256_setzero_ps();
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 t2;
	                  register __m256 t3
	                  register __m256 metric;  
	                  t0 = _mm256_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t3 = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_max_ps(_mm256_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm256_div_ps(_mm256_min_ps(t0,t1),
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,
	                                     _mm256_sub_ps(t2,t3));
	                  return (metric);
	                                                    
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_load_stlb_hit_ymm8r4(const float * __restrict pDTLB_LOAD_MISSES_STLB_HIT_c1,
	                                          const float * __restrict pDTLB_LOAD_MISSES_WALK_ACTIVE,
	                                          const float * __restrict pCYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                          const float * __restrict pMEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                          const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 DTLB_LOAD_MISSES_STLB_HIT_c1 = 
	                                        _mm256_loadu_ps(&pDTLB_LOAD_MISSES_STLB_HIT_c1[0]);
	                  const __m256 DTLB_LOAD_MISSES_WALK_ACTIVE =
	                                        _mm256_loadu_ps(&pDTLB_LOAD_MISSES_WALK_ACTIVE[0]);
	                  const __m256 CYCLE_ACTIVITY_CYCLES_MEM_ANY =
	                                        _mm256_loadu_ps(&pCYCLE_ACTIVITY_CYCLES_MEM_ANY[0]);
	                  const __m256 MEMORY_ACTIVITY_CYCLES_L1D_MISS=
	                                        _mm256_loadu_ps(&pMEMORY_ACTIVITY_CYCLES_L1D_MISS[0]);
	                  const __m256 CPU_CLK_UNHALTED_THREAD        =
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C7 = _mm256_set1_ps(7.0f);
	                  const __m256 C0 = _mm256_setzero_ps();
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 t2;
	                  register __m256 t3
	                  register __m256 metric;  
	                  t0 = _mm256_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t3 = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_max_ps(_mm256_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm256_div_ps(_mm256_min_ps(t0,t1),
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,
	                                     _mm256_sub_ps(t2,t3));
	                  return (metric);
	                                                    
	         } 
	         
/*
    "MetricName": "Load_STLB_Miss",
      "LegacyName": "metric_TMA_........Load_STLB_Miss(%)",
      "ParentCategory": "DTLB_Load",
      "Level": 5,
      "BriefDescription": "This metric estimates the fraction of cycles where the Second-level TLB (STLB) was missed by load accesses, performing a hardware page walk",
      "UnitOfMeasure": "percent",
*/	

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_load_stlb_miss_ymm8r4(const __m256 DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                            const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                                 
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_load_stlb_miss_ymm8r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_ACTIVE,
	                                            const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m256 DTLB_LOAD_MISSES_WALK_ACTIVE = 
	                                           _mm256_loadu_ps(&pDTLB_LOAD_MISSES_WALK_ACTIVE[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD      =
	                                           _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,t0);
	                  return (metric);                                 
	         }  
	         
/*
   "MetricName": "Store_Fwd_Blk",
      "LegacyName": "metric_TMA_......Store_Fwd_Blk(%)",
      "ParentCategory": "L1_Bound",
      "Level": 4,
      "BriefDescription": "This metric roughly estimates fraction of cycles when the memory subsystem had loads blocked since they could not forward data from earlier (in program order) overlapping stores. To streamline memory operations in the pipeline; a load can avoid waiting for memory if a prior in-flight store is writing the data that the load wants to read (store forwarding process). However; in some cases the load may be blocked for a significant time pending the store forward. For example; when the prior store is writing a smaller region than the load is reading.",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 store_fwd_blk_ymm8r4(const __m256 LD_BLOCKS_STORE_FORWARD,
	                                       const __m256 CPU_CLK_UNHALTED_THREAD) {
	                                       
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C13  = _mm256_set1_ps(13.0f);
	                  const __m256 C1   = _mm256_set1_ps(1.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 metric;                   
	                  t0 = _mm256_div_ps(LD_BLOCKS_STORE_FORWARD,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_min_ps(_mm256_mul_ps(C13,t0),C1);
	                  metric = _mm256_mul_ps(C100,t1);
	                  return (metric);
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 store_fwd_blk_ymm8r4(const float * __restrict pLD_BLOCKS_STORE_FORWARD,
	                                       const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                                       
	                  register __m256 LD_BLOCKS_STORE_FORWARD = 
	                                           _mm256_loadu_ps(&pLD_BLOCKS_STORE_FORWARD[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD =
	                                           _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C13  = _mm256_set1_ps(13.0f);
	                  const __m256 C1   = _mm256_set1_ps(1.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 metric;                   
	                  t0 = _mm256_div_ps(LD_BLOCKS_STORE_FORWARD,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_min_ps(_mm256_mul_ps(C13,t0),C1);
	                  metric = _mm256_mul_ps(C100,t1);
	                  return (metric);
	         }	
	         
/*
      "MetricName": "Split_Loads",
      "LegacyName": "metric_TMA_......Split_Loads(%)",
      "ParentCategory": "L1_Bound",
      "Level": 4,
      "BriefDescription": "This metric estimates fraction of cycles handling memory load split accesses - load that cross 64-byte cache line boundary. ",
      "UnitOfMeasure": "percent",
*/   

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_split_loads_ymm8r4(const __m256 L1D_PEND_MISS_PENDING,
	                                         const __m256 MEM_LOAD_COMPLETED_L1_MISS_ANY,
	                                         const __m256 LD_BLOCKS_NO_SR,
	                                         const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C1   = _mm256_set1_ps(1.0f);  
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 t2;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(L1_PEND_MISS_PENDING,
	                                     MEM_LOAD_COMPLETED_L1_MISS_ANY);
	                  t1 = _mm256_div_ps(LD_BLOCKS_NO_SR,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t2 = _mm256_min_ps(_mm256_mul_ps(t0,t1),C1);
	                  metric = _mm256_mul_ps(C100,t2);
	                  return (metric);                          
	         }    
	         
	              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_split_loads_ymm8r4(const float * __restrict pL1D_PEND_MISS_PENDING,
	                                         const float * __restrict pMEM_LOAD_COMPLETED_L1_MISS_ANY,
	                                         const float * __restrict pLD_BLOCKS_NO_SR,
	                                         const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m256 L1D_PEND_MISS_PENDING = 
	                                          _mm256_loadu_ps(&pL1D_PEND_MISS_PENDING[0]);
	                  register __m256 MEM_LOAD_COMPLETED_L1_MISS_ANY = 
	                                          _mm256_loadu_ps(&pMEM_LOAD_COMPLETED_L1_MISS_ANY[0]);
	                  register __m256 LD_BLOCKS_NO_SR        = 
	                                          _mm256_loadu_ps(&pLD_BLOCKS_NO_SR[0]);
	                  register __m256 CPU_CLK_UNHALTED_THREAD=
	                                          _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  const __m256 C1   = _mm256_set1_ps(1.0f);  
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 t2;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(L1_PEND_MISS_PENDING,
	                                     MEM_LOAD_COMPLETED_L1_MISS_ANY);
	                  t1 = _mm256_div_ps(LD_BLOCKS_NO_SR,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t2 = _mm256_min_ps(_mm256_mul_ps(t0,t1),C1);
	                  metric = _mm256_mul_ps(C100,t2);
	                  return (metric);                          
	         } 
	         
/*
        "MetricName": "FB_Full",
      "LegacyName": "metric_TMA_......FB_Full(%)",
      "ParentCategory": "L1_Bound",
      "Level": 4,
      "BriefDescription": "This metric does a *rough estimation* of how often L1D Fill Buffer unavailability limited additional L1D miss memory access requests to proceed. The higher the metric value; the deeper the memory hierarchy level the misses are satisfied from (metric values >1 are valid). Often it hints on approaching bandwidth limits (to L2 cache; L3 cache or external memory).",
      "UnitOfMeasure": "percent",
*/

                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m256 spr_fb_full_ymm8r4(const __m256 L1D_PEND_MISS_FB_FULL,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	                  
	                   const __m256 C100 = _mm256_set1_ps(100.0f);
	                   register __m256 t0;
	                   register __m256 metric;
	                   t0 = _mm256_div_ps(L1D_PEND_MISS_FB_FULL,
	                                      CPU_CLK_UNHALTED_THREAD);
	                   metric = _mm256_mul_ps(C100,t0);
	                   return (metric);                              
	       }  
	       
	                                               
	       	__ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_fb_full_ymm8r4(const float * __restrict pL1D_PEND_MISS_FB_FULL,
	                                  const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                   register __m256 L1D_PEND_MISS_FB_FULL = 
	                                           _mm256_loadu_ps(&pL1D_PEND_MISS_FB_FULL[0]);
	                   register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                           _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                   const __m256 C100 = _mm256_set1_ps(100.0f);
	                   register __m256 t0;
	                   register __m256 metric;
	                   t0 = _mm256_div_ps(L1D_PEND_MISS_FB_FULL,
	                                      CPU_CLK_UNHALTED_THREAD);
	                   metric = _mm256_mul_ps(C100,t0);
	                   return (metric);                              
	       }  
	                
/*
    "MetricName": "L2_Bound",
      "LegacyName": "metric_TMA_....L2_Bound(%)",
      "ParentCategory": "Memory_Bound",
      "Level": 3,
      "BriefDescription": "This metric estimates how often the CPU was stalled due to L2 cache accesses by loads.  Avoiding cache misses (i.e. L1 misses/L2 hits) can improve the latency and increase performance.",
      "UnitOfMeasure": "percent",   
*/

          	__ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256	 spr_l2_bound_ymm8r4(const __m256 MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                     const __m256 MEMORY_ACTIVITY_STALLS_L2_MISS ,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	                 
	                 const __m256 C100 = _mm256_set1_ps(100.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 metric;
	                 t0 = _mm256_sub_ps(MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L2_MISS);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm256_mul_ps(C100,t1);
	                 return (metric);                    
	       }  
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256	 spr_l2_bound_ymm8r4(const float * __restrict pMEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                     const float * __restrict pMEMORY_ACTIVITY_STALLS_L2_MISS ,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                 register __m256 MEMORY_ACTIVITY_STALLS_L1D_MISS = 
	                                                _mm256_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L1D_MISS[0]);
	                 register __m256 MEMORY_ACTIVITY_STALLS_L2_MISS  =
	                                                _mm256_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L2_MISS[0]);
	                 register __m256 CPU_CLK_UNHALTED_THREAD         =
	                                                _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                 const __m256 C100 = _mm256_set1_ps(100.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 metric;
	                 t0 = _mm256_sub_ps(MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L2_MISS);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm256_mul_ps(C100,t1);
	                 return (metric);                    
	       }   
	      
/*
    "MetricName": "L3_Bound",
      "LegacyName": "metric_TMA_....L3_Bound(%)",
      "ParentCategory": "Memory_Bound",
      "Level": 3,
      "BriefDescription": "This metric estimates how often the CPU was stalled due to loads accesses to L3 cache or contended with a sibling Core.  Avoiding cache misses (i.e. L2 misses/L3 hits) can improve the latency and increase performance.",
      "UnitOfMeasure": "percent",
*/   

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256	 spr_l3_bound_ymm8r4(const __m256 MEMORY_ACTIVITY_STALLS_L2_MISS,
	                                     const __m256 MEMORY_ACTIVITY_STALLS_L3_MISS ,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	                 
	                 const __m256 C100 = _mm256_set1_ps(100.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 metric;
	                 t0 = _mm256_sub_ps(MEMORY_ACTIVITY_STALLS_L2_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L3_MISS);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm256_mul_ps(C100,t1);
	                 return (metric);                    
	       }  
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256	 spr_l3_bound_ymm8r4(const float * __restrict pMEMORY_ACTIVITY_STALLS_L2_MISS,
	                                     const float * __restrict pMEMORY_ACTIVITY_STALLS_L3_MISS ,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                 register __m256 MEMORY_ACTIVITY_STALLS_L2_MISS = 
	                                                _mm256_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L2_MISS[0]);
	                 register __m256 MEMORY_ACTIVITY_STALLS_L3_MISS  =
	                                                _mm256_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L3_MISS[0]);
	                 register __m256 CPU_CLK_UNHALTED_THREAD         =
	                                                _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                 const __m256 C100 = _mm256_set1_ps(100.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 metric;
	                 t0 = _mm256_sub_ps(MEMORY_ACTIVITY_STALLS_L2_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L3_MISS);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm256_mul_ps(C100,t1);
	                 return (metric);                    
	       }   
	       
/*
     "MetricName": "SQ_Full",
      "LegacyName": "metric_TMA_......SQ_Full(%)",
      "ParentCategory": "L3_Bound",
      "Level": 4,
      "BriefDescription": "This metric measures fraction of cycles where the Super Queue (SQ) was full taking into account all request-types and both hardware SMT threads (Logical Processors). The Super Queue is used for requests to access the L2 cache or to go out to the Uncore.",
      "UnitOfMeasure": "percent",
*/	

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_sq_full_ymm8r4(const __m256 XQ_FULL_CYCLES,
	                                  const __m256 L1D_PEND_MISS_L2_STALLS,
	                                  const __m256 CPU_CLK_UNHALTED_THREAD) {
	               
	                 const __m256 C100 = _mm256_set1_ps(100.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 metric;
	                 t0 = _mm256_add_ps(XQ_FULL_CYCLES,
	                                    L1D_PEND_MISS_L2_STALLS);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm256_mul_ps(C100,t1);
	                 return (metric);                       
	       } 
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_sq_full_ymm8r4(const float * __restrict pXQ_FULL_CYCLES,
	                                  const float * __restrict pL1D_PEND_MISS_L2_STALLS,
	                                  const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                 register __m256 XQ_FULL_CYCLES = 
	                                   _mm256_loadu_ps(&pXQ_FULL_CYCLES[0]);
	                 register __m256 L1D_PEND_MISS_L2_STALLS = 
	                                   _mm256_loadu_ps(&pL1D_PEND_MISS_L2_STALLS[0]);
	                 register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                   _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                 const __m256 C100 = _mm256_set1_ps(100.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 metric;
	                 t0 = _mm256_add_ps(XQ_FULL_CYCLES,
	                                    L1D_PEND_MISS_L2_STALLS);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm256_mul_ps(C100,t1);
	                 return (metric);                       
	       }  
	       
/*
     "MetricName": "MEM_Bandwidth",
      "LegacyName": "metric_TMA_......MEM_Bandwidth(%)",
      "ParentCategory": "DRAM_Bound",
      "Level": 4,
      "BriefDescription": "This metric estimates fraction of cycles where the core's performance was likely hurt due to approaching bandwidth limits of external memory (DRAM).  The underlying heuristic assumes that a similar off-core traffic is generated by all IA cores. This metric does not aggregate non-data-read requests by this logical processor; requests from other IA Logical Processors/Physical Cores/sockets; or other non-IA devices like GPU; hence the maximum external memory bandwidth limits may or may not be approached when this metric is flagged (see Uncore counters for that).",
      "UnitOfMeasure": "percent",
*/ 

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mem_bw_ymm8r4(const __m256 CPU_CLK_UNHALTED_THREAD,
	                                 const __m256 OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4) {
	               
	               const __m256 C100 = _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                  OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4);
	               t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t1);
	               return (metric);                         
	       }  
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mem_bw_ymm8r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                 const float * __restrict pOFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4) {
	               
	               register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                      _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               register __m256 OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4 = 
	                                      _mm256_loadu_ps(&pOFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4[0]);
	               const __m256 C100 = _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                  OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4);
	               t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t1);
	               return (metric);                         
	       } 
	       
/*
     "MetricName": "MBA_Stalls",
      "LegacyName": "metric_TMA_........MBA_Stalls(%)",
      "ParentCategory": "MEM_Bandwidth",
      "Level": 5,
      "BriefDescription": "This metric estimates fraction of cycles where the core's performance was likely hurt due to memory bandwidth Allocation feature (RDT's memory bandwidth throttling).",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mba_stalls_ymm8r4(const __m256 INT_MISC_MBA_STALLS,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	               
	               const __m256 C100 = _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps(INT_MISC_MBA_STALLS,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	      }  
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mba_stalls_ymm8r4(const float * __restrict pINT_MISC_MBA_STALLS,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	               
	               register __m256 INT_MISC_MBA_STALLS = 
	                                       _mm256_loadu_ps(&pINT_MISC_MBA_STALLS[0]);
	               register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m256 C100 = _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps(INT_MISC_MBA_STALLS,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	      } 
	      
/*
      "MetricName": "MEM_Latency",
      "LegacyName": "metric_TMA_......MEM_Latency(%)",
      "ParentCategory": "DRAM_Bound",
      "Level": 4,
      "BriefDescription": "This metric estimates fraction of cycles where the performance was likely hurt due to latency from external memory (DRAM).  This metric does not aggregate requests from other Logical Processors/Physical Cores/sockets (see Uncore counters for that).",
      "UnitOfMeasure": "percent",
*/   

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mem_latency_ymm8r4(const __m256 CPU_CLK_UNHALTED_THREAD,
	                                      const __m256 OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DATA_RD,
	                                      const __m256 OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4) {
	               
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(_mm256_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                                   OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DATA_RD),
	                                                   CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_div_ps(_mm256_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                                   OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4),
	                                                   CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,
	                                 _mm256_sub_ps(t0,t1));
	                  return (metric);                             
	       }
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mem_latency_ymm8r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                      const float * __restrict pOFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DATA_RD,
	                                      const float * __restrict pOFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4) {
	               
	                  register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                         _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  register __m256 OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DATA_RD = 
	                                         _mm256_loadu_ps(&pOFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DATA_RD[0]);
	                  register __m256 OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4 = 
	                                         _mm256_loadu_ps(&pOFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4[0]);
	                  const __m256 C100 = _mm256_set1_ps(100.0f);
	                  register __m256 t0;
	                  register __m256 t1;
	                  register __m256 metric;
	                  t0 = _mm256_div_ps(_mm256_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                                   OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DATA_RD),
	                                                   CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm256_div_ps(_mm256_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                                   OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4),
	                                                   CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm256_mul_ps(C100,
	                                 _mm256_sub_ps(t0,t1));
	                  return (metric);                             
	       }
	       
/*
    "MetricName": "Store_Bound",
      "LegacyName": "metric_TMA_....Store_Bound(%)",
      "ParentCategory": "Memory_Bound",
      "Level": 3,
      "BriefDescription": "This metric estimates how often CPU was stalled  due to RFO store memory accesses; RFO store issue a read-for-ownership request before the write. Even though store accesses do not typically stall out-of-order CPUs; there are few cases where stores can lead to actual stalls. This metric will be flagged should RFO stores be a bottleneck.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_bound_ymm8r4(const __m256 EXE_ACTIVITY_BOUND_ON_STORES,
	                                      const __m256 CPU_CLK_UNHALTED_THREAD) {
	               
	               const __m256 C100 = _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps( EXE_ACTIVITY_BOUND_ON_STORES,
	                                   CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	       }
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_bound_ymm8r4(const float * __restrict pEXE_ACTIVITY_BOUND_ON_STORES,
	                                      const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	               
	               register __m256 EXE_ACTIVITY_BOUND_ON_STORES = 
	                                           _mm256_loadu_ps(&pEXE_ACTIVITY_BOUND_ON_STORES[0]);
	               register __m256 CPU_CLK_UNHALTED_THREAD      =
	                                           _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m256 C100 = _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps( EXE_ACTIVITY_BOUND_ON_STORES,
	                                   CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	       }

/*
   "MetricName": "Store_Latency",
      "LegacyName": "metric_TMA_......Store_Latency(%)",
      "ParentCategory": "Store_Bound",
      "Level": 4,
      "BriefDescription": "This metric estimates fraction of cycles the CPU spent handling L1D store misses. Store accesses usually less impact out-of-order core performance; however; holding resources for longer time can lead into undesired implications (e.g. contention on L1D fill-buffer entries - see FB_Full)",
      "UnitOfMeasure": "percent",    
*/	

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_latency_ymm8r4(const __m256 MEM_STORE_RETIRED_L2_HIT,
	                                        const __m256 MEM_INST_RETIRED_LOCK_LOADS,
	                                        const __m256 MEM_INST_RETIRED_ALL_STORES,
	                                        const __m256 CPU_CLK_UNHALTED_THREAD,
	                                        const __m256 OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO) {
	                                        
	                const __m256 C100 = _mm256_set1_ps(100.0f);
	                const __m256 C10  = _mm256_set1_ps(10.0f);
	                const __m256 C1   = _mm256_set1_ps(1.0f);
	                register __m256 t0;
	                register __m256 t1;
	                register __m256 t2;
	                register __m256 t3;
	                register __m256 metric;
	                t0 = _mm256_min_ps( CPU_CLK_UNHALTED_THREAD,
	                                    OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO);
	                t1 = _mm256_sub_ps(C1,
	                               _mm256_div_ps(MEM_INST_RETIRED_LOCK_LOADS,
	                                             MEM_INST_RETIRED_ALL_STORES));
	                t2 = _mm256_fmadd_ps(_mm256_mul_ps(
	                                 MEM_STORE_RETIRED_L2_HIT,C10),t1,t1);
	                t3 = _mm256_div_ps(_mm256_mul_ps(t2,t0),
	                                         CPU_CLK_UNHALTED_THREAD);
	                metric = _mm256_mul_ps(C100,t3);
	                return (metric);                              
	      }
	      
	      
              
                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_latency_ymm8r4(const float * __restrict pMEM_STORE_RETIRED_L2_HIT,
	                                        const float * __restrict pMEM_INST_RETIRED_LOCK_LOADS,
	                                        const float * __restrict pMEM_INST_RETIRED_ALL_STORES,
	                                        const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                        const float * __restrict pOFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO) {
	                        
	                register __m256  MEM_STORE_RETIRED_L2_HIT = 
	                                          _mm256_loadu_ps(&pMEM_STORE_RETIRED_L2_HIT[0]);
	                register __m256  MEM_INST_RETIRED_LOCK_LOADS = 
	                                          _mm256_loadu_ps(&pMEM_INST_RETIRED_LOCK_LOADS[0]);
	                register __m256  MEM_INST_RETIRED_ALL_STORES = 
	                                          _mm256_loadu_ps(&pMEM_INST_RETIRED_ALL_STORES[0]);
	                register __m256  CPU_CLK_UNHALTED_THREAD     =
	                                          _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                register __m256  OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO = 
	                                          _mm256_loadu_ps(&pOFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO[0]);
	                const __m256 C100 = _mm256_set1_ps(100.0f);
	                const __m256 C10  = _mm256_set1_ps(10.0f);
	                const __m256 C1   = _mm256_set1_ps(1.0f);
	                register __m256 t0;
	                register __m256 t1;
	                register __m256 t2;
	                register __m256 t3;
	                register __m256 metric;
	                t0 = _mm256_min_ps( CPU_CLK_UNHALTED_THREAD,
	                                    OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO);
	                t1 = _mm256_sub_ps(C1,
	                               _mm256_div_ps(MEM_INST_RETIRED_LOCK_LOADS,
	                                             MEM_INST_RETIRED_ALL_STORES));
	                t2 = _mm256_fmadd_ps(_mm256_mul_ps(
	                                 MEM_STORE_RETIRED_L2_HIT,C10),t1,t1);
	                t3 = _mm256_div_ps(_mm256_mul_ps(t2,t0),
	                                         CPU_CLK_UNHALTED_THREAD);
	                metric = _mm256_mul_ps(C100,t3);
	                return (metric);                              
	      }
	      
/*
       "MetricName": "False_Sharing",
      "LegacyName": "metric_TMA_......False_Sharing(%)",
      "ParentCategory": "Store_Bound",
      "Level": 4,
      "BriefDescription": "This metric roughly estimates how often CPU was handling synchronizations due to False Sharing. False Sharing is a multithreading hiccup; where multiple Logical Processors contend on different data-elements mapped into the same cache line. ",
      "UnitOfMeasure": "percent",  
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256  spr_false_sharing_ymm8r4(const __m256 CPU_CLK_UNHALTED_THREAD,
	                                         const __m256 CPU_CLK_UNHALTED_REF_TSC,
	                                         const __m256 OCR_DEMAND_RFO_L3_HIT_SNOOP_HITM,
	                                         const float SYSTEM_TSC_FREQ,
	                                         const float dur_time) {
	                
	                const __m256 C0000000001 = 
	                             _mm256_set1_ps(0.000000001f);
	                const __m256 C001        =
	                             _mm256_set1_ps(0.001f);
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                const __m256 C80         =
	                             _mm256_set1_ps(80.0f);
	                const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	                register __m256 t0;
	                register __m256 t1;
	                register __m256 t2;
	                register __m256 t3;
	                register __m256 v1;
	                register __m256 v2;
	                register __m256 metric;
	                v1 = _mm256_set1_ps(SYSTEM_TSC_FREQ);
	                t0 = mm256_mul_ps(v1,C000000001);
	             	v2 = _mm256_set1_ps(dur_time);
	                t1 = _mm256_div_ps(CPU_CLK_UNHALTED_THREAD,
	                                   CPU_CLK_UNHALTED_REF_TSC);
	                t2 = _mm256_mul_ps(v2,C001);
	                t2 = _mm256_mul_ps(t0,
	                               _mm256_mul_ps(t1,t2));
	                t3 = _mm256_div_ps(OCR_DEMAND_RFO_L3_HIT_SNOOP_HITM,
	                                   CPU_CLK_UNHALTED_THREAD);
	                v1 = _mm256_mul_ps(C80,
	                           _mm256_mul_ps(t2,t3)); 
	                metric = _mm256_mul_ps(C100,
	                               _mm256_min_ps(v1,C1));
	                return (metric);         
	      }
	           
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256  spr_false_sharing_ymm8r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                         const float * __restrict pCPU_CLK_UNHALTED_REF_TSC,
	                                         const float * __restrict pOCR_DEMAND_RFO_L3_HIT_SNOOP_HITM,
	                                         const float SYSTEM_TSC_FREQ,
	                                         const float dur_time) {
	                
	                register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                register __m256 CPU_CLK_UNHALTED_REF_TSC = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_REF_TSC[0]);
	                register __m256 OCR_DEMAND_RFO_L3_HIT_SNOOP_HITM = 
	                                       _mm256_loadu_ps(&pOCR_DEMAND_RFO_L3_HIT_SNOOP_HITM[0]);
	                const __m256 C0000000001 = 
	                             _mm256_set1_ps(0.000000001f);
	                const __m256 C001        =
	                             _mm256_set1_ps(0.001f);
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                const __m256 C80         =
	                             _mm256_set1_ps(80.0f);
	                const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	                register __m256 t0;
	                register __m256 t1;
	                register __m256 t2;
	                register __m256 t3;
	                register __m256 v1;
	                register __m256 v2;
	                register __m256 metric;
	                v1 = _mm256_set1_ps(SYSTEM_TSC_FREQ);
	                t0 = mm256_mul_ps(v1,C000000001);
	             	v2 = _mm256_set1_ps(dur_time);
	                t1 = _mm256_div_ps(CPU_CLK_UNHALTED_THREAD,
	                                   CPU_CLK_UNHALTED_REF_TSC);
	                t2 = _mm256_mul_ps(v2,C001);
	                t2 = _mm256_mul_ps(t0,
	                               _mm256_mul_ps(t1,t2));
	                t3 = _mm256_div_ps(OCR_DEMAND_RFO_L3_HIT_SNOOP_HITM,
	                                   CPU_CLK_UNHALTED_THREAD);
	                v1 = _mm256_mul_ps(C80,
	                           _mm256_mul_ps(t2,t3)); 
	                metric = _mm256_mul_ps(C100,
	                               _mm256_min_ps(v1,C1));
	                return (metric);         
	      }   
	      
/*
      "MetricName": "Split_Stores",
      "LegacyName": "metric_TMA_......Split_Stores(%)",
      "ParentCategory": "Store_Bound",
      "Level": 4,
      "BriefDescription": "This metric represents rate of split store accesses.  Consider aligning your data to the 64-byte cache line granularity.",
      "UnitOfMeasure": "percent",
*/    

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_split_stores_ymm8r4(const __m256 MEM_INST_RETIRED_SPLIT_STORES,
	                                   const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	               
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 metric;
	                t0 = _mm256_div_ps(MEM_INST_RETIRED_SPLIT_STORES,
	                                   CPU_CLK_UNHALTED_DISTRIBUTED);
	                metric = _mm256_mul_ps(C100,t0);
	                return (metric);                            
	      } 
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_split_stores_ymm8r4(const float * __restrict pMEM_INST_RETIRED_SPLIT_STORES,
	                                       const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	               
	                register __m256 MEM_INST_RETIRED_SPLIT_STORES = 
	                                        _mm256_loadu_ps(&pMEM_INST_RETIRED_SPLIT_STORES[0]);
	                register __m256 CPU_CLK_UNHALTED_DISTRIBUTED  =
	                                        _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 metric;
	                t0 = _mm256_div_ps(MEM_INST_RETIRED_SPLIT_STORES,
	                                   CPU_CLK_UNHALTED_DISTRIBUTED);
	                metric = _mm256_mul_ps(C100,t0);
	                return (metric);                            
	      } 
	      
	      
	      
/*
      "MetricName": "Streaming_Stores",
      "LegacyName": "metric_TMA_......Streaming_Stores(%)",
      "ParentCategory": "Store_Bound",
      "Level": 4,
      "BriefDescription": "This metric estimates how often CPU was stalled  due to Streaming store memory accesses; Streaming store optimize out a read request required by RFO stores. Even though store accesses do not typically stall out-of-order CPUs; there are few cases where stores can lead to actual stalls. This metric will be flagged should Streaming stores be a bottleneck.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_stream_stores_ymm8r4(const __m256 OCR_STREAMING_WR_ANY_RESPONSE,
	                                        const __m256 CPU_CLK_UNHALTED_THREAD) {
	               
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               const __m256 C9          =
	                             _mm256_set1_ps(9.0f);
	               const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_div_ps(_mm512_mul_ps(C9,
	                                      OCR_STREAMING_WR_ANY_RESPONSE),
	                                                       CPU_CLK_UNHALTED_THREAD);
	               t2 = _mm256_min_ps(t0,C1);
	               metric = _mm256_mul_ps(C100,t2);
	               return (metric);                               
	      } 
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_stream_stores_ymm8r4(const __m256 OCR_STREAMING_WR_ANY_RESPONSE,
	                                        const __m256 CPU_CLK_UNHALTED_THREAD) {
	               
	               const __m256 OCR_STREAMING_WR_ANY_RESPONSE = 
	                                         _mm256_loadu_ps(&pOCR_STREAMING_WR_ANY_RESPONSE[0]);
	               const __m256 CPU_CLK_UNHALTED_THREAD       =
	                                         _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               const __m256 C9          =
	                             _mm256_set1_ps(9.0f);
	               const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_div_ps(_mm512_mul_ps(C9,
	                                      OCR_STREAMING_WR_ANY_RESPONSE),
	                                                       CPU_CLK_UNHALTED_THREAD);
	               t2 = _mm256_min_ps(t0,C1);
	               metric = _mm256_mul_ps(C100,t2);
	               return (metric);                               
	      } 
	      
/*
       "MetricName": "DTLB_Store",
      "LegacyName": "metric_TMA_......DTLB_Store(%)",
      "ParentCategory": "Store_Bound",
      "Level": 4,
      "BriefDescription": "This metric roughly estimates the fraction of cycles spent handling first-level data TLB store misses.  As with ordinary data caching; focus on improving data locality and reducing working-set size to reduce DTLB overhead.  Additionally; consider using profile-guided optimization (PGO) to collocate frequently-used data on the same page.  Try using larger page sizes for large amounts of frequently-used data.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_dtlb_store_ymm8r4(const __m256 DTLB_STORE_MISSES_STLB_HIT_c1,
	                                     const __m256 DTLB_STORE_MISSES_WALK_ACTIVE,
	                                     const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                                     
	                 const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                 const __m256 C7          = 
	                             _mm256_set1_ps(7.0f);
	                 const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 t2;
	                 register __m256 metric;
	                 t0 = _mm256_fmadd_ps(C7,DTLB_STORE_MISSES_STLB_HIT_c1,
	                                      DTLB_STORE_MISSES_WALK_ACTIVE);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_DISTRIBUTED);
	                 t2 = _mm256_min_ps(t1,C1);
	                 metric = _mm256_mul_ps(C100,t2);
	                 return (metric);        
	       }
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_dtlb_store_ymm8r4(const float * __restrict pDTLB_STORE_MISSES_STLB_HIT_c1,
	                                     const float * __restrict pDTLB_STORE_MISSES_WALK_ACTIVE,
	                                     const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                     
	                 register __m256 DTLB_STORE_MISSES_STLB_HIT_c1 = 
	                                           _mm256_loadu_ps(&pDTLB_STORE_MISSES_STLB_HIT_c1[0]);
	                 register __m256 DTLB_STORE_MISSES_WALK_ACTIVE =
	                                           _mm256_loadu_ps(&pDTLB_STORE_MISSES_WALK_ACTIVE[0]);
	                 register __m256 CPU_CLK_UNHALTED_DISTRIBUTED  =
	                                           _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);                
	                 const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                 const __m256 C7          = 
	                             _mm256_set1_ps(7.0f);
	                 const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 t2;
	                 register __m256 metric;
	                 t0 = _mm256_fmadd_ps(C7,DTLB_STORE_MISSES_STLB_HIT_c1,
	                                      DTLB_STORE_MISSES_WALK_ACTIVE);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_DISTRIBUTED);
	                 t2 = _mm256_min_ps(t1,C1);
	                 metric = _mm256_mul_ps(C100,t2);
	                 return (metric);        
	       }   
	       
/*
     "MetricName": "Store_STLB_Hit",
      "LegacyName": "metric_TMA_........Store_STLB_Hit(%)",
      "ParentCategory": "DTLB_Store",
      "Level": 5,
      "BriefDescription": "This metric roughly estimates the fraction of cycles where the TLB was missed by store accesses, hitting in the second-level TLB (STLB)",
      "UnitOfMeasure": "percent",  
*/   

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_stlb_hit_ymm8r4(const __m256 DTLB_STORE_MISSES_STLB_HIT_c1,
	                                     const __m256 DTLB_STORE_MISSES_WALK_ACTIVE,
	                                     const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                                     
	                 const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                 const __m256 C7          = 
	                             _mm256_set1_ps(7.0f);
	                 const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 t2;
	                 register __m256 t3;
	                 register __m256 metric;
	                 t0 = _mm256_fmadd_ps(C7,DTLB_STORE_MISSES_STLB_HIT_c1,
	                                      DTLB_STORE_MISSES_WALK_ACTIVE);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_DISTRIBUTED);
	                 t2 = _mm256_min_ps(t1,C1);
	                 t3 = _mm256_div_ps(DTLB_STORE_MISSES_WALK_ACTIVE,
	                                    CPU_CLK_UNHALTED_DISTRIBUTED);
	                 metric = _mm256_mul_ps(C100,
	                                    _mm256_sub_ps(t2,t3));
	                 return (metric);
	     }
	     
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_stlb_hit_ymm8r4(const float * __restrict pDTLB_STORE_MISSES_STLB_HIT_c1,
	                                     const float * __restrict pDTLB_STORE_MISSES_WALK_ACTIVE,
	                                     const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                          
	                 register __m256 DTLB_STORE_MISSES_STLB_HIT_c1 = 
	                                           _mm256_loadu_ps(&pDTLB_STORE_MISSES_STLB_HIT_c1[0]);
	                 register __m256 DTLB_STORE_MISSES_WALK_ACTIVE =
	                                           _mm256_loadu_ps(&pDTLB_STORE_MISSES_WALK_ACTIVE[0]);
	                 register __m256 CPU_CLK_UNHALTED_DISTRIBUTED  =
	                                           _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);                  
	                 const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                 const __m256 C7          = 
	                             _mm256_set1_ps(7.0f);
	                 const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	                 register __m256 t0;
	                 register __m256 t1;
	                 register __m256 t2;
	                 register __m256 t3;
	                 register __m256 metric;
	                 t0 = _mm256_fmadd_ps(C7,DTLB_STORE_MISSES_STLB_HIT_c1,
	                                      DTLB_STORE_MISSES_WALK_ACTIVE);
	                 t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_DISTRIBUTED);
	                 t2 = _mm256_min_ps(t1,C1);
	                 t3 = _mm256_div_ps(DTLB_STORE_MISSES_WALK_ACTIVE,
	                                    CPU_CLK_UNHALTED_DISTRIBUTED);
	                 metric = _mm256_mul_ps(C100,
	                                    _mm256_sub_ps(t2,t3));
	                 return (metric);
	     }
	     
/*
        "MetricName": "Store_STLB_Miss",
      "LegacyName": "metric_TMA_........Store_STLB_Miss(%)",
      "ParentCategory": "DTLB_Store",
      "Level": 5,
      "BriefDescription": "This metric estimates the fraction of cycles where the STLB was missed by store accesses, performing a hardware page walk",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_stlb_miss_ymm8r4(const __m256 DTLB_STORE_MISSES_WALK_ACTIVE,
	                                          const __m256 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                                          
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 metric;
	                t0 = _mm256_div_ps(DTLB_STORE_MISSES_WALK_ACTIVE,
	                                   CPU_CLK_UNHALTED_DISTRIBUTED);
	                metric = _mm256_mul_ps(C100,t0);
	                return (metric);                                  
	       }
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_store_stlb_miss_ymm8r4(const float * __restrict pDTLB_STORE_MISSES_WALK_ACTIVE,
	                                          const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                             
	                register __m256 DTLB_STORE_MISSES_WALK_ACTIVE = 
	                                          _mm256_loadu_ps(&pDTLB_STORE_MISSES_WALK_ACTIVE[0]);
	                register __m256 CPU_CLK_UNHALTED_DISTRIBUTED  =
	                                          _mm256_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);             
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 metric;
	                t0 = _mm256_div_ps(DTLB_STORE_MISSES_WALK_ACTIVE,
	                                   CPU_CLK_UNHALTED_DISTRIBUTED);
	                metric = _mm256_mul_ps(C100,t0);
	                return (metric);                                  
	       }
	       
/*
    "MetricName": "Core_Bound",
      "LegacyName": "metric_TMA_..Core_Bound(%)",
      "ParentCategory": "Backend_Bound",
      "Level": 2,
      "BriefDescription": "This metric represents fraction of slots where Core non-memory issues were of a bottleneck.  Shortage in hardware compute resources; or dependencies in software's instructions are both categorized under Core Bound. Hence it may indicate the machine ran out of an out-of-order resource; certain execution units are overloaded or dependencies in program's data- or instruction-flow are limiting the performance (e.g. FP-chained long-latency arithmetic operations).",
      "UnitOfMeasure": "percent",
*/	       

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
               	static inline
	        float spr_core_bound_r4(const float PERF_METRICS_BACKEND_BOUND,
	                                const float PERF_METRICS_FRONTEND_BOUND,
	                                const float PERF_METRICS_BAD_SPECULATION,
	                                const float PERF_METRICS_RETIRING,
	                                const float PERF_METRICS_MEMORY_BOUND) {
	                                
	              float bcda,t0,t1;
	              float metric;
	              bcda = PERF_METRICS_FRONTEND_BOUND+
	                     PERF_METRICS_BAD_SPECULATION+
	                     PERF_METRICS_RETIRING+
	                     PERF_METRICS_BACKEND_BOUND;
	              t0 =   PERF_METRICS_BACKEND_BOUND/
	                     bcda;
	              t1 =   PERF_METRICS_MEMORY_BOUND/
	                     bcda;
	              metric = 100.0f*(std::max(0.0f,t0-t1));
	              return (metric);                      
	     }
	     
/*
        "MetricName": "Divider",
      "LegacyName": "metric_TMA_....Divider(%)",
      "ParentCategory": "Core_Bound",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of cycles where the Divider unit was active. Divide and square root instructions are performed by the Divider unit and can take considerably longer latency than integer or Floating Point addition; subtraction; or multiplication.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_divider_active_ymm8r4(const __m256 ARITH_DIV_ACTIVE,
	                                         const __m256 CPU_CLK_UNHALTED_THREAD) {
	               
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 metric;
	                t0 = _mm256_div_ps(ARITH_DIV_ACTIVE,
	                                   CPU_CLK_UNHALTED_THREAD);
	                metric = _mm256_mul_ps(C100,t0);
	                return (metric);                                 
	      } 
	      
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_divider_active_ymm8r4(const float * __restrict pARITH_DIV_ACTIVE,
	                                         const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	               
	                register __m256 ARITH_DIV_ACTIVE = 
	                                     _mm256_loadu_ps(&pARITH_DIV_ACTIVE[0]);
	                register __m256 CPU_CLK_UNHALTED_THREAD = 
	                                     _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 metric;
	                t0 = _mm256_div_ps(ARITH_DIV_ACTIVE,
	                                   CPU_CLK_UNHALTED_THREAD);
	                metric = _mm256_mul_ps(C100,t0);
	                return (metric);                                 
	      } 
	    
/*
     "MetricName": "Ports_Utilized_0",
      "LegacyName": "metric_TMA_......Ports_Utilized_0(%)",
      "ParentCategory": "Ports_Utilization",
      "Level": 4,
      "BriefDescription": "This metric represents fraction of cycles CPU executed no uops on any execution port (Logical Processor cycles since ICL, Physical Core cycles otherwise). Long-latency instructions like divides may contribute to this metric.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_ports_utils_0_ymm8r4(const __m256 EXE_ACTIVITY_3_PORTS_UTIL_u0x80,
	                                        const __m256 CPU_CLK_UNHALTED_THREAD,
	                                        const __m256 RESOURCE_STALLS_SCOREBOARD,
	                                        const __m256 CYCLE_ACTIVITY_STALLS_TOTAL,
	                                        const __m256 EXE_ACTIVITY_BOUND_ON_LOADS) {
	                
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 t1;
	                register __m256 t2;
	                register __m256 metric;  
	                t0  = _mm256_div_ps(EXE_ACTIVITY_3_PORTS_UTIL_u0x80,
	                                    CPU_CLK_UNHALTED_THREAD);
	                t1  = _mm256_sub_ps(CYCLE_ACTIVITY_STALLS_TOTAL,
	                                    EXE_ACTIVITY_BOUND_ON_LOADS);
	                t2  = _mm256_add_ps(t0,
	                                _mm256_div_ps(RESOURCE_STALLS_SCOREBOARD,
	                                              CPU_CLK_UNHALTED_THREAD));
	                t0  = _mm256_div_ps(t1,CPU_CLK_UNHALTED_THREAD);
	                metric = _mm256_mul_ps(C100,
	                               _mm256_mul_ps(t2,t0));
	                return (metric);                              
	      }
	      
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_ports_utils_0_ymm8r4(const float * __restrict pEXE_ACTIVITY_3_PORTS_UTIL_u0x80,
	                                        const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                        const float * __restrict pRESOURCE_STALLS_SCOREBOARD,
	                                        const float * __restrict pCYCLE_ACTIVITY_STALLS_TOTAL,
	                                        const float * __restrict pEXE_ACTIVITY_BOUND_ON_LOADS) {
	                
	                register __m256 EXE_ACTIVITY_3_PORTS_UTIL_u0x80 = 
	                                            _mm256_loadu_ps(&pEXE_ACTIVITY_3_PORTS_UTIL_u0x80[0]);
	                register __m256 CPU_CLK_UNHALTED_THREAD         =
	                                            _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                register __m256 RESOURCE_STALLS_SCOREBOARD      =
	                                            _mm256_loadu_ps(&pRESOURCE_STALLS_SCOREBOARD[0]);
	                register __m256 CYCLE_ACTIVITY_STALLS_TOTAL     =
	                                            _mm256_loadu_ps(&pCYCLE_ACTIVITY_STALLS_TOTAL[0]);
	                register __m256 EXE_ACTIVITY_BOUND_ON_LOADS     =
	                                            _mm256_loadu_ps(&pEXE_ACTIVITY_BOUND_ON_LOADS[0]);
	                const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	                register __m256 t0;
	                register __m256 t1;
	                register __m256 t2;
	                register __m256 metric;  
	                t0  = _mm256_div_ps(EXE_ACTIVITY_3_PORTS_UTIL_u0x80,
	                                    CPU_CLK_UNHALTED_THREAD);
	                t1  = _mm256_sub_ps(CYCLE_ACTIVITY_STALLS_TOTAL,
	                                    EXE_ACTIVITY_BOUND_ON_LOADS);
	                t2  = _mm256_add_ps(t0,
	                                _mm256_div_ps(RESOURCE_STALLS_SCOREBOARD,
	                                              CPU_CLK_UNHALTED_THREAD));
	                t0  = _mm256_div_ps(t1,CPU_CLK_UNHALTED_THREAD);
	                metric = _mm256_mul_ps(C100,
	                               _mm256_mul_ps(t2,t0));
	                return (metric);                              
	      } 
	      
/*
     "MetricName": "Serializing_Operation",
      "LegacyName": "metric_TMA_........Serializing_Operation(%)",
      "ParentCategory": "Ports_Utilized_0",
      "Level": 5,
      "BriefDescription": "This metric represents fraction of cycles the CPU issue-pipeline was stalled due to serializing operations. Instructions like CPUID; WRMSR or LFENCE serialize the out-of-order execution which may limit performance.",
      "UnitOfMeasure": "percent",
*/	

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_serial_ops_ymm8r4(const __m256 RESOURCE_STALLS_SCOREBOARD,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	                                     
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps(RESOURCE_STALLS_SCOREBOARD,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	      }   
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_serial_ops_ymm8r4(const float * __restrict pRESOURCE_STALLS_SCOREBOARD,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                                     
	               register __m256 RESOURCE_STALLS_SCOREBOARD = 
	                                       _mm256_loadu_ps(&pRESOURCE_STALLS_SCOREBOARD[0]);
	               register __m256 CPU_CLKUNHALTED_THREAD     = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps(RESOURCE_STALLS_SCOREBOARD,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	      }   
	     
/*
    "MetricName": "Slow_Pause",
      "LegacyName": "metric_TMA_..........Slow_Pause(%)",
      "ParentCategory": "Serializing_Operation",
      "Level": 6,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to PAUSE Instructions.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_slow_pause_ymm8r4(const __m256 CPU_CLK_UNHALTED_PAUSE,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	                                     
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps(CPU_CLK_UNHALTED_PAUSE,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	      }   
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_slow_pause_ymm8r4(const float * __restrict pCPU_CLK_UNHALTED_PAUSE,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                                     
	               register __m256 CPU_CLK_UNHALTED_PAUSE = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_PAUSE[0]);
	               register __m256 CPU_CLKUNHALTED_THREAD     = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               register __m256 t0;
	               register __m256 metric;
	               t0 = _mm256_div_ps(CPU_CLK_UNHALTED_PAUSE,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,t0);
	               return (metric);                              
	      }  
	      
/*
     "MetricName": "Memory_Fence",
      "LegacyName": "metric_TMA_..........Memory_Fence(%)",
      "ParentCategory": "Serializing_Operation",
      "Level": 6,
      "BriefDescription": "This metric represents fraction of cycles the CPU was stalled due to LFENCE Instructions.",
      "UnitOfMeasure": "percent",
*/ 

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mem_fence_ymm8r4( const __m256 MISC2_RETIRED_LFENCE,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	                                     
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               const __m256 C13         = 
	                             _mm256_set1_ps(13.0f);
	               const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_mul_ps(C13,MISC2_RETIRED_LFENCE);
	               t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,
	                                  _mm256_min_ps(t1,C1));
	               return (metric);                              
	      }   
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mem_fence_ymm8r4(const float * __restrict pMISC2_RETIRED_LFENCE,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                                     
	               register __m256 MISC2_RETIRED_LFENCE = 
	                                       _mm256_loadu_ps(&pMISC2_RETIRED_LFENCE[0]);
	               register __m256 CPU_CLKUNHALTED_THREAD     = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               const __m256 C13         = 
	                             _mm256_set1_ps(13.0f);
	               const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_mul_ps(C13,MISC2_RETIRED_LFENCE);
	               t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,
	                                  _mm256_min_ps(t1,C1));
	               return (metric);                              
	      }  
	      
/*
     "MetricName": "Mixing_Vectors",
      "LegacyName": "metric_TMA_........Mixing_Vectors(%)",
      "ParentCategory": "Ports_Utilized_0",
      "Level": 5,
      "BriefDescription": "The Mixing_Vectors metric gives the percentage of injected blend uops out of all uops issued. Usually a Mixing_Vectors over 5% is worth investigating. Read more in Appendix B1 of the Optimizations Guide for this topic.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mix_vectors_ymm8r4( const __m256 ASSISTS_SSE_AVX_MIX,
	                                     const __m256 CPU_CLK_UNHALTED_THREAD) {
	                                     
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               const __m256 C160         = 
	                             _mm256_set1_ps(160.0f);
	               const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_mul_ps(C160,ASSISTS_SSE_AVX_MIX);
	               t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,
	                                  _mm256_min_ps(t1,C1));
	               return (metric);                              
	      }   
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m256 spr_mix_vectors_ymm8r4(const float * __restrict pASSISTS_SSE_AVX_MIX,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                                     
	               register __m256 ASSISTS_SSE_AVX_MIX = 
	                                       _mm256_loadu_ps(&pASSISTS_SSE_AVX_MIX[0]);
	               register __m256 CPU_CLKUNHALTED_THREAD     = 
	                                       _mm256_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m256 C100        =
	                             _mm256_set1_ps(100.0f);
	               const __m256 C160         = 
	                             _mm256_set1_ps(160.0f);
	               const __m256 C1          =
	                             _mm256_set1_ps(1.0f);
	               register __m256 t0;
	               register __m256 t1;
	               register __m256 metric;
	               t0 = _mm256_mul_ps(C160,ASSISTS_SSE_AVX_MIX);
	               t1 = _mm256_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm256_mul_ps(C100,
	                                  _mm256_min_ps(t1,C1));
	               return (metric);                              
	      }  
	      
	      


} // gms








#endif /*__GMS_SPR_METRICS_YMM8R4*/
