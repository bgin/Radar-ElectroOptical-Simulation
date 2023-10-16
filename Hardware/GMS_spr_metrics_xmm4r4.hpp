

/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!*/

#ifndef __GMS_SPR_METRICS_XMM4R4_HPP__
#define __GMS_SPR_METRICS_XMM4R4_HPP__



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
	           __m128 spr_cpu_op_freq_xmm4r4(const __m128 CPU_CLK_UNHALTED_THREAD,
	                                         const __m128 CPU_CLK_UNHALTED_REF_TSC,
	                                         const __m128 SYSTEM_TSC_FREQ) {
	                  
	                  const __m128 C000000001 = _mm128_set1_ps(1.0e-09f);
	                  register __m128 t0,metric;
	                  t0 = _mm128_div_ps(CPU_CLK_UNHALTED_THREAD,
	                                     CPU_CLK_UNHALTED_REF_TSC);
	                  metric = _mm128_mul_ps(C000000001,
	                                     _mm128_mul_ps(t0,SYSTEM_TSC_FREQ));
	                  return (metric);
	                                             
	          }
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_cpu_op_freq_xmm4r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                         const float * __restrict pCPU_CLK_UNHALTED_REF_TSC,
	                                         const float * __restrict pSYSTEM_TSC_FREQ) {
	                  
	                  register __m128 CPU_CLK_UNHALTED_THREAD   =  
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  register __m128 CPU_CLK_UNHALTED_REF_TSC  =  
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_REF_TSC[0]);
	                  register __m128 SYSTEM_TSC_FREQ           =  
	                                        _mm128_loadu_ps(&pSYSTEM_TSC_FREQ[0]);
	                  const __m128 C000000001 = _mm128_set1_ps(1.0e-09f);
	                  register __m128 t0,metric;
	                  t0 = _mm128_div_ps(CPU_CLK_UNHALTED_THREAD,
	                                     CPU_CLK_UNHALTED_REF_TSC);
	                  metric = _mm128_mul_ps(C000000001,
	                                     _mm128_mul_ps(t0,SYSTEM_TSC_FREQ));
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
	           __m128 spr_cpu_util_xmm4r4(const __m128 CPU_CLK_UNHALTED_REF_TSC,
	                                      const __m128 TSC_STAMP) {
	                                      
	                  register __m128 metric;
	                  result = _mm128_div_ps(CPU_CLK_UNHALTED_REF_TSC,TSC_STAMP);
	                  return (metric);                             
	        }  
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_cpu_util_xmm4r4(const float * __restrict pCPU_CLK_UNHALTED_REF_TSC,
	                                      const float * __restrict pTSC_STAMP) {
	                                      
	                  register __m128 CPU_CLK_UNHALTED_REF_TSC = 
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_REF_TSC[0]);
	                  register __m128 TSC_STAMP                =
	                                        _mm128_loadu_ps(&pTSC_STAMP[0]);   
	                  register __m128 metric;
	                  result = _mm128_div_ps(CPU_CLK_UNHALTED_REF_TSC,TSC_STAMP);
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
	           __m128 spr_cpi_xmm4r4(const __m128 CPU_CLK_UNHALTED_THREAD,
	                                 const __m128 INST_RETIRED_ANY) {
	                                 
	                  register __m128 metric;
	                  metric = _mm128_div_ps(CPU_CLK_UNHALTED_THREAD,
	                                         INST_RETIRED_ANY);
	                  return (metric);                       
	         }
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_cpi_xmm4r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                 const float * __restrict pINST_RETIRED_ANY) {
	                       
	                  register __m128 CPU_CLK_UNHALTED_THREAD = 
	                                         _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  register __m128 INST_RETIRED_ANY        = 
	                                         _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);          
	                  register __m128 metric;
	                  metric = _mm128_div_ps(CPU_CLK_UNHALTED_THREAD,
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
	           __m128 spr_loads_per_instr_xmm4r4(const __m128 MEM_INST_RETIRED_ALL_LOADS,
	                                             const __m128 INST_RETIRED_ANY) {
	                   
	                   register __m128 metric;
	                   metric = _mm128_div_ps(MEM_INST_RETIRED_ALL_LOADS,
	                                          INST_RETIRED_ANY);
	                   return (metric);
	                                                                       
	          }
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_loads_per_instr_xmm4r4(const float * __restrict pMEM_INST_RETIRED_ALL_LOADS,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                   
	                   register __m128 MEM_INST_RETIRED_ALL_LOADS = 
	                                           _mm128_loadu_ps(&pMEM_INST_RETIRED_ALL_LOADS[0]);
	                   register __m128 INST_RETIRED_ANY           =
	                                           _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                   register __m128 metric;
	                   metric = _mm128_div_ps(MEM_INST_RETIRED_ALL_LOADS,
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
	           __m128 spr_stores_per_instr_xmm4r4(const __m128 MEM_INST_RETIRED_ALL_STORES,
	                                              const __m128 INST_RETIRED_ANY) {
	                   
	                   register __m128 metric;
	                   metric = _mm128_div_ps(MEM_INST_RETIRED_ALL_STORES,
	                                          INST_RETIRED_ANY);
	                   return (metric);
	                                                                       
	          }
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_stores_per_instr_xmm4r4(const float * __restrict pMEM_INST_RETIRED_ALL_STORES,
	                                              const float * __restrict pINST_RETIRED_ANY) {
	                   
	                   register __m128 MEM_INST_RETIRED_ALL_STORES = 
	                                           _mm128_loadu_ps(&pMEM_INST_RETIRED_ALL_STORES[0]);
	                   register __m128 INST_RETIRED_ANY            =
	                                           _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                   register __m128 metric;
	                   metric = _mm128_div_ps(MEM_INST_RETIRED_ALL_STORES,
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
	           __m128 spr_l1d_mpi_xmm4r4(const __m128 L1D_REPLACEMENT,
	                                     const __m128 INST_RETIRED_ANY) {
	                                     
	                  register __m128 metric;
	                  metric = _mm128_div_ps(L1D_REPLACEMENT,
	                                         INST_RETIRED_ANY);
	                  return (metric);                            
	        }  
	        
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l1d_mpi_xmm4r4(const float * __restrict pL1D_REPLACEMENT,
	                                     const float * __restrict pINST_RETIRED_ANY) {
	                                     
	                  register __m128 L1D_REPLACEMENT  = 
	                                     _mm128_loadu_ps(&pL1D_REPLACEMENT[0]);
	                  register __m128 INST_RETIRED_ANY =
	                                     _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(L1D_REPLACEMENT,
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
	           __m128 spr_l1d_ddrh_pi_xmm4r4(const __m128 MEM_LOAD_RETIRED_L1_HIT,
	                                         const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(MEM_LOAD_RETIRED_L1_HIT,
	                                         INST_RETIRED_ANY);
	                  return (metric);                       
	                                                      
	         }
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l1d_ddrh_pi_xmm4r4(const float * __restrict pMEM_LOAD_RETIRED_L1_HIT,
	                                         const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 MEM_LOAD_RETIRED_L1_HIT = 
	                                        _mm128_loadu_ps(&pMEM_LOAD_RETIRED_L1_HIT[0]);
	                  register __m128 INST_RETIRED_ANY        =
	                                        _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(MEM_LOAD_RETIRED_L1_HIT,
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
	           __m128 spr_l1i_crmp_pi_xmm4r4(const __m128 L2_RQSTS_ALL_CODE_RD,
	                                         const __m128 INST_RETIRED_ANY) {
	                                         
	                   register __m128 metric;
	                   metric = _mm128_div_ps(L2_RQSTS_ALL_CODE_RD,
	                                          INST_RETIRED_ANY);
	                   return (metric);                               
	         }
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l1i_crmp_pi_xmm4r4(const float * __restrict pL2_RQSTS_ALL_CODE_RD,
	                                         const float * __restrict pINST_RETIRED_ANY) {
	                                 
	                   register __m128 L2_RQSTS_ALL_CODE_RD = 
	                                           _mm128_loadu_ps(&pL2_RQSTS_ALL_CODE_RD[0]);
	                   register __m128 INST_RETIRED_ANY     =
	                                           _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);        
	                   register __m128 metric;
	                   metric = _mm128_div_ps(L2_RQSTS_ALL_CODE_RD,
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
	           __m128 spr_l2_ddrh_pi_xmm4r4(const __m128 MEM_LOAD_RETIRED_L2_HIT,
	                                        const __m128 INST_RETIRED_ANY) {
	                   
	                  register __m128 metric;
	                  metric = _mm128_div_ps(MEM_LOAD_RETIRED_L2_HIT,
	                                         INST_RETIRED_ANY);
	                  return (metric);                             
	          }
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l2_ddrh_pi_xmm4r4(const float * __restrict pMEM_LOAD_RETIRED_L2_HIT,
	                                        const float * __restrict pINST_RETIRED_ANY) {
	                   
	                  register __m128 MEM_LOAD_RETIRED_L2_HIT = 
	                                          _mm128_loadu_ps(&pMEM_LOAD_RETIRED_L2_HIT[0]);
	                  register __m128 INST_RETIRED_ANY        = 
	                                          _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(MEM_LOAD_RETIRED_L2_HIT,
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
	           __m128 spr_l2_mpi_xmm4r4(const __m128 L2_LINES_IN_ALL,
	                                    const __m128 INST_RETIRED_ANY) {
	                                    
	                  register __m128 metric;
	                  metric = _mm128_div_ps(L2_LINES_IN_ALL,
	                                         INST_RETIRED_ANY);
	                  return (metric);                           
	        }
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l2_mpi_xmm4r4(const float * __restrict pL2_LINES_IN_ALL,
	                                    const float * __restrict pINST_RETIRED_ANY) {
	                     
	                  register __m128 L2_LINES_IN_ALL  = 
	                                      _mm128_loadu_ps(&pL2_LINES_IN_ALL[0]);
	                  register __m128 INST_RETIRED_ANY =
	                                      _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);               
	                  register __m128 metric;
	                  metric = _mm128_div_ps(L2_LINES_IN_ALL,
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
	           __m128 spr_l2_ddr_mpi_xmm4r4(const __m128 MEM_LOAD_RETIRED_L2_MISS,
	                                        const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(MEM_LOAD_RETIRED_L2_MISS,
	                                         INST_RETIRED_ANY);
	                  return (metric);                               
	         }  
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l2_ddr_mpi_xmm4r4(const float * __restrict pMEM_LOAD_RETIRED_L2_MISS,
	                                        const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 MEM_LOAD_RETIRED_L2_MISS = 
	                                          _mm128_loadu_ps(&pMEM_LOAD_RETIRED_L2_MISS[0]);
	                  register __m128 INST_RETIRED_ANY         =
	                                          _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(MEM_LOAD_RETIRED_L2_MISS,
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
	           __m128 spr_l2_dc_mpi_xmm4r4(const __m128 L2_RQSTS_CODE_RD_MISS,
	                                       const __m128 INST_RETIRED_ANY) {
	                          
	                  register __m128 metric;
	                  metric = _mm128_div_ps(L2_RQSTS_CODE_RD_MISS,
	                                         INST_RETIRED_ANY);
	                  return (metric);             
	          }    
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l2_dc_mpi_xmm4r4(const float * __restrict pL2_RQSTS_CODE_RD_MISS,
	                                       const float * __restrict pINST_RETIRED_ANY) {
	                          
	                  register __m128 L2_RQSTS_CODE_RD_MISS = 
	                                          _mm128_loadu_ps(&pL2_RQSTS_CODE_RD_MISS[0]);
	                  register __m128 INST_RETIRED_ANY      =
	                                          _mm128_loadu_ps(pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(L2_RQSTS_CODE_RD_MISS,
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
	           __m128 spr_llc_drmpi_dpp_xmm4r4(const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA,
	                                           const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD,
	                                           const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF,
	                                           const __m128 INST_RETIRED_ANY) {
	                                           
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA,
	                                 _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD,
	                                               UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF));
	                  metric = _mm128_div_ps(t0,INST_RETIRED_ANY);
	                  return (metric);                               
	         } 
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_drmpi_dpp_xmm4r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA,
	                                           const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD,
	                                           const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                            
	                  restrict __m128 UNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA = 
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA[0]);
	                  restrict __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD         =
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF    =
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF[0]);
	                  register __m128 INST_RETIRED_ANY                        =
	                                             _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);              
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_LLCPREFDATA,
	                                 _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD,
	                                               UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF));
	                  metric = _mm128_div_ps(t0,INST_RETIRED_ANY);
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
	           __m128 spr_llc_crmpi_dpp_xmm4r4(const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_CRD,
	                                           const __m128 INST_RETIRED_ANY) {
	                   
	                   register __m128 metric;
	                   metric = _mm128_div_ps(UNC_CHA_TOR_INSERTS_IA_MISS_CRD,
	                                          INST_RETIRED_ANY);
	                   return (metric);                                 
	         }
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_crmpi_dpp_xmm4r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_CRD,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                   
	                   register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_CRD = 
	                                               _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_CRD[0]);
	                   register __m128 INST_RETIRED_ANY                = 
	                                               _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                   register __m128 metric;
	                   metric = _mm128_div_ps(UNC_CHA_TOR_INSERTS_IA_MISS_CRD,
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
	           __m128 spr_llc_ddrm_lat_xmm4r4(const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD,
                                                  const __m128 UNC_CHA_CLOCKTICKS,
                                                  const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const float duration_time) {
                                                 
                          const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;
                        
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
                          return (metric);
                  }    
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_ddrm_lat_xmm4r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD,
                                                  const float * __restrict pUNC_CHA_CLOCKTICKS,
                                                  const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                  const float duration_time) {
                          
                          register __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD = 
                                                     _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD[0]);
                          register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD   =
                                                     _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD[0]);
                          register __m128 UNC_CHA_CLOCKTICKS                =
                                                     _mm128_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
                          register __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD =
                                                     _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD[0]);                      
                          const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;
                        
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
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
	           __m128 spr_llc_ddr_mllr_xmm4r4(const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                          const __m128 UNC_CHA_CLOCKTICKS,
	                                          const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const float duration_time) {
	                                          
	                  const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;  
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
                          return (metric);                            
	          } 
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_ddr_mllr_xmm4r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                          const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                          const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
	                                          const float duration_time) {
	                      
	                  register __m128  UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL = 
	                                              _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL[0]);
	                  register __m128  UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL   =
	                                              _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL[0]);
	                  register __m128  UNC_CHA_CLOCKTICKS                      =
	                                              _mm128_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
	                  register __m128  UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL =
	                                              _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL[0]);                   
	                  const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;  
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_LOCAL));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
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
	           __m128 spr_llc_ddr_mlrr_xmm4r4(const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                          const __m128 UNC_CHA_CLOCKTICKS,
	                                          const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const float duration_time) {
	                     
	                                    
	                  const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;  
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
                          return (metric);  
                }
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_ddr_mlrr_xmm4r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                          const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                          const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
	                                          const float duration_time) {
	                     
	                  register __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE = 
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE   =
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE[0]);
	                  register __m128 UNC_CHA_CLOCKTICKS                       =
	                                             _mm128_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
	                  register __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE =
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE[0]);
	                  const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;  
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_REMOTE));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
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
	           __m128 spr_llc_ddr_mdlat_xmm4r4(const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR,
	                                           const __m128 UNC_CHA_CLOCKTICKS,
	                                           const __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const float duration_time) {
	                                          
	                  const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;  
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
                          return (metric);  
                }  
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_ddr_mdlat_xmm4r4(const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR,
	                                           const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                           const float * __restrict pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
	                                           const float duration_time) {
	                          
	                  register __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR = 
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR   =
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR[0]);
	                  register __m128 UNC_CHA_CLOCKTICKS                    =
	                                             _mm128_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
	                  register __m128 UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR =
	                                             _mm128_loadu_ps(&pUNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR[0]);                
	                  const __m128 C1E9 = _mm128_set1_ps(1.0e+09f); 
                          register __m128 t0;
                          register __m128 t1;
                          register __m128 vdt;
                          register __m128 metric;  
                          t0  = _mm128_mul_ps(C1E9,
                                          _mm128_div_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
                                                        UNC_CHA_TOR_INSERTS_IA_MISS_DRD_DDR)); 
                          vdt = _mm128_set1_ps(duration_time);
                          t1  = _mm128_div_ps(UNC_CHA_CLOCKTICKS,
                                          _mm128_add_ps(UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR,
                                                        UNC_CHA_TOR_OCCUPANCY_IA_MISS_DRD_DDR));
                          metric  = _mm128_mul_ps(_mm128_div_ps(t0,t1),vdt);
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
	           __m128 spr_itlb_2l_mpi_xmm4r4(const __m128 ITLB_MISSES_WALK_COMPLETED,
	                                         const __m128 INST_RETIRED_ANY) {
	                                         
	                  register __m128 metric;
	                  metric = _mm128_div_ps(ITLB_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                              
	          }
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_itlb_2l_mpi_xmm4r4(const float * __restrict pITLB_MISSES_WALK_COMPLETED,
	                                         const float * __restrict pINST_RETIRED_ANY) {
	                          
	                  register __m128 ITLB_MISSES_WALK_COMPLETED = 
	                                             _mm128_loadu_ps(&pITLB_MISSES_WALK_COMPLETED[0]);
	                  register __m128 INST_RETIRED_ANY           =
	                                             _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);               
	                  register __m128 metric;
	                  metric = _mm128_div_ps(ITLB_MISSES_WALK_COMPLETED,
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
	           __m128 spr_itlb_2l_llp_mpi_xmm4r4(const __m128 ITLB_MISSES_WALK_COMPLETED_2M_4M,
	                                             const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(ITLB_MISSES_WALK_COMPLETED_2M_4M,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }
	          
	            __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_itlb_2l_llp_mpi_xmm4r4(const float * __restrict pITLB_MISSES_WALK_COMPLETED_2M_4M,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 ITLB_MISSES_WALK_COMPLETED_2M_4M = 
	                                             _mm128_loadu_ps(&pITLB_MISSES_WALK_COMPLETED_2M_4M[0]);
	                  register __m128 INST_RETIRED_ANY                 =
	                                             _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(ITLB_MISSES_WALK_COMPLETED_2M_4M,
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
	           __m128 spr_dtlb_2l_llmpi_xmm4r4(const __m128 DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dtlb_2l_llmpi_xmm4r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 DTLB_LOAD_MISSES_WALK_COMPLETED = 
	                                           _mm128_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED[0]);
	                  register __m128 INST_RETIRED_ANY                =
	                                           _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED,
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
	           __m128 spr_dtlb_2l_2mblp_mpi_xmm4r4(const __m128 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                               const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }
	          
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dtlb_2l_2mblp_mpi_xmm4r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                               const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M = 
	                                           _mm128_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M[0]);
	                  register __m128 INST_RETIRED_ANY                      =
	                                           _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
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
	           __m128 spr_dtlb_2l_llmpi_xmm4r4(const __m128 DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }  
	          
	            __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dtlb_2l_llmpi_xmm4r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED,
	                                           const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 DTLB_LOAD_MISSES_WALK_COMPLETED = 
	                                           _mm128_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED[0]);
	                  register __m128 INST_RETIRED_ANY                =
	                                           _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED,
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
	           __m128 spr_dtlb_2l_2mllmpi_xmm4r4(const __m128 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                             const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }  
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dtlb_2l_2mllmpi_xmm4r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M = 
	                                           _mm128_loadu_ps(&pDTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M[0]);
	                  register __m128 INST_RETIRED_ANY                      =
	                                           _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
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
	           __m128 spr_dtlb_2l_store_mpi_xmm4r4(const __m128 DTLB_STORE_MISSES_WALK_COMPLETED,
	                                               const __m128 INST_RETIRED_ANY) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_STORE_MISSES_WALK_COMPLETED,
	                                         INST_RETIRED_ANY);
	                  return (metric);                                  
	          }  
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dtlb_2l_store_mpi_xmm4r4(const float * __restrict pDTLB_STORE_MISSES_WALK_COMPLETED,
	                                               const float * __restrict pINST_RETIRED_ANY) {
	                  
	                  register __m128 DTLB_STORE_MISSES_WALK_COMPLETED = 
	                                            _mm128_loadu_ps(&pDTLB_STORE_MISSES_WALK_COMPLETED[0]);
	                  register __m128 INST_RETIRED_ANY                 = 
	                                            _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(DTLB_STORE_MISSES_WALK_COMPLETED,
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
	           __m128 spr_numa_reads_ald_xmm4r4(const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL,
	                                            const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE) {
	               
	                  register __m128 t0,t1,t2;
	                  register __m128 metric;
	                  t0 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL);
	                  t1 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE);
	                  t2 = _mm128_add_ps(t0,t1);
	                  metric = _mm128_div_ps(t0,t2);
	                  return (metric);                             
	         }
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_numa_reads_ald_xmm4r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE) {
	               
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL      = 
	                                         _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL = 
	                                         _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE     =
	                                         _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE=
	                                         _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE[0]);
	                  register __m128 t0,t1,t2;
	                  register __m128 metric;
	                  t0 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL);
	                  t1 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE);
	                  t2 = _mm128_add_ps(t0,t1);
	                  metric = _mm128_div_ps(t0,t2);
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
	           __m128 spr_numa_reads_ard_xmm4r4(const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE,
	                                            const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL) {
	               
	                  register __m128 t0,t1,t2;
	                  register __m128 metric;
	                  t0 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE);
	                  t1 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL);
	                  t2 = _mm128_add_ps(t0,t1);
	                  metric = _mm128_div_ps(t0,t2);
	                  return (metric);                             
	         }   
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_numa_reads_ard_xmm4r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                            const float * __restrict pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL) {
	               
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE      = 
	                                                     _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE =
	                                                     _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL       =
	                                                     _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL  = 
	                                                     _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL[0]);
	                  register __m128 t0,t1,t2;
	                  register __m128 metric;
	                  t0 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_REMOTE,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_REMOTE);
	                  t1 = _mm128_add_ps(UNC_CHA_TOR_INSERTS_IA_MISS_DRD_LOCAL,
	                                     UNC_CHA_TOR_INSERTS_IA_MISS_DRD_PREF_LOCAL);
	                  t2 = _mm128_add_ps(t0,t1);
	                  metric = _mm128_div_ps(t0,t2);
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
	           __m128 spr_uncore_freq_xmm4r4(const __m128 UNC_CHA_CLOCKTICKS,
	                                     const float n_scks,
	                                     const float idt) { // shall be inverse value passed.
	                 
	                 const __m128 C000000001 = _mm128_set1_ps(1.0e-09f);
	                 register  __m128 vscks;
	                 register  __m128 vidt;
	                 register  __m128 t0;
	                 register  __m128 t1;
	                 register  __m128 metric;
	                 vscks = _mm128_set1_ps(n_scks);
	                 vdt   = _mm128_set1_ps(dt);
	                 t0    = _mm128_mul_ps(C000000001,
	                                   _mm128_mul_ps(vscks,UNC_CHA_CLOCKTICKS));
	                 t1    = _mm128_div_ps(UNC_CHA_CLOCKTICKS,t0);
	                 metric= _mm128_mul_ps(t1,vidt);
	                 return (metric); 
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_uncore_freq_xmm4r4(const float * __restrict pUNC_CHA_CLOCKTICKS,
	                                         const float n_scks,
	                                         const float idt) { // shall be inverse value passed.
	                 
	                 register __m128 UNC_CHA_CLOCKTICKS = _mm128_loadu_ps(&pUNC_CHA_CLOCKTICKS[0]);
	                 const __m128 C000000001 = _mm128_set1_ps(1.0e-09f);
	                 register  __m128 vscks;
	                 register  __m128 vidt;
	                 register  __m128 t0;
	                 register  __m128 t1;
	                 register  __m128 metric;
	                 vscks = _mm128_set1_ps(n_scks);
	                 vdt   = _mm128_set1_ps(dt);
	                 t0    = _mm128_mul_ps(C000000001,
	                                   _mm128_mul_ps(vscks,UNC_CHA_CLOCKTICKS));
	                 t1    = _mm128_div_ps(UNC_CHA_CLOCKTICKS,t0);
	                 metric= _mm128_mul_ps(t1,vidt);
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
	           __m128 spr_upi_data_tx_bw_xmm4r4(const __m128 UNC_UPI_TxL_FLITS_ALL_DATA,
	                                            const float idt) {
	             
	                  const __m128 C000001                     =
	                                    _mm128_set1_ps(1.0e-06f);
	                  const __m128 C711111111111111111111111   =
	                                    _mm128_set1_ps(7.11111111111111111111111f);
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(UNC_UPI_TxL_FLITS_ALL_DATA,
	                                   _mm128_mul_ps(C711111111111111111111111,
	                                                                        C000001));
	                  metric = _mm128_mul_ps(t0,vidt);
	                  return (metric);                      
	          }
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_upi_data_tx_bw_xmm4r4(const float * __restrict pUNC_UPI_TxL_FLITS_ALL_DATA,
	                                            const float idt) {
	             
	                  register __m128 UNC_UPI_TxL_FLITS_ALL_DATA = 
	                                         _mm128_loadu_ps(&pUNC_UPI_TxL_FLITS_ALL_DATA[0]);
	                  const __m128 C000001                     =
	                                    _mm128_set1_ps(1.0e-06f);
	                  const __m128 C711111111111111111111111   =
	                                    _mm128_set1_ps(7.11111111111111111111111f);
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(UNC_UPI_TxL_FLITS_ALL_DATA,
	                                   _mm128_mul_ps(C711111111111111111111111,
	                                                                        C000001));
	                  metric = _mm128_mul_ps(t0,vidt);
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
	           __m128 spr_mem_bw_rd_xmm4r4(const __m128 UNC_M_CAS_COUNT_RD,
	                                       const float idt) {
	                
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                            UNC_M_CAS_COUNT_RD,C640),
	                                                                 C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_mem_bw_rd_xmm4r4(const float * __restrict pUNC_M_CAS_COUNT_RD,
	                                       const float idt) {
	                
	                  register __m128 UNC_M_CAS_COUNT_RD = 
	                                       _mm128_loadu_ps(&pUNC_M_CAS_COUNT_RD[0]);
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                            UNC_M_CAS_COUNT_RD,C640),
	                                                                 C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
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
	           __m128 spr_mem_bw_wr_xmm4r4(const __m128 UNC_M_CAS_COUNT_WR,
	                                       const float idt) {
	                
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                            UNC_M_CAS_COUNT_WR,C640),
	                                                                 C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_mem_bw_wr_xmm4r4(const float * __restrict pUNC_M_CAS_COUNT_WR,
	                                       const float idt) {
	                
	                  register __m128 UNC_M_CAS_COUNT_WR = 
	                                       _mm128_loadu_ps(&pUNC_M_CAS_COUNT_WR[0]);
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                            UNC_M_CAS_COUNT_WR,C640),
	                                                                 C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
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
	           __m128 spr_mem_bw_total_xmm4r4(const __m128 UNC_M_CAS_COUNT_RD,
	                                          const __m128 UNC_M_CAS_COUNT_WR,
	                                          const float idt) {
	                  
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;  
	                  register __m128 t1;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0 =  _mm128_add_ps(UNC_M_CAS_COUNT_RD,
	                                      UNC_M_CAS_COUNT_WR);
	                  t1 =  _mm128_mul_ps(
	                                _mm128_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm128_mul_ps(t1,vidt);
	                  return (metric);                                                
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_mem_bw_total_xmm4r4(const float * __restrict pUNC_M_CAS_COUNT_RD,
	                                          const float * __restrict pUNC_M_CAS_COUNT_WR,
	                                          const float idt) {
	                  
	                  register __m128 UNC_M_CAS_COUNT_RD = 
	                                       _mm128_loadu_ps(&pUNC_M_CAS_COUNT_RD[0]);
	                  register __m128 UNC_M_CAS_COUNT_WR =
	                                       _mm128_loadu_ps(&pUNC_M_CAS_COUNT_WR[0]);
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;  
	                  register __m128 t1;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0 =  _mm128_add_ps(UNC_M_CAS_COUNT_RD,
	                                      UNC_M_CAS_COUNT_WR);
	                  t1 =  _mm128_mul_ps(
	                                _mm128_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm128_mul_ps(t1,vidt);
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
	           __m128 spr_io_bw_read_xmm4r4(const __m128 UNC_CHA_TOR_INSERTS_IO_PCIRDCUR,
	                                       const float idt) {
	                
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                         UNC_CHA_TOR_INSERTS_IO_PCIRDCUR,C640),
	                                                                        C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_io_bw_read_xmm4r4(const float * __restrict pUNC_CHA_TOR_INSERTS_IO_PCIRDCUR,
	                                        const float idt) {
	                  
	                  register __m128 UNC_CHA_TOR_INSERTS_IO_PCIRDCUR = 
	                                    _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IO_PCIRDCUR[0]);
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                         UNC_CHA_TOR_INSERTS_IO_PCIRDCUR,C640),
	                                                                        C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
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
	           __m128 spr_io_bw_read_xmm4r4(  const __m128 UNC_CHA_TOR_INSERTS_IO_ITOM,
	                                          const __m128 UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR ,
	                                          const float idt) {
	                  
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;  
	                  register __m128 t1;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0 =  _mm128_add_ps(UNC_CHA_TOR_INSERTS_IO_ITOM,
	                                      UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR);
	                  t1 =  _mm128_mul_ps(
	                                _mm128_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm128_mul_ps(t1,vidt);
	                  return (metric);                                                
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_io_bw_read_xmm4r4(  const float * __restrict pUNC_CHA_TOR_INSERTS_IO_ITOM,
	                                          const float * __restrict pUNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR ,
	                                          const float idt) {
	                  
	                  register __m128 UNC_CHA_TOR_INSERTS_IO_ITOM          = 
	                                         _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IO_ITOM[0]);
	                  register __m128 UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR =
	                                         _mm128_loadu_ps(&pUNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR[0]);
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;  
	                  register __m128 t1;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0 =  _mm128_add_ps(UNC_CHA_TOR_INSERTS_IO_ITOM,
	                                      UNC_CHA_TOR_INSERTS_IO_ITOMCACHENEAR);
	                  t1 =  _mm128_mul_ps(
	                                _mm128_mul_ps(t0,vidt),
	                                                 C000001);
	                  metric = _mm128_mul_ps(t1,vidt);
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
	           __m128 spr_uops_ddic_xmm4r4(const __m128  IDQ_DSB_UOPS,
	                                       const __m128  IDQ_MITE_UOPS,
	                                       const __m128  IDQ_MS_UOPS,
	                                       const __m128  LSD_UOPS) {
	             
	                  register __m128 t0;
	               	  register __m128 metric;
	                  t0 = _mm128_add_ps(_mm128_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm128_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm128_div_ps(IDQ_DSB_UOPS,t0);
	                  return (metric);                             
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_uops_ddic_xmm4r4(const float * __restrict  pIDQ_DSB_UOPS,
	                                       const float * __restrict  pIDQ_MITE_UOPS,
	                                       const float * __restrict  pIDQ_MS_UOPS,
	                                       const float * __restrict  pLSD_UOPS) {
	                  
	                  register __m128 IDQ_DSB_UOPS  = 
	                                     _mm128_loadu_ps(&pIDQ_DSB_UOPS[0]);
	                  register __m128 IDQ_MITE_UOPS =
	                                     _mm128_loadu_ps(&IDQ_MITE_UOPS[0]);
	                  register __m128 IDQ_MS_UOPS   =
	                                     _mm128_loadu_ps(&IDQ_MS_UOPS[0]);
	                  register __m128 LSD_UOPS      = 
	                                     _mm128_loadu_ps(&LSD_UOPS[0]);
	                  register __m128 t0;
	               	  register __m128 metric;
	                  t0 = _mm128_add_ps(_mm128_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm128_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm128_div_ps(IDQ_DSB_UOPS,t0);
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
	           __m128 spr_uops_dldp_xmm4r4(const __m128  IDQ_DSB_UOPS,
	                                       const __m128  IDQ_MITE_UOPS,
	                                       const __m128  IDQ_MS_UOPS,
	                                       const __m128  LSD_UOPS) {
	             
	                  register __m128 t0;
	               	  register __m128 metric;
	                  t0 = _mm128_add_ps(_mm128_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm128_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm128_div_ps(IDQ_MITE_UOPS,t0);
	                  return (metric);                             
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_uops_dldp_xmm4r4(const float * __restrict  pIDQ_DSB_UOPS,
	                                       const float * __restrict  pIDQ_MITE_UOPS,
	                                       const float * __restrict  pIDQ_MS_UOPS,
	                                       const float * __restrict  pLSD_UOPS) {
	             
	                  register __m128 IDQ_DSB_UOPS  = 
	                                     _mm128_loadu_ps(&pIDQ_DSB_UOPS[0]);
	                  register __m128 IDQ_MITE_UOPS =
	                                     _mm128_loadu_ps(&IDQ_MITE_UOPS[0]);
	                  register __m128 IDQ_MS_UOPS   =
	                                     _mm128_loadu_ps(&IDQ_MS_UOPS[0]);
	                  register __m128 LSD_UOPS      = 
	                                     _mm128_loadu_ps(&LSD_UOPS[0]);
	                  register __m128 t0;
	               	  register __m128 metric;
	                  t0 = _mm128_add_ps(_mm128_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm128_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm128_div_ps(IDQ_MITE_UOPS,t0);
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
	           __m128 spr_uops_dmseq_xmm4r4(const __m128  IDQ_DSB_UOPS,
	                                       const __m128  IDQ_MITE_UOPS,
	                                       const __m128  IDQ_MS_UOPS,
	                                       const __m128  LSD_UOPS) {
	             
	                  register __m128 t0;
	               	  register __m128 metric;
	                  t0 = _mm128_add_ps(_mm128_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm128_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm128_div_ps(IDQ_MS_UOPS,t0);
	                  return (metric);                             
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_uops_dmseq_xmm4r4(const float * __restrict  pIDQ_DSB_UOPS,
	                                       const  float * __restrict  pIDQ_MITE_UOPS,
	                                       const  float * __restrict  pIDQ_MS_UOPS,
	                                       const  float * __restrict  pLSD_UOPS) {
	             
	                  register __m128 IDQ_DSB_UOPS  = 
	                                     _mm128_loadu_ps(&pIDQ_DSB_UOPS[0]);
	                  register __m128 IDQ_MITE_UOPS =
	                                     _mm128_loadu_ps(&IDQ_MITE_UOPS[0]);
	                  register __m128 IDQ_MS_UOPS   =
	                                     _mm128_loadu_ps(&IDQ_MS_UOPS[0]);
	                  register __m128 LSD_UOPS      = 
	                                     _mm128_loadu_ps(&LSD_UOPS[0]);
	                  register __m128 t0;
	               	  register __m128 metric;
	                  t0 = _mm128_add_ps(_mm128_add_ps(IDQ_DSB_UOPS,
	                                                   IDQ_MITE_UOPS),
	                                     _mm128_add_ps(IDQ_MS_UOPS,
	                                                   LSD_UOPS));
	                  metric = _mm128_div_ps(IDQ_MS_UOPS,t0);
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
	           __m128 spr_llc_mlmbw_rd_xmm4r4(const __m128 UNC_CHA_REQUESTS_READS_LOCAL,
	                                          const float idt) {
	                
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_LOCAL,C640),
	                                                                        C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_mlmbw_rd_xmm4r4(const float * __restrict pUNC_CHA_REQUESTS_READS_LOCAL,
	                                          const float idt) {
	                
	                  register __m128 UNC_CHA_REQUESTS_READS_LOCAL = 
	                                         _mm128_loadu_ps(&pUNC_CHA_REQUESTS_READS_LOCAL[0]);
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_LOCAL,C640),
	                                                                        C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
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
	           __m128 spr_llc_mrmbw_rd_xmm4r4(const __m128 UNC_CHA_REQUESTS_READS_REMOTE,
	                                          const float idt) {
	                
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_REMOTE,C640),
	                                                                        C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
	                  return (metric);
	          }  
	          
	            __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_llc_mrmbw_rd_xmm4r4(const float * __restrict pUNC_CHA_REQUESTS_READS_REMOTE,
	                                          const float idt) {
	                
	                  register __m128 UNC_CHA_REQUESTS_READS_REMOTE = 
	                                         _mm128_loadu_ps(&pUNC_CHA_REQUESTS_READS_REMOTE[0]);
	                  const __m128 C640    = _mm128_mul_ps(64.0f);
	                  const __m128 C000001 = _mm128_set1_ps(1.0e-06f);              
	                  register __m128 vidt;
	                  register __m128 t0;
	                  register __m128 metric;
	                  vidt = _mm128_set1_ps(idt);
	                  t0   = _mm128_mul_ps(_mm128_mul_ps(
	                                         UNC_CHA_REQUESTS_READS_REMOTE,C640),
	                                                                        C000001);
	                  metric = _mm128_mul_ps(t0,vidt);
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
	           __m128 spr_icache_misses_xmm4r4(const __m128 ICACHE_DATA_STALLS,
	                                           const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(ICACHE_DATA_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
	                  return (metric);                                 
	         }   
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_icache_misses_xmm4r4(const float * __restrict pICACHE_DATA_STALLS,
	                                           const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m128 ICACHE_DATA_STALLS = 
	                                        _mm128_loadu_ps(&pICACHE_DATA_STALLS[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD = 
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(ICACHE_DATA_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
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
	           __m128 spr_itlb_misses_xmm4r4(const __m128 ICACHE_TAG_STALLS,
	                                         const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(ICACHE_TAG_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
	                  return (metric);                                 
	         }   
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_itlb_misses_xmm4r4(const float * __restrict pICACHE_TAG_STALLS,
	                                         const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m128 ICACHE_TAG_STALLS = 
	                                        _mm128_loadu_ps(&pICACHE_TAG_STALLS[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD = 
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(ICACHE_TAG_STALLS,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
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
	           __m128 spr_branch_resteers_xmm4r4(const __m128 INT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const __m128 CPU_CLK_UNHALTED_THREAD,
	                                             const __m128 INT_MISC_UNKNOWN_BRANCH_CYCLES) {
	                  
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 metric;  
	                  t0  = _mm128_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1  = _mm128_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,
	                                     _mm128_add_ps(t0,t1));
	                  return (metric);                                
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_branch_resteers_xmm4r4(const float * __restrict pINT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                             const float * __restrict pINT_MISC_UNKNOWN_BRANCH_CYCLES) {
	                  
	                  register __m128 INT_MISC_CLEAR_RESTEER_CYCLES = 
	                                          _mm128_loadu_ps(&pINT_MISC_CLEAR_RESTEER_CYCLES[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  register __m128 INT_MISC_UNKNOWN_BRANCH_CYCLES=
	                                          _mm128_loadu_ps(&pINT_MISC_UNKNOWN_BRANCH_CYCLES[0]);
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 metric;  
	                  t0  = _mm128_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1  = _mm128_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                      CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,
	                                     _mm128_add_ps(t0,t1));
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
	           __m128 spr_unknown_branches_xmm4r4(const __m128 INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                              const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
	                  return (metric);
	                                                      
	         }   
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_unknown_branches_xmm4r4(const float * __restrict pINT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                              const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m128 INT_MISC_UNKNOWN_BRANCH_CYCLES = 
	                                       _mm128_loadu_ps(&pINT_MISC_UNKNOWN_BRANCH_CYCLES[0]);
	                   register __m128 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
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
	           __m128 spr_dsb_switches_xmm4r4(const __m128 DSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                          const __m128 CPU_CLK_UNHALTED_THREAD) {
	                 
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;  
	                  t0 = _mm128_div_ps(DSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
	                  return (metric);                              
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dsb_switches_xmm4r4(const float * __restrict pDSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                          const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                  register __m128 DSB2MITE_SWITCHES_PENALTY_CYCLES = 
	                                          _mm128_loadu_ps(&pDSB2MITE_SWITCHES_PENALTY_CYCLES[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;  
	                  t0 = _mm128_div_ps(DSB2MITE_SWITCHES_PENALTY_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
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
	           __m128 spr_lcp_xmm4r4(const __m128 DECODE_LCP,
	                                 const __m128 CPU_CLK_UNHALTED_THREAD) {
	                 
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;  
	                  t0 = _mm128_div_ps(DECODE_LCP,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
	                  return (metric);                            
	         }    
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_lcp_xmm4r4(const float * __restrict pDECODE_LCP,
	                                 const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                  register __m128 DECODE_LCP              =
	                                  _mm128_loadu_ps(&pDECODE_LCP[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD =
	                                          _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C1000 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;  
	                  t0 = _mm128_div_ps(DECODE_LCP,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C1000,t0);
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
	           __m128 spr_info_mem_l2mpki_xmm4r4(const __m128 MEM_LOAD_RETIRED_L2_MISS,
	                                             const __m128 INST_RETIRED_ANY) {
	                                             
	                  const __m128 C10000 = _mm128_set1_ps(1000.0f);
	                  register __m128 t0;
	                  register __m128 metric; 
	                  t0 = _mm128_div_ps(MEM_LOAD_RETIRED_L2_MISS,
	                                     INST_RETIRED_ANY);
	                  metric = _mm128_mul_ps(C10000,t0);
	                  return (metric);                                                       
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_info_mem_l2mpki_xmm4r4(const float * __restrict pMEM_LOAD_RETIRED_L2_MISS,
	                                             const float * __restrict pINST_RETIRED_ANY) {
	                                  
	                  register __m128 MEM_LOAD_RETIRED_L2_MISS = 
	                                          _mm128_loadu_ps(&pMEM_LOAD_RETIRED_L2_MISS[0]);
	                  register __m128 INST_RETIRED_ANY         =
	                                          _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);          
	                  const __m128 C10000 = _mm128_set1_ps(1000.0f);
	                  register __m128 t0;
	                  register __m128 metric; 
	                  t0 = _mm128_div_ps(MEM_LOAD_RETIRED_L2_MISS,
	                                     INST_RETIRED_ANY);
	                  metric = _mm128_mul_ps(C10000,t0);
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
	           __m128 spr_bad_spec_ipm_xmm4r4(const __m128 INST_RETIRED_ANY,
	                                          const __m128 BR_MISP_RETIRED_ALL_BRANCHES) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(INST_RETIRED_ANY,
	                                     BR_MISP_RETIRED_ALL_BRANCHES);
	                  return (metric);                          
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_bad_spec_ipm_xmm4r4(const float * __restrict pINST_RETIRED_ANY,
	                                          const float * __restrict pBR_MISP_RETIRED_ALL_BRANCHES) {
	                  
	                  register __m128 INST_RETIRED_ANY = 
	                                      _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128  BR_MISP_RETIRED_ALL_BRANCHES = 
	                                      _mm128_loadu_ps(&pBR_MISP_RETIRED_ALL_BRANCHES[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(INST_RETIRED_ANY,
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
	           __m128 spr_inst_mix_iptb_xmm4r4(const __m128 INST_RETIRED_ANY,
	                                           const __m128 BR_INST_RETIRED_NEAR_TAKEN) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(INST_RETIRED_ANY,
	                                         BR_INST_RETIRED_NEAR_TAKEN);
	                  return (metric);                          
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_inst_mix_iptb_xmm4r4(const float * __restrict pINST_RETIRED_ANY,
	                                           const float * __restrict pBR_INST_RETIRED_NEAR_TAKEN) {
	                  
	                  register __m128 INST_RETIRED_ANY = 
	                                      _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128  BR_INST_RETIRED_NEAR_TAKEN = 
	                                      _mm128_loadu_ps(&pBR_INST_RETIRED_NEAR_TAKEN[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(INST_RETIRED_ANY,
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
	           __m128 spr_core_coreipc_xmm4r4( const __m128 INST_RETIRED_ANY,
	                                           const __m128 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                  register __m128 metric;
	                  metric = _mm128_div_ps(INST_RETIRED_ANY,
	                                         CPU_CLK_UNHALTED_DISTRIBUTED);
	                  return (metric);                          
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_core_coreipc_xmm4r4( const float * __restrict pINST_RETIRED_ANY,
	                                           const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                  register __m128 INST_RETIRED_ANY = 
	                                      _mm128_loadu_ps(&pINST_RETIRED_ANY[0]);
	                  register __m128  CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                      _mm128_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                  register __m128 metric;
	                  metric = _mm128_div_ps(INST_RETIRED_ANY,
	                                         CPU_CLK_UNHALTED_DISTRIBUTED);
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
	           __m128 spr_avx_assists_xmm4r4(const __m128 ASSISTS_SSE_AVX_MIX,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C63  = _mm128_set1_ps(63.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_mul_ps(C63,
	                             _mm128_div_ps(ASSISTS_SSE_AVX_MIX,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_avx_assists_xmm4r4(const float * __restrict pASSISTS_SSE_AVX_MIX,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m128 ASSISTS_SSE_AVX_MIX = 
	                                         _mm128_loadu_ps(&pASSISTS_SSE_AVX_MIX[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C63  = _mm128_set1_ps(63.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_mul_ps(C63,
	                             _mm128_div_ps(ASSISTS_SSE_AVX_MIX,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
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
	           __m128 spr_fp_assists_xmm4r4(const __m128 ASSISTS_FP,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C30  = _mm128_set1_ps(30.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_mul_ps(C30,
	                             _mm128_div_ps(ASSISTS_FP,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_fp_assists_xmm4r4(const float * __restrict pASSISTS_FP,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m128 ASSISTS_FP = 
	                                         _mm128_loadu_ps(&pASSISTS_FP[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C30  = _mm128_set1_ps(30.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_mul_ps(C30,
	                             _mm128_div_ps(ASSISTS_FP,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
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
	           __m128 spr_page_faults_xmm4r4(const __m128 ASSISTS_PAGE_FAULTS,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C99  = _mm128_set1_ps(99.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_mul_ps(C99,
	                             _mm128_div_ps(ASSISTS_PAGE_FAULTS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_page_faults_xmm4r4(const float * __restrict pASSISTS_PAGE_FAULTS,
	                                         const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m128 ASSISTS_PAGE_FAULTS = 
	                                         _mm128_loadu_ps(&pASSISTS_PAGE_FAULTS[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C99  = _mm128_set1_ps(99.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_mul_ps(C99,
	                             _mm128_div_ps(ASSISTS_PAGE_FAULTS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
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
	           __m128 spr_micrcode_seq_xmm4r4(const __m128 UOPS_RETIRED_MS,
	                                          const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_div_ps(UOPS_RETIRED_MS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
	                  return (metric);                               
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_micrcode_seq_xmm4r4(const float * __restrict pUOPS_RETIRED_MS,
	                                          const float TOPDOWN_SLOTS_perf_metrics) {
	                  
	                  register __m128 UOPS_RETIRED_MS = 
	                                      _mm128_loadu_ps(&pUOPS_RETIRED_MS[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 vx;
	                  register __m128 t0;
	                  register __metric;
	                  vx = _mm128_set1_ps(TOPDOWN_SLOTS_perf_metrics);
	                  t0 = _mm128_div_ps(UOPS_RETIRED_MS,
	                                           TOPDOWN_SLOTS_perf_metrics));
	                  metric = _mm128_mul_ps(C100,t0);
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
	           __m128 spr_branch_resteers_xmm4r4(const __m128 INT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const __m128 CPU_CLK_UNHALTED_THREAD,
	                                             const __m128 INT_MISC_UNKNOWN_BRANCH_CYCLES) {
	           
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0,t1;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C100,
	                                     _mm128_add_ps(t0,t1));
	                  return (metric);                                  
	          }  
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_branch_resteers_xmm4r4(const float * __restrict pINT_MISC_CLEAR_RESTEER_CYCLES,
	                                             const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                             const float * __restrict pINT_MISC_UNKNOWN_BRANCH_CYCLES) {
	            
	                  register __m128 INT_MISC_CLEAR_RESTEER_CYCLES = 
	                                          _mm128_loadu_ps(&pINT_MISC_CLEAR_RESTEER_CYCLES[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD       =
	                                          _mm128_loadu_ps(&pINT_MISC_CLEAR_RESTEER_CYCLES[0]);
	                  register __m128 INT_MISC_UNKNOWN_BRANCH_CYCLES=
	                                          _mm128_loadu_ps(&pINT_MISC_UNKNOWN_BRANCH_CYCLES[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0,t1;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(INT_MISC_CLEAR_RESTEER_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_div_ps(INT_MISC_UNKNOWN_BRANCH_CYCLES,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C100,
	                                     _mm128_add_ps(t0,t1));
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
	           __m128 spr_mite_xmm4r4(const __m128 IDQ_MITE_CYCLES_ANY,
	                                  const __m128 IDQ_MITE_CYCLES_OK,
	                                  const __m128 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                    const __m128 C100 = _mm128_set1_ps(100.0f);
	                    const __m128 C05  = _mm128_set1_ps(0.5f);
	                    register __m128 t0,t1;
	                    register __m128 metric;
	                    t0 = _mm128_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm128_sub_ps(IDQ_MITE_CYCLES_ANY,
	                                       IDQ_MITE_CYCLES_OK);
	                    metric = _mm128_mul_ps(C100,
	                                   _mm128_div_ps(t1,t0));
	                    return (metric);                        
	         }  
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_mite_xmm4r4(const float * __restrict pIDQ_MITE_CYCLES_ANY,
	                                  const float * __restrict pIDQ_MITE_CYCLES_OK,
	                                  const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                  
	                    register __m128 IDQ_MITE_CYCLES_ANY  =
	                                            _mm128_loadu_ps(&pIDQ_MITE_CYCLES_ANY[0]);
	                    register __m128 IDQ_MITE_CYCLES_OK   =
	                                            _mm128_loadu_ps(&pIDQ_MITE_CYCLES_OK[0]);
	                    register __m128 CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                            _mm128_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                    const __m128 C100 = _mm128_set1_ps(100.0f);
	                    const __m128 C05  = _mm128_set1_ps(0.5f);
	                    register __m128 t0,t1;
	                    register __m128 metric;
	                    t0 = _mm128_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm128_sub_ps(IDQ_MITE_CYCLES_ANY,
	                                       IDQ_MITE_CYCLES_OK);
	                    metric = _mm128_mul_ps(C100,
	                                   _mm128_div_ps(t1,t0));
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
	           __m128 spr_decoder0_only_xmm4r4(const __m128 INST_DECODED_DECODERS_c1,
	                                           const __m128 INST_DECODED_DECODERS_c2,
	                                           const __m128 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    const __m128 C100 = _mm128_set1_ps(100.0f);
	                    const __m128 C05  = _mm128_set1_ps(0.5f);
	                    register __m128 t0,t1;
	                    register __m128 metric;
	                    t0 = _mm128_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm128_sub_ps(INST_DECODED_DECODERS_c1,
	                                       INST_DECODED_DECODERS_c2);
	                    metric = _mm128_mul_ps(C100,
	                                   _mm128_div_ps(t1,t0));
	                    return (metric);                                        
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_decoder0_only_xmm4r4(const float * __restrict pINST_DECODED_DECODERS_c1,
	                                           const float * __restrict pINST_DECODED_DECODERS_c2,
	                                           const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    register __m128 INST_DECODED_DECODERS_c1 = 
	                                        _mm128_loadu_ps(&pINST_DECODED_DECODERS_c1[0]);
	                    register __m128 INST_DECODED_DECODERS_c2 = 
	                                        _mm128_loadu_ps(&pINST_DECODED_DECODERS_c2[0]);
	                    register __m128 CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                    const __m128 C100 = _mm128_set1_ps(100.0f);
	                    const __m128 C05  = _mm128_set1_ps(0.5f);
	                    register __m128 t0,t1;
	                    register __m128 metric;
	                    t0 = _mm128_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm128_sub_ps(INST_DECODED_DECODERS_c1,
	                                       INST_DECODED_DECODERS_c2);
	                    metric = _mm128_mul_ps(C100,
	                                   _mm128_div_ps(t1,t0));
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
	           __m128 spr_dsb_xmm4r4(          const __m128 IDQ_DSB_CYCLES_ANY,
	                                           const __m128 IDQ_DSB_CYCLES_OK,
	                                           const __m128 CPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    const __m128 C100 = _mm128_set1_ps(100.0f);
	                    const __m128 C05  = _mm128_set1_ps(0.5f);
	                    register __m128 t0,t1;
	                    register __m128 metric;
	                    t0 = _mm128_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm128_sub_ps(IDQ_DSB_CYCLES_ANY,
	                                       IDQ_DSB_CYCLES_OK);
	                    metric = _mm128_mul_ps(C100,
	                                   _mm128_div_ps(t1,t0));
	                    return (metric);                                        
	        }
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dsb_only_xmm4r4(     const float * __restrict pIDQ.DSB_CYCLES_ANY,
	                                           const float * __restrict pIDQ.DSB_CYCLES_OK,
	                                           const float * __restrict pCPU_CLK_UNHALTED_DISTRIBUTED) {
	                                           
	                    register __m128 IDQ_DSB_CYCLES_ANY = 
	                                        _mm128_loadu_ps(&pIDQ_DSB_CYCLES_ANY[0]);
	                    register __m128 IDQ_DSB_CYCLES_OK = 
	                                        _mm128_loadu_ps(&pIDQ_DSB_CYCLES_OK[0]);
	                    register __m128 CPU_CLK_UNHALTED_DISTRIBUTED = 
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_DISTRIBUTED[0]);
	                    const __m128 C100 = _mm128_set1_ps(100.0f);
	                    const __m128 C05  = _mm128_set1_ps(0.5f);
	                    register __m128 t0,t1;
	                    register __m128 metric;
	                    t0 = _mm128_mul_ps(CPU_CLK_UNHALTED_DISTRIBUTED,C05);
	                    t1 = _mm128_sub_ps(IDQ.DSB_CYCLES_ANY,
	                                       IDQ.DSB_CYCLES_OK);
	                    metric = _mm128_mul_ps(C100,
	                                   _mm128_div_ps(t1,t0));
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
	           __m128 spr_l1_bound_xmm4r4(const __m128 EXE_ACTIVITY_BOUND_ON_LOADS,
	                                      const __m128 MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                      const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0,t1;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(_mm128_sub_ps(EXE_ACTIVITY_BOUND_ON_LOADS,
	                                                   MEMORY_ACTIVITY_STALLS_L1D_MISS),
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_max_ps(t0,_mm128_setzero_ps());
	                  metric = _mm128_mul_ps(C100,t1);
	                  return (metric);                           
	          }
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_l1_bound_xmm4r4(const float * __restrict pEXE_ACTIVITY_BOUND_ON_LOADS,
	                                      const float * __restrict pMEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                      const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m128 EXE_ACTIVITY_BOUND_ON_LOADS = 
	                                              _mm128_loadu_ps(&pEXE_ACTIVITY_BOUND_ON_LOADS[0]);
	                  register __m128 MEMORY_ACTIVITY_STALLS_L1D_MISS =
	                                              _mm128_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L1D_MISS[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD     = 
	                                              _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0,t1;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(_mm128_sub_ps(EXE_ACTIVITY_BOUND_ON_LOADS,
	                                                   MEMORY_ACTIVITY_STALLS_L1D_MISS),
	                                      CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_max_ps(t0,_mm128_setzero_ps());
	                  metric = _mm128_mul_ps(C100,t1);
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
	           __m128 spr_dtlb_load_xmm4r4(const __m128 DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                       const __m128 DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                       const __m128 CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                       const __m128 MEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                       const __m128 CPU_CLK_UNHALTED_THREAD) {
	                            
	                  const __m128 C7 = _mm128_set1_ps(7.0f);
	                  const __m128 C0 = _mm128_setzero_ps();
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 t2;
	                  register __m128 metric;
	                  t0 = _mm128_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t1 = _mm128_max_ps(_mm128_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm128_min_ps(t0,t1);
	                  metric = _mm128_mul_ps(C100,
	                                     _mm128_div_ps(t2,CPU_CLK_UNHALTED_THREAD));
	                  return (metric);      
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_dtlb_load_xmm4r4(const float * __restrict pDTLB_LOAD_MISSES_STLB_HIT_c1,
	                                       const float * __restrict pDTLB_LOAD_MISSES_WALK_ACTIVE,
	                                       const float * __restrict pCYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                       const float * __restrict pMEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                       const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                            
	                  const __m128 DTLB_LOAD_MISSES_STLB_HIT_c1 = 
	                                        _mm128_loadu_ps(&pDTLB_LOAD_MISSES_STLB_HIT_c1[0]);
	                  const __m128 DTLB_LOAD_MISSES_WALK_ACTIVE =
	                                        _mm128_loadu_ps(&pDTLB_LOAD_MISSES_WALK_ACTIVE[0]);
	                  const __m128 CYCLE_ACTIVITY_CYCLES_MEM_ANY =
	                                        _mm128_loadu_ps(&pCYCLE_ACTIVITY_CYCLES_MEM_ANY[0]);
	                  const __m128 MEMORY_ACTIVITY_CYCLES_L1D_MISS=
	                                        _mm128_loadu_ps(&pMEMORY_ACTIVITY_CYCLES_L1D_MISS[0]);
	                  const __m128 CPU_CLK_UNHALTED_THREAD        =
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C7 = _mm128_set1_ps(7.0f);
	                  const __m128 C0 = _mm128_setzero_ps();
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 t2;
	                  register __m128 metric;
	                  t0 = _mm128_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t1 = _mm128_max_ps(_mm128_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm128_min_ps(t0,t1);
	                  metric = _mm128_mul_ps(C100,
	                                     _mm128_div_ps(t2,CPU_CLK_UNHALTED_THREAD));
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
	           __m128 spr_load_stlb_hit_xmm4r4(const __m128 DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                          const __m128 DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                          const __m128 CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                          const __m128 MEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                          const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 C7 = _mm128_set1_ps(7.0f);
	                  const __m128 C0 = _mm128_setzero_ps();
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 t2;
	                  register __m128 t3
	                  register __m128 metric;  
	                  t0 = _mm128_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t3 = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_max_ps(_mm128_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm128_div_ps(_mm128_min_ps(t0,t1),
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C100,
	                                     _mm128_sub_ps(t2,t3));
	                  return (metric);
	                                                    
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_load_stlb_hit_xmm4r4(const float * __restrict pDTLB_LOAD_MISSES_STLB_HIT_c1,
	                                          const float * __restrict pDTLB_LOAD_MISSES_WALK_ACTIVE,
	                                          const float * __restrict pCYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                          const float * __restrict pMEMORY_ACTIVITY_CYCLES_L1D_MISS,
	                                          const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 DTLB_LOAD_MISSES_STLB_HIT_c1 = 
	                                        _mm128_loadu_ps(&pDTLB_LOAD_MISSES_STLB_HIT_c1[0]);
	                  const __m128 DTLB_LOAD_MISSES_WALK_ACTIVE =
	                                        _mm128_loadu_ps(&pDTLB_LOAD_MISSES_WALK_ACTIVE[0]);
	                  const __m128 CYCLE_ACTIVITY_CYCLES_MEM_ANY =
	                                        _mm128_loadu_ps(&pCYCLE_ACTIVITY_CYCLES_MEM_ANY[0]);
	                  const __m128 MEMORY_ACTIVITY_CYCLES_L1D_MISS=
	                                        _mm128_loadu_ps(&pMEMORY_ACTIVITY_CYCLES_L1D_MISS[0]);
	                  const __m128 CPU_CLK_UNHALTED_THREAD        =
	                                        _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C7 = _mm128_set1_ps(7.0f);
	                  const __m128 C0 = _mm128_setzero_ps();
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 t2;
	                  register __m128 t3
	                  register __m128 metric;  
	                  t0 = _mm128_fmad_ps(C7,DTLB_LOAD_MISSES_STLB_HIT_c1,
	                                      DTLB_LOAD_MISSES_WALK_ACTIVE);
	                  t3 = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_max_ps(_mm128_sub_ps(CYCLE_ACTIVITY_CYCLES_MEM_ANY,
	                                                   MEMORY_ACTIVITY_CYCLES_L1D_MISS),C0); 
	                  t2 = _mm128_div_ps(_mm128_min_ps(t0,t1),
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C100,
	                                     _mm128_sub_ps(t2,t3));
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
	           __m128 spr_load_stlb_miss_xmm4r4(const __m128 DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                            const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C100,t0);
	                  return (metric);                                 
	         } 
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_load_stlb_miss_xmm4r4(const float * __restrict pDTLB_LOAD_MISSES_WALK_ACTIVE,
	                                            const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m128 DTLB_LOAD_MISSES_WALK_ACTIVE = 
	                                           _mm128_loadu_ps(&pDTLB_LOAD_MISSES_WALK_ACTIVE[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD      =
	                                           _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  register __m128 t0;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(DTLB_LOAD_MISSES_WALK_ACTIVE,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  metric = _mm128_mul_ps(C100,t0);
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
	           __m128 store_fwd_blk_xmm4r4(const __m128 LD_BLOCKS_STORE_FORWARD,
	                                       const __m128 CPU_CLK_UNHALTED_THREAD) {
	                                       
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C13  = _mm128_set1_ps(13.0f);
	                  const __m128 C1   = _mm128_set1_ps(1.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 metric;                   
	                  t0 = _mm128_div_ps(LD_BLOCKS_STORE_FORWARD,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_min_ps(_mm128_mul_ps(C13,t0),C1);
	                  metric = _mm128_mul_ps(C100,t1);
	                  return (metric);
	         }
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 store_fwd_blk_xmm4r4(const float * __restrict pLD_BLOCKS_STORE_FORWARD,
	                                       const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                                       
	                  register __m128 LD_BLOCKS_STORE_FORWARD = 
	                                           _mm128_loadu_ps(&pLD_BLOCKS_STORE_FORWARD[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD =
	                                           _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C13  = _mm128_set1_ps(13.0f);
	                  const __m128 C1   = _mm128_set1_ps(1.0f);
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 metric;                   
	                  t0 = _mm128_div_ps(LD_BLOCKS_STORE_FORWARD,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t1 = _mm128_min_ps(_mm128_mul_ps(C13,t0),C1);
	                  metric = _mm128_mul_ps(C100,t1);
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
	           __m128 spr_split_loads_xmm4r4(const __m128 L1D_PEND_MISS_PENDING,
	                                         const __m128 MEM_LOAD_COMPLETED_L1_MISS_ANY,
	                                         const __m128 LD_BLOCKS_NO_SR,
	                                         const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C1   = _mm128_set1_ps(1.0f);  
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 t2;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(L1_PEND_MISS_PENDING,
	                                     MEM_LOAD_COMPLETED_L1_MISS_ANY);
	                  t1 = _mm128_div_ps(LD_BLOCKS_NO_SR,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t2 = _mm128_min_ps(_mm128_mul_ps(t0,t1),C1);
	                  metric = _mm128_mul_ps(C100,t2);
	                  return (metric);                          
	         }    
	         
	              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 spr_split_loads_xmm4r4(const float * __restrict pL1D_PEND_MISS_PENDING,
	                                         const float * __restrict pMEM_LOAD_COMPLETED_L1_MISS_ANY,
	                                         const float * __restrict pLD_BLOCKS_NO_SR,
	                                         const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                  register __m128 L1D_PEND_MISS_PENDING = 
	                                          _mm128_loadu_ps(&pL1D_PEND_MISS_PENDING[0]);
	                  register __m128 MEM_LOAD_COMPLETED_L1_MISS_ANY = 
	                                          _mm128_loadu_ps(&pMEM_LOAD_COMPLETED_L1_MISS_ANY[0]);
	                  register __m128 LD_BLOCKS_NO_SR        = 
	                                          _mm128_loadu_ps(&pLD_BLOCKS_NO_SR[0]);
	                  register __m128 CPU_CLK_UNHALTED_THREAD=
	                                          _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                  const __m128 C100 = _mm128_set1_ps(100.0f);
	                  const __m128 C1   = _mm128_set1_ps(1.0f);  
	                  register __m128 t0;
	                  register __m128 t1;
	                  register __m128 t2;
	                  register __m128 metric;
	                  t0 = _mm128_div_ps(L1_PEND_MISS_PENDING,
	                                     MEM_LOAD_COMPLETED_L1_MISS_ANY);
	                  t1 = _mm128_div_ps(LD_BLOCKS_NO_SR,
	                                     CPU_CLK_UNHALTED_THREAD);
	                  t2 = _mm128_min_ps(_mm128_mul_ps(t0,t1),C1);
	                  metric = _mm128_mul_ps(C100,t2);
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
	           __m128 spr_fb_full_xmm4r4(const __m128 L1D_PEND_MISS_FB_FULL,
	                                     const __m128 CPU_CLK_UNHALTED_THREAD) {
	                  
	                   const __m128 C100 = _mm128_set1_ps(100.0f);
	                   register __m128 t0;
	                   register __m128 metric;
	                   t0 = _mm128_div_ps(L1D_PEND_MISS_FB_FULL,
	                                      CPU_CLK_UNHALTED_THREAD);
	                   metric = _mm128_mul_ps(C100,t0);
	                   return (metric);                              
	       }  
	       
	                                               
	       	__ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m128 spr_fb_full_xmm4r4(const float * __restrict pL1D_PEND_MISS_FB_FULL,
	                                  const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                  
	                   register __m128 L1D_PEND_MISS_FB_FULL = 
	                                           _mm128_loadu_ps(&pL1D_PEND_MISS_FB_FULL[0]);
	                   register __m128 CPU_CLK_UNHALTED_THREAD = 
	                                           _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                   const __m128 C100 = _mm128_set1_ps(100.0f);
	                   register __m128 t0;
	                   register __m128 metric;
	                   t0 = _mm128_div_ps(L1D_PEND_MISS_FB_FULL,
	                                      CPU_CLK_UNHALTED_THREAD);
	                   metric = _mm128_mul_ps(C100,t0);
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
	        __m128	 spr_l2_bound_xmm4r4(const __m128 MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                     const __m128 MEMORY_ACTIVITY_STALLS_L2_MISS ,
	                                     const __m128 CPU_CLK_UNHALTED_THREAD) {
	                 
	                 const __m128 C100 = _mm128_set1_ps(100.0f);
	                 register __m128 t0;
	                 register __m128 t1;
	                 register __m128 metric;
	                 t0 = _mm128_sub_ps(MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L2_MISS);
	                 t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm128_mul_ps(C100,t1);
	                 return (metric);                    
	       }  
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m128	 spr_l2_bound_xmm4r4(const float * __restrict pMEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                     const float * __restrict pMEMORY_ACTIVITY_STALLS_L2_MISS ,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                 register __m128 MEMORY_ACTIVITY_STALLS_L1D_MISS = 
	                                                _mm128_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L1D_MISS[0]);
	                 register __m128 MEMORY_ACTIVITY_STALLS_L2_MISS  =
	                                                _mm128_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L2_MISS[0]);
	                 register __m128 CPU_CLK_UNHALTED_THREAD         =
	                                                _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                 const __m128 C100 = _mm128_set1_ps(100.0f);
	                 register __m128 t0;
	                 register __m128 t1;
	                 register __m128 metric;
	                 t0 = _mm128_sub_ps(MEMORY_ACTIVITY_STALLS_L1D_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L2_MISS);
	                 t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm128_mul_ps(C100,t1);
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
	        __m128	 spr_l3_bound_xmm4r4(const __m128 MEMORY_ACTIVITY_STALLS_L2_MISS,
	                                     const __m128 MEMORY_ACTIVITY_STALLS_L3_MISS ,
	                                     const __m128 CPU_CLK_UNHALTED_THREAD) {
	                 
	                 const __m128 C100 = _mm128_set1_ps(100.0f);
	                 register __m128 t0;
	                 register __m128 t1;
	                 register __m128 metric;
	                 t0 = _mm128_sub_ps(MEMORY_ACTIVITY_STALLS_L2_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L3_MISS);
	                 t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm128_mul_ps(C100,t1);
	                 return (metric);                    
	       }  
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m128	 spr_l3_bound_xmm4r4(const float * __restrict pMEMORY_ACTIVITY_STALLS_L2_MISS,
	                                     const float * __restrict pMEMORY_ACTIVITY_STALLS_L3_MISS ,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                 register __m128 MEMORY_ACTIVITY_STALLS_L2_MISS = 
	                                                _mm128_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L2_MISS[0]);
	                 register __m128 MEMORY_ACTIVITY_STALLS_L3_MISS  =
	                                                _mm128_loadu_ps(&pMEMORY_ACTIVITY_STALLS_L3_MISS[0]);
	                 register __m128 CPU_CLK_UNHALTED_THREAD         =
	                                                _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                 const __m128 C100 = _mm128_set1_ps(100.0f);
	                 register __m128 t0;
	                 register __m128 t1;
	                 register __m128 metric;
	                 t0 = _mm128_sub_ps(MEMORY_ACTIVITY_STALLS_L2_MISS,
	                                    MEMORY_ACTIVITY_STALLS_L3_MISS);
	                 t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm128_mul_ps(C100,t1);
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
	        __m128 spr_sq_full_xmm4r4(const __m128 XQ_FULL_CYCLES,
	                                  const __m128 L1D_PEND_MISS_L2_STALLS,
	                                  const __m128 CPU_CLK_UNHALTED_THREAD) {
	               
	                 const __m128 C100 = _mm128_set1_ps(100.0f);
	                 register __m128 t0;
	                 register __m128 t1;
	                 register __m128 metric;
	                 t0 = _mm128_add_ps(XQ_FULL_CYCLES,
	                                    L1D_PEND_MISS_L2_STALLS);
	                 t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm128_mul_ps(C100,t1);
	                 return (metric);                       
	       } 
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m128 spr_sq_full_xmm4r4(const float * __restrict pXQ_FULL_CYCLES,
	                                  const float * __restrict pL1D_PEND_MISS_L2_STALLS,
	                                  const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	                 
	                 register __m128 XQ_FULL_CYCLES = 
	                                   _mm128_loadu_ps(&pXQ_FULL_CYCLES[0]);
	                 register __m128 L1D_PEND_MISS_L2_STALLS = 
	                                   _mm128_loadu_ps(&pL1D_PEND_MISS_L2_STALLS[0]);
	                 register __m128 CPU_CLK_UNHALTED_THREAD = 
	                                   _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	                 const __m128 C100 = _mm128_set1_ps(100.0f);
	                 register __m128 t0;
	                 register __m128 t1;
	                 register __m128 metric;
	                 t0 = _mm128_add_ps(XQ_FULL_CYCLES,
	                                    L1D_PEND_MISS_L2_STALLS);
	                 t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	                 metric = _mm128_mul_ps(C100,t1);
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
	        __m128 spr_mem_bw_xmm4r4(const __m128 CPU_CLK_UNHALTED_THREAD,
	                                 const __m128 OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4) {
	               
	               const __m128 C100 = _mm128_set1_ps(100.0f);
	               register __m128 t0;
	               register __m128 t1;
	               register __m128 metric;
	               t0 = _mm128_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                  OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4);
	               t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm128_mul_ps(C100,t1);
	               return (metric);                         
	       }  
	       
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m128 spr_mem_bw_xmm4r4(const float * __restrict pCPU_CLK_UNHALTED_THREAD,
	                                 const float * __restrict pOFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4) {
	               
	               register __m128 CPU_CLK_UNHALTED_THREAD = 
	                                      _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               register __m128 OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4 = 
	                                      _mm128_loadu_ps(&pOFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4[0]);
	               const __m128 C100 = _mm128_set1_ps(100.0f);
	               register __m128 t0;
	               register __m128 t1;
	               register __m128 metric;
	               t0 = _mm128_min_ps(CPU_CLK_UNHALTED_THREAD,
	                                  OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD_c4);
	               t1 = _mm128_div_ps(t0,CPU_CLK_UNHALTED_THREAD);
	               metric = _mm128_mul_ps(C100,t1);
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
	        __m128 spr_mba_stalls_xmm4r4(const __m128 INT_MISC_MBA_STALLS,
	                                     const __m128 CPU_CLK_UNHALTED_THREAD) {
	               
	               const __m128 C100 = _mm128_set1_ps(100.0f);
	               register __m128 t0;
	               register __m128 metric;
	               t0 = _mm128_div_ps(INT_MISC_MBA_STALLS,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm128_mul_ps(C100,t0);
	               return (metric);                              
	      }  
	      
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        static inline
	        __m128 spr_mba_stalls_xmm4r4(const float * __restrict pINT_MISC_MBA_STALLS,
	                                     const float * __restrict pCPU_CLK_UNHALTED_THREAD) {
	               
	               register __m128 INT_MISC_MBA_STALLS = 
	                                       _mm128_loadu_ps(&pINT_MISC_MBA_STALLS[0]);
	               register __m128 CPU_CLK_UNHALTED_THREAD = 
	                                       _mm128_loadu_ps(&pCPU_CLK_UNHALTED_THREAD[0]);
	               const __m128 C100 = _mm128_set1_ps(100.0f);
	               register __m128 t0;
	               register __m128 metric;
	               t0 = _mm128_div_ps(INT_MISC_MBA_STALLS,
	                                  CPU_CLK_UNHALTED_THREAD);
	               metric = _mm128_mul_ps(C100,t0);
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

               
	       
	                


} // gms








#endif /*__GMS_SPR_METRICS_XMM4R4*/
