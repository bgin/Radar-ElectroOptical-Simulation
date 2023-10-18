
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

#ifndef __GMS_SPR_METRICS_R4_HPP__
#define __GMS_SPR_METRICS_R4_HPP__


#include <immintrin.h>
#include "GMS_config.h"


namespace gms {


*
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
      "MetricName": "Retiring",
      "LegacyName": "metric_TMA_Retiring(%)",
      "Level": 1,
      "BriefDescription": "This category represents fraction of slots utilized by useful work i.e. issued uops that eventually get retired. Ideally; all pipeline slots would be attributed to the Retiring category.  Retiring of 100% would indicate the maximum Pipeline_Width throughput was achieved.  Maximizing Retiring typically increases the Instructions-per-cycle (see IPC metric). Note that a high Retiring value does not necessary mean there is no room for more performance.  For example; Heavy-operations or Microcode Assists are categorized under Retiring. They often indicate suboptimal performance and can often be optimized or avoided. ",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_retiring_r4(const float PERF_METRICS_RETIRING,
	                              const float PERF_METRICS_FRONTEND_BOUND,
	                              const float PERF_METRICS_BAD_SPECULATION,
	                              const float PERF_METRICS_BACKEND_BOUND) {
	                              
	               constexpr float C100 = 100.0f;
	               float bcad;
	               float t0;
	               float metric;
	               bcad = PERF_METRICS_FRONTEND_BOUND+
	                      PERF_METRICS_BAD_SPECULATION+
	                      PERF_METRICS_RETIRING+
	                      PERF_METRICS_BACKEND_BOUND;
	               t0 =   PERF_METRICS_RETIRING/bcad;
	               metric = C100*t0;
	               return (metric);                     
	      }
	      
/*
    "MetricName": "Light_Operations",
      "LegacyName": "metric_TMA_..Light_Operations(%)",
      "ParentCategory": "Retiring",
      "Level": 2,
      "BriefDescription": "This metric represents fraction of slots where the CPU was retiring light-weight operations -- instructions that require no more than one uop (micro-operation). This correlates with total number of instructions used by the program. A uops-per-instruction (see UPI metric) ratio of 1 or less should be expected for decently optimized software running on Intel Core/Xeon products. While this often indicates efficient X86 instructions were executed; high value does not necessarily mean better performance cannot be achieved.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_light_ops_r4(const float PERF_METRICS_RETIRING,
	                               const float PERF_METRICS_FRONTEND_BOUND,
	                               const float PERF_METRICS_BAD_SPECULATION,
	                               const float PERF_METRICS_BACKEND_BOUND,
	                               const float PERF_METRICS_HEAVY_OPERATIONS) {
	                               
	              constexpr float C100 = 100.0f;
	              float bcad;
	              float t0;
	              float t1;
	              float metric;
	              bcad = PERF_METRICS_FRONTEND_BOUND+
	                     PERF_METRICS_BAD_SPECULATION+
	                     PERF_METRICS_RETIRING+
	                     PERF_METRICS_BACKEND_BOUND;
	              t0 =  PERF_METRICS_RETIRING/bcad;
	              t1 =  PERF_METRICS_HEAVY_OPERATIONS/bcad;
	              metric = C100*std::max(0.0f,t0-t1);
	              return (metric);                       
	     }
	     
/*
    "MetricName": "X87_Use",
      "LegacyName": "metric_TMA_......X87_Use(%)",
      "ParentCategory": "FP_Arith",
      "Level": 4,
      "BriefDescription": "This metric serves as an approximation of legacy x87 usage. It accounts for instructions beyond X87 FP arithmetic operations; hence may be used as a thermometer to avoid X87 high usage and preferably upgrade to modern ISA. Tip: consider compiler flags to generate newer AVX (or SSE) instruction sets; which typically perform better and feature vectors.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_x87_use_r4(const float PERF_METRICS_RETIRING,
	                             const float PERF_METRICS_FRONTEND_BOUND,
	                             const float PERF_METRICS_BAD_SPECULATION,
	                             const float PERF_METRICS_BACKEND_BOUND,
	                             const float UOPS_EXECUTED_X87,
	                             const float UOPS_EXECUTED_THREAD) {
	              
	              constexpr float C100 = 100.0f;
	              float bcad;
	              float t0;
	              float t1;
	              float metric;  
	              bcad = PERF_METRICS_FRONTEND_BOUND+
	                     PERF_METRICS_BAD_SPECULATION+
	                     PERF_METRICS_RETIRING+
	                     PERF_METRICS_BACKEND_BOUND;
	              t0   = PERF_METRICS_RETIRING/bcad;
	              t1   = t0*UOPS_EXECUTED_X87/
	                     UOPS_EXECUTED_THREAD;
	              metric = C100*t1;
	              return (metric);
	       }
	       
	      
/*
    "  MetricName": "FP_Scalar",
      "LegacyName": "metric_TMA_......FP_Scalar(%)",
      "ParentCategory": "FP_Arith",
      "Level": 4,
      "BriefDescription": "This metric approximates arithmetic floating-point (FP) scalar uops fraction the CPU has executed (retired). May overcount due to FMA double counting.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_fp_scalar_r4(const float FP_ARITH_INST_RETIRED_SCALAR_SINGLE_u0x03,
	                               const float FP_ARITH_INST_RETIRED2_SCALAR,
	                               const float PERF_METRICS_RETIRING,
	                               const float PERF_METRICS_FRONTEND_BOUND,
	                               const float PERF_METRICS_BAD_SPECULATION,
	                               const float PERF_METRICS_BACKEND_BOUND,
	                               const float TOPDOWN_SLOTS_perf_metrics) {
	              
	              constexpr float C100 = 100.0f;
	              float ab;
	              float decf;
	              float t0;
	              float t1;
	              float metric;
	              decf = PERF_METRICS_FRONTEND_BOUND+
	                     PERF_METRICS_BAD_SPECULATION+
	                     PERF_METRICS_RETIRING+
	                     PERF_METRICS_BACKEND_BOUND;
	              ab   = FP_ARITH_INST_RETIRED_SCALAR_SINGLE_u0x03+
	                     FP_ARITH_INST_RETIRED2_SCALAR;
	              t0   = PERF_METRICS_RETIRING/decf;
	              t1   = ab/t0;
	              metric = C100*t1*TOPDOWN_SLOTS_perf_metrics;
	              return (metric);                     
	       }
	       
/*
      "MetricName": "FP_Vector",
      "LegacyName": "metric_TMA_......FP_Vector(%)",
      "ParentCategory": "FP_Arith",
      "Level": 4,
      "BriefDescription": "This metric approximates arithmetic floating-point (FP) vector uops fraction the CPU has executed (retired) aggregated across all vector widths. May overcount due to FMA double counting.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_fp_vector_r4(const float FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE_u0x3c,
	                               const float FP_ARITH_INST_RETIRED2_VECTOR,
	                               const float PERF_METRICS_RETIRING,
	                               const float PERF_METRICS_FRONTEND_BOUND,
	                               const float PERF_METRICS_BAD_SPECULATION,
	                               const float PERF_METRICS_BACKEND_BOUND,
	                               const float TOPDOWN_SLOTS_perf_metrics) {
	                               
	              constexpr float C100 = 100.0f;
	              float ab;
	              float decf;
	              float t0;
	              float t1;
	              float metric;
	              decf = PERF_METRICS_FRONTEND_BOUND+
	                     PERF_METRICS_BAD_SPECULATION+
	                     PERF_METRICS_RETIRING+
	                     PERF_METRICS_BACKEND_BOUND;
	              ab   = FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE_u0x3c+
	                     FP_ARITH_INST_RETIRED2_VECTOR;
	              t0   = PERF_METRICS_RETIRING/decf; 
	              t1   = ab/t0;
	              metric = C100*std::min(t1*TOPDOWN_SLOTS_perf_metrics,1.0f);
	              return (metric);                        
	     }
	     
/*
      "MetricName": "FP_Vector_128b",
      "LegacyName": "metric_TMA_........FP_Vector_128b(%)",
      "ParentCategory": "FP_Vector",
      "Level": 5,
      "BriefDescription": "This metric approximates arithmetic FP vector uops fraction the CPU has retired for 128-bit wide vectors. May overcount due to FMA double counting.",
      "UnitOfMeasure": "percent",   
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_fp_vector_128b_r4(const float FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE,
	                                    const float FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE,
	                                    const float FP_ARITH_INST_RETIRED2_128B_PACKED_HALF,
	                                    const float PERF_METRICS_RETIRING,
	                                    const float PERF_METRICS_FRONTEND_BOUND,
	                                    const float PERF_METRICS_BAD_SPECULATION,
	                                    const float PERF_METRICS_BACKEND_BOUND,
	                                    const float TOPDOWN_SLOTS_perf_metrics) {
	                                    
	                 constexpr float C100 = 100.0f;
	                 float abc;
	                 float efdg;
	                 float t0;
	                 float t1;
	                 float metric;
	                 abc = FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE+
	                       FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE+
	                       FP_ARITH_INST_RETIRED2_128B_PACKED_HALF;
	                 efdg = PERF_METRICS_FRONTEND_BOUND+
	                        PERF_METRICS_BAD_SPECULATION+
	                        PERF_METRICS_RETIRING+
	                        PERF_METRICS_BACKEND_BOUND;
	                 t0   = PERF_METRICS_RETIRING/efdg;
	                 t1   = abc/t0*TOPDOWN_SLOTS_perf_metrics;
	                 metric = C100*std::min(t1,1.0f);
	                 return (metric);                           
	    }
	    
/*
   "MetricName": "FP_Vector_256b",
      "LegacyName": "metric_TMA_........FP_Vector_256b(%)",
      "ParentCategory": "FP_Vector",
      "Level": 5,
      "BriefDescription": "This metric approximates arithmetic FP vector uops fraction the CPU has retired for 256-bit wide vectors. May overcount due to FMA double counting.",
      "UnitOfMeasure": "percent",
*/
	        
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_fp_vector_256b_r4(const float FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE,
	                                    const float FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE,
	                                    const float FP_ARITH_INST_RETIRED2_256B_PACKED_HALF,
	                                    const float PERF_METRICS_RETIRING,
	                                    const float PERF_METRICS_FRONTEND_BOUND,
	                                    const float PERF_METRICS_BAD_SPECULATION,
	                                    const float PERF_METRICS_BACKEND_BOUND,
	                                    const float TOPDOWN_SLOTS_perf_metrics) {
	                                    
	                 constexpr float C100 = 100.0f;
	                 float abc;
	                 float efdg;
	                 float t0;
	                 float t1;
	                 float metric;
	                 abc = FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE+
	                       FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE+
	                       FP_ARITH_INST_RETIRED2_256B_PACKED_HALF;
	                 efdg = PERF_METRICS_FRONTEND_BOUND+
	                        PERF_METRICS_BAD_SPECULATION+
	                        PERF_METRICS_RETIRING+
	                        PERF_METRICS_BACKEND_BOUND;
	                 t0   = PERF_METRICS_RETIRING/efdg;
	                 t1   = abc/t0*TOPDOWN_SLOTS_perf_metrics;
	                 metric = C100*std::min(t1,1.0f);
	                 return (metric);                           
	    }  
	    
/*
    "MetricName": "FP_Vector_512b",
      "LegacyName": "metric_TMA_........FP_Vector_512b(%)",
      "ParentCategory": "FP_Vector",
      "Level": 5,
      "BriefDescription": "This metric approximates arithmetic FP vector uops fraction the CPU has retired for 512-bit wide vectors. May overcount due to FMA double counting.",
      "UnitOfMeasure": "percent",
*/   

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_fp_vector_512b_r4(const float FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE,
	                                    const float FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE,
	                                    const float FP_ARITH_INST_RETIRED2_512B_PACKED_HALF,
	                                    const float PERF_METRICS_RETIRING,
	                                    const float PERF_METRICS_FRONTEND_BOUND,
	                                    const float PERF_METRICS_BAD_SPECULATION,
	                                    const float PERF_METRICS_BACKEND_BOUND,
	                                    const float TOPDOWN_SLOTS_perf_metrics) {
	                                    
	                 constexpr float C100 = 100.0f;
	                 float abc;
	                 float efdg;
	                 float t0;
	                 float t1;
	                 float metric;
	                 abc = FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE+
	                       FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE+
	                       FP_ARITH_INST_RETIRED2_512B_PACKED_HALF;
	                 efdg = PERF_METRICS_FRONTEND_BOUND+
	                        PERF_METRICS_BAD_SPECULATION+
	                        PERF_METRICS_RETIRING+
	                        PERF_METRICS_BACKEND_BOUND;
	                 t0   = PERF_METRICS_RETIRING/efdg;
	                 t1   = abc/t0*TOPDOWN_SLOTS_perf_metrics;
	                 metric = C100*std::min(t1,1.0f);
	                 return (metric);                           
	    }  
         
/*
         "MetricName": "FP_AMX",
      "LegacyName": "metric_TMA_......FP_AMX(%)",
      "ParentCategory": "FP_Arith",
      "Level": 4,
      "BriefDescription": "This metric approximates arithmetic floating-point (FP) matrix uops fraction the CPU has retired (aggregated across all supported FP datatypes in AMX engine). Refer to AMX_Busy and GFLOPs metrics for actual AMX utilization and FP performance, resp.",
      "UnitOfMeasure": "percent",
*/  

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_fp_amx_r4(  const float AMX_OPS_RETIRED_BF16_c1,
	                              const float PERF_METRICS_RETIRING,
	                              const float PERF_METRICS_FRONTEND_BOUND,
	                              const float PERF_METRICS_BAD_SPECULATION,
	                              const float PERF_METRICS_BACKEND_BOUND,
	                              const float TOPDOWN_SLOTS_perf_metrics) {
	                              
	              constexpr float C100 = 100.0f;
	              float cdbe;
	              float t0;   
	              float t1;
	              float metric;
	              cdbe = PERF_METRICS_FRONTEND_BOUND+
	                     PERF_METRICS_BAD_SPECULATION+
	                     PERF_METRICS_RETIRING+
	                     PERF_METRICS_BACKEND_BOUND;
	              t0   = PERF_METRICS_RETIRING/cdbe;
	              t1   = AMX_OPS_RETIRED_BF16_c1/(t0*
	                     TOPDOWN_SLOTS_perf_metrics);
	              metric = C100*t1;
	              return (metric);                      
	     }   
	     
/*
      "MetricName": "Int_Vector_128b",
      "LegacyName": "metric_TMA_......Int_Vector_128b(%)",
      "ParentCategory": "Int_Operations",
      "Level": 4,
      "BriefDescription": "This metric represents 128-bit vector Integer ADD/SUB/SAD or VNNI (Vector Neural Network Instructions) uops fraction the CPU has retired.",
      "UnitOfMeasure": "percent",
*/   

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_int_vector_128b_r4( const float INT_VEC_RETIRED_ADD_128,
	                                      const float INT_VEC_RETIRED_VNNI_128,
	                                      const float PERF_METRICS_RETIRING,
	                                      const float PERF_METRICS_FRONTEND_BOUND,
	                                      const float PERF_METRICS_BAD_SPECULATION,
	                                      const float PERF_METRICS_BACKEND_BOUND,
	                                      const float TOPDOWN_SLOTS_perf_metrics) {
	               
	                constexpr float C100 = 100.0f;
	                float decf;
	                float ab;
	                float t0;
	                float t1;
	                float metric;
	                ab   =  INT_VEC_RETIRED_ADD_128+
	                        INT_VEC_RETIRED_VNNI_128;
	                decf =  PERF_METRICS_FRONTEND_BOUND+
	                        PERF_METRICS_BAD_SPECULATION+
	                        PERF_METRICS_RETIRING+
	                        PERF_METRICS_BACKEND_BOUND;
	                t0   =  PERF_METRICS_RETIRING/decf;
	                t1   =  ab/(t0*TOPDOWN_SLOTS_perf_metrics);  
	                metric = C100*t1;
	                return (metric);                       
	       } 
	       
/*
       "MetricName": "Int_Vector_256b",
      "LegacyName": "metric_TMA_......Int_Vector_256b(%)",
      "ParentCategory": "Int_Operations",
      "Level": 4,
      "BriefDescription": "This metric represents 256-bit vector Integer ADD/SUB/SAD or VNNI (Vector Neural Network Instructions) uops fraction the CPU has retired.",
      "UnitOfMeasure": "percent",
*/ 

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_int_vector_256b_r4( const float INT_VEC_RETIRED_ADD_256,
	                                      const float INT_VEC_RETIRED_MUL_256,
	                                      const float INT_VEC_RETIRED_VNNI_256,
	                                      const float PERF_METRICS_RETIRING,
	                                      const float PERF_METRICS_FRONTEND_BOUND,
	                                      const float PERF_METRICS_BAD_SPECULATION,
	                                      const float PERF_METRICS_BACKEND_BOUND,
	                                      const float TOPDOWN_SLOTS_perf_metrics) {
	                                      
	                  constexpr float C100 = 100.0f;
	                  float abc;
	                  float efdg;
	                  float t0;
	                  float t1;
	                  float metric;
	                  abc = INT_VEC_RETIRED_ADD_256+
	                        INT_VEC_RETIRED_MUL_256+
	                        INT_VEC_RETIRED_VNNI_256;
	                  efdg = PERF_METRICS_FRONTEND_BOUND+
	                         PERF_METRICS_BAD_SPECULATION+
	                         PERF_METRICS_RETIRING+
	                         PERF_METRICS_BACKEND_BOUND;
	                  t0  =  PERF_METRICS_RETIRING/efdg*
	                         TOPDOWN_SLOTS_perf_metrics;
	                  t1  =  abc/t0;
	                  metric = C100*t1;
	                  return (metric);                           
	    }
	    
/*
     "MetricName": "Int_AMX",
      "LegacyName": "metric_TMA_......Int_AMX(%)",
      "ParentCategory": "Int_Operations",
      "Level": 4,
      "BriefDescription": "This metric approximates arithmetic Integer (Int) matrix uops fraction the CPU has retired (aggregated across all supported Int datatypes in AMX engine). Refer to AMX_Busy and TIOPs metrics for actual AMX utilization and Int performance, resp.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_int_amx_r4( const float AMX_OPS_RETIRED_INT8_c1,
	                              const float PERF_METRICS_RETIRING,
	                              const float PERF_METRICS_FRONTEND_BOUND,
	                              const float PERF_METRICS_BAD_SPECULATION,
	                              const float PERF_METRICS_BACKEND_BOUND,
	                              const float TOPDOWN_SLOTS_perf_metrics) {
	                              
	              constexpr float C100 = 100.0f;
	              float cdbe;
	              float t0;
	              float t1;
	              float metric;
	              cdbe =  PERF_METRICS_FRONTEND_BOUND+
	                      PERF_METRICS_BAD_SPECULATION+
	                      PERF_METRICS_RETIRING+
	                      PERF_METRICS_BACKEND_BOUND;
	              t0   =  PERF_METRICS_RETIRING/cdbe*  
	                      TOPDOWN_SLOTS_perf_metrics;
	              t1   =  AMX_OPS_RETIRED_INT8_c1/t0;
	              metric = C100*t1;
	              return (metric);              
	      }
	      
/*
      "MetricName": "Shuffles",
      "LegacyName": "metric_TMA_......Shuffles(%)",
      "ParentCategory": "Int_Operations",
      "Level": 4,
      "BriefDescription": "This metric represents Shuffle (cross \"vector lane\" data transfers) uops fraction the CPU has retired.",
      "UnitOfMeasure": "percent",
*/

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_shuffles_r4(const float INT_VEC_RETIRED_SHUFFLES, 
	                              const float PERF_METRICS_RETIRING,
	                              const float PERF_METRICS_FRONTEND_BOUND,
	                              const float PERF_METRICS_BAD_SPECULATION,
	                              const float PERF_METRICS_BACKEND_BOUND,
	                              const float TOPDOWN_SLOTS_perf_metrics) {
	          
	              constexpr float C100 = 100.0f;
	              float cdbe;
	              float t0;
	              float t1;
	              float metric;
	              cdbe =  PERF_METRICS_FRONTEND_BOUND+
	                      PERF_METRICS_BAD_SPECULATION+
	                      PERF_METRICS_RETIRING+
	                      PERF_METRICS_BACKEND_BOUND;
	              t0   =  PERF_METRICS_RETIRING/cdbe*  
	                      TOPDOWN_SLOTS_perf_metrics;
	              t1   =  INT_VEC_RETIRED_SHUFFLES/t0;
	              metric = C100*t1;
	              return (metric);                             
	      }
	      
/*
      "MetricName": "Memory_Operations",
      "LegacyName": "metric_TMA_....Memory_Operations(%)",
      "ParentCategory": "Light_Operations",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of slots where the CPU was retiring memory operations -- uops for memory load or store accesses.",
      "UnitOfMeasure": "percent",
*/

        
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_mem_ops_r4( const float PERF_METRICS_RETIRING,
	                              const float PERF_METRICS_FRONTEND_BOUND,
	                              const float PERF_METRICS_BAD_SPECULATION,
	                              const float PERF_METRICS_BACKEND_BOUND,
	                              const float PERF_METRICS_HEAVY_OPERATIONS,
	                              const float MEM_UOP_RETIRED_ANY,
	                              const float TOPDOWN_SLOTS_perf_metrics) {
	           
	                constexpr float C100 = 100.0f;
	                float bcad;
	                float t0;
	                float t1;
	                float t2;
	                float t3;
	                float metric;
	                bcad = PERF_METRICS_FRONTEND_BOUND+
	                       PERF_METRICS_BAD_SPECULATION+
	                       PERF_METRICS_RETIRING+
	                       PERF_METRICS_BACKEND_BOUND;
	                t0   = PERF_METRICS_RETIRING/bcad;
	                t1   = PERF_METRICS_HEAVY_OPERATIONS/
	                       bcad;
	                t2   = std::max(0.0f,t0-t1);
	                t3   = MEM_UOP_RETIRED_ANY/t0;
	                metric = C100*t2*t3*TOPDOWN_SLOTS_perf_metrics;
	                return (metric);                    
	     }    
	     
/*
     "MetricName": "Fused_Instructions",
      "LegacyName": "metric_TMA_....Fused_Instructions(%)",
      "ParentCategory": "Light_Operations",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of slots where the CPU was retiring fused instructions -- where one uop can represent multiple contiguous instructions. The instruction pairs of CMP+JCC or DEC+JCC are commonly used examples.",
      "UnitOfMeasure": "percent",
*/	

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_fused_instr_r4( const float PERF_METRICS_RETIRING,
	                                  const float PERF_METRICS_FRONTEND_BOUND,
	                                  const float PERF_METRICS_BAD_SPECULATION,
	                                  const float PERF_METRICS_BACKEND_BOUND,
	                                  const float PERF_METRICS_HEAVY_OPERATIONS,
	                                  const float INST_RETIRED_MACRO_FUSED,
	                                  const float TOPDOWN_SLOTS_perf_metrics) {
	           
	                constexpr float C100 = 100.0f;
	                float bcad;
	                float t0;
	                float t1;
	                float t2;
	                float t3;
	                float metric;
	                bcad = PERF_METRICS_FRONTEND_BOUND+
	                       PERF_METRICS_BAD_SPECULATION+
	                       PERF_METRICS_RETIRING+
	                       PERF_METRICS_BACKEND_BOUND;
	                t0   = PERF_METRICS_RETIRING/bcad;
	                t1   = PERF_METRICS_HEAVY_OPERATIONS/
	                       bcad;
	                t2   = std::max(0.0f,t0-t1);
	                t3   = INST_RETIRED_MACRO_FUSED/t0;
	                metric = C100*t2*t3*TOPDOWN_SLOTS_perf_metrics;
	                return (metric);                    
	     }   
	     
/*
      "MetricName": "Non_Fused_Branches",
      "LegacyName": "metric_TMA_....Non_Fused_Branches(%)",
      "ParentCategory": "Light_Operations",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of slots where the CPU was retiring branch instructions that were not fused. Non-conditional branches like direct JMP or CALL would count here. Can be used to examine fusible conditional jumps that were not fused.",
      "UnitOfMeasure": "percent",
*/ 

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_non_fused_br_r4(const float PERF_METRICS_RETIRING,
	                                  const float PERF_METRICS_FRONTEND_BOUND,
	                                  const float PERF_METRICS_BAD_SPECULATION,
	                                  const float PERF_METRICS_BACKEND_BOUND,
	                                  const float PERF_METRICS_HEAVY_OPERATIONS,
	                                  const float BR_INST_RETIRED_ALL_BRANCHES,
	                                  const float INST_RETIRED_MACRO_FUSED,
	                                  const float TOPDOWN_SLOTS_perf_metrics) {
	                                  
	                constexpr float C100 = 100.0f;
	                float bcad;
	                float t0;
	                float t1;
	                float t2;
	                float t3;
	                float metric;
	                bcad = PERF_METRICS_FRONTEND_BOUND+
	                       PERF_METRICS_BAD_SPECULATION+
	                       PERF_METRICS_RETIRING+
	                       PERF_METRICS_BACKEND_BOUND;  
	                t0   = PERF_METRICS_RETIRING/bcad;
	                t1   = PERF_METRICS_HEAVY_OPERATIONS/
	                       bcad;
	                t2   = BR_INST_RETIRED_ALL_BRANCHES-
	                       INST_RETIRED_MACRO_FUSED;
	                t3   = std::max(0.0f,t0-t1)*(t2/t0);
	                metric = C100*t3*TOPDOWN_SLOTS_perf_metrics;
	                return (metric);                       
	      }   
	      
/*
      "MetricName": "Nop_Instructions",
      "LegacyName": "metric_TMA_....Nop_Instructions(%)",
      "ParentCategory": "Light_Operations",
      "Level": 3,
      "BriefDescription": "This metric represents fraction of slots where the CPU was retiring NOP (no op) instructions. Compilers often use NOPs for certain address alignments - e.g. start address of a function or loop body.",
      "UnitOfMeasure": "percent", 
*/      

                __ATTR_HOT__
	        __ATTR_ALIGN__(32)
             	static inline
	        float spr_nop_instr_r4(   const float PERF_METRICS_RETIRING,
	                                  const float PERF_METRICS_FRONTEND_BOUND,
	                                  const float PERF_METRICS_BAD_SPECULATION,
	                                  const float PERF_METRICS_BACKEND_BOUND,
	                                  const float PERF_METRICS_HEAVY_OPERATIONS,
	                                  const float INST_RETIRED_NOP,
	                                  const float TOPDOWN_SLOTS_perf_metrics) {
	           
	                constexpr float C100 = 100.0f;
	                float bcad;
	                float t0;
	                float t1;
	                float t2;
	                float t3;
	                float metric;
	                bcad = PERF_METRICS_FRONTEND_BOUND+
	                       PERF_METRICS_BAD_SPECULATION+
	                       PERF_METRICS_RETIRING+
	                       PERF_METRICS_BACKEND_BOUND;
	                t0   = PERF_METRICS_RETIRING/bcad;
	                t1   = PERF_METRICS_HEAVY_OPERATIONS/
	                       bcad;
	                t2   = std::max(0.0f,t0-t1);
	                t3   = INST_RETIRED_NOP/t0;
	                metric = C100*t2*t3*TOPDOWN_SLOTS_perf_metrics;
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
	        float  spr_heavy_ops_r4(const float PERF_METRICS_HEAVY_OPERATIONS,
	                                const float PERF_METRICS_FRONTEND_BOUND,
	                                const float PERF_METRICS_BAD_SPECULATION,
	                                const float PERF_METRICS_RETIRING,
	                                const float PERF_METRICS_BACKEND_BOUND) {
	                                
	               constexpr float C100 = 100.0f;
	               float bcde;
	               float t0;
	               float metric;
	               bcde = PERF_METRICS_FRONTEND_BOUND+
	                      PERF_METRICS_BAD_SPECULATION+
	                      PERF_METRICS_RETIRING+
	                      PERF_METRICS_BACKEND_BOUND;
	               t0   = PERF_METRICS_HEAVY_OPERATIONS/t0;
	               metric = C100*t0;
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
	        float  spr_few_uops_instr_r4( const float PERF_METRICS_HEAVY_OPERATIONS,
	                                      const float PERF_METRICS_FRONTEND_BOUND,
	                                      const float PERF_METRICS_BAD_SPECULATION,
	                                      const float PERF_METRICS_RETIRING,
	                                      const float PERF_METRICS_BACKEND_BOUND,
	                                      const float UOPS_RETIRED_MS,
	                                      const float TOPDOWN_SLOTS_perf_metrics) {
	               
	               constexpr float C100 = 100.0f;
	               float bcde;
	               float t0;
	               float t1;
	               float metric;
	               bcde = PERF_METRICS_FRONTEND_BOUND+
	                      PERF_METRICS_BAD_SPECULATION+
	                      PERF_METRICS_RETIRING+
	                      PERF_METRICS_BACKEND_BOUND; 
	               t0   = PERF_METRICS_HEAVY_OPERATIONS/bcde;
	               t1   = UOPS_RETIRED_MS/TOPDOWN_SLOTS_perf_metrics;
	               metric = C100*std::max(0.0f,t0-t1);
	               return (metric);                            
	      }                                 	         


}

#endif /*_GMS_SPR_METRICS_R4_HPP*/
