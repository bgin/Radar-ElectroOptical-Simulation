

#ifndef __GMS_HASWELL_MARCH_EVENTS_H__
#define __GMS_HASWELL_MARCH_EVENTS_H__


namespace file_info {

#if defined _WIN64
   #include "../GMS_version.h"
#elif defined __linux
   #include "GMS_version.h"
#endif
     const unsigned int gGMS_HASWELL_MARCH_EVENTS_MAJOR = gms::common::gVersionInfo.m_VersionMajor;
     const unsigned int gGMS_HASWELL_MARCH_EVENTS_MINOR = gms::common::gVersionInfo.m_VersionMinor;
     const unsigned int gGMS_HASWELL_MARCH_EVENTS_MICRO = gms::common::gVersionInfo.m_VersionMicro;
     const unsigned int gGMS_HASWELL_MARCH_EVENTS_FULLVER =
       1000U*gGMS_HASWELL_MARCH_EVENTS_MAJOR+100U*gGMS_HASWELL_MARCH_EVENTS_MINOR+10U*gGMS_HASWELL_MARCH_EVENTS_MICRO;
     const char * const pgGMS_HASWELL_MARCH_EVENTS_CREATE_DATE = "26-10-2019 18:31 +00200 (SAT 26 OCT 2019 GMT+2)";
     const char * const pgGMS_HASWELL_MARCH_EVENTS_BUILD_DATE  = "00-00-0000 00:00";
     const char * const pgGMS_HASWELL_MARCH_EVENTS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
     const char * const pgGMS_HASWELL_MARCH_EVENTS_SYNOPSIS    = "Intel Haswell family micro-architectural events (raw codes)";
     
}


/*
     Haswell PMC raw encoding definitions
     Based on haswell_core_v28.json file
*/
#if !defined (INST_RETIRED_ANY)
     #define INST_RETIRED_ANY                          0x100
#endif
#if !defined (CPU_CLK_UNHALTED_THREAD)
     #define CPU_CLK_UNHALTED_THREAD                   0x200
#endif
#if !defined (CPU_CLK_UNHALTED_THREAD_ANY)
     #define CPU_CLK_UNHALTED_THREAD_ANY               0x200
#endif
#if !defined (CPU_CLK_UNHALTED_REF_TSC)
     #define CPU_CLK_UNHALTED_REF_TSC                  0x300
#endif
#if !defined (LD_BLOCKS_STORE_FORWARD)
    #define LD_BLOCKS_STORE_FORWARD                    0x203
#endif
#if !defined (LD_BLOCKS_NO_SR)
    #define LD_BLOCKS_NO_SR                            0x803
#endif
#if !defined (MISALIGN_MEM_REF_LOADS)
    #define MISALIGN_MEM_REF_LOADS                     0x105
#endif
#if !defined (MISALIGN_MEM_REF_STORES)
    #define MISALIGN_MEM_REF_STORES                    0x205
#endif
#if !defined (LD_BLOCKS_PARTIAL_ADDRESS_ALIAS)
    #define LD_BLOCKS_PARTIAL_ADDRESS_ALIAS            0x107
#endif
#if !defined (DTLB_LOAD_MISSES_MISS_CAUSES_A_WALK)
    #define DTLB_LOAD_MISSES_MISS_CAUSES_A_WALK        0x108
#endif
#if !defined (DTLB_LOAD_MISSES_WALK_COMPLETED_4K)
    #define DTLB_LOAD_MISSES_WALK_COMPLETED_4K         0x208
#endif
#if !defined (DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M)
    #define DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M      0x408
#endif
#if !defined (DTLB_LOAD_MISSES_WALK_COMPLETED_1G)
    #define DTLB_LOAD_MISSES_WALK_COMPLETED_1G         0x608
#endif
#if !defined (DTLB_LOAD_MISSES_WALK_COMPLETED)
    #define DTLB_LOAD_MISSES_WALK_COMPLETED            0x808
#endif
#if !defined (DTLB_LOAD_MISSES_WALK_DURATION)
    #define DTLB_LOAD_MISSES_WALK_DURATION             0x1008
#endif
#if !defined (DTLB_LOAD_MISSES_STLB_HIT_4K)
    #define DTLB_LOAD_MISSES_STLB_HIT_4K               0x2008
#endif
#if !defined (DTLB_LOAD_MISSES_STLB_HIT_2M)
    #define DTLB_LOAD_MISSES_STLB_HIT_2M               0x4008
#endif
#if !defined (DTLB_LOAD_MISSES_STLB_HIT)
    #define DTLB_LOAD_MISSES_STLB_HIT                  0x6008
#endif
#if !defined (DTLB_LOAD_MISSES_PDE_CACHE_MISS)
    #define DTLB_LOAD_MISSES_PDE_CACHE_MISS            0x8008
#endif
#if !defined (INT_MISC_RECOVERY_CYCLES)
    #define INT_MISC_RECOVERY_CYCLES                   0x30D
#endif
#if !defined (UOPS_ISSUED_ANY)
    #define UOPS_ISSUED_ANY                            0x100E
#endif
#if !defined (UOPS_ISSUED_STALL_CYCLES)
    #define UOPS_ISSUED_STALL_CYCLES                   0x100E
#endif
#if !defined (UOPS_ISSUED_CORE_STALL_CYCLES)
    #define UOPS_ISSUED_CORE_STALL_CYCLES              0x100E
#endif
#if !defined (UOPS_ISSUED_FLAGS_MERGE)
    #define UOPS_ISSUED_FLAGS_MERGE                    0x100E
#endif
#if !defined (UOPS_ISSUED_SLOW_LEA)
    #define UOPS_ISSUED_SLOW_LEA                       0x200E
#endif
#if !defined (UOPS_ISSUED_SINGLE_MUL)
    #define UOPS_ISSUED_SINGLE_MUL                     0x400E
#endif
#if !defined (ARITH_DIVIDER_UOPS)
    #define ARITH_DIVIDER_UOPS                         0x0214
#endif
#if !defined (L2_RQSTS_DEMAND_DATA_RD_MISS)
    #define L2_RQSTS_DEMAND_DATA_RD_MISS               0x2124
#endif
#if !defined (L2_RQSTS_RFO_MISS)
    #define L2_RQSTS_RFO_MISS                          0x2224
#endif
#if !defined (L2_RQSTS_CODE_RD_MISS)
    #define L2_RQSTS_CODE_RD_MISS)                     0x2424
#endif
#if !defined (L2_RQSTS_ALL_DEMAND_MISS)
    #define L2_RQSTS_ALL_DEMAND_MISS                   0x2724
#endif
#if !defined (L2_RQSTS_L2_PF_MISS)
    #define L2_RQSTS_L2_PF_MISS                        0x3024
#endif
#if !defined (L2_RQSTS_MISS)
    #define L2_RQSTS_MISS                              0x3F24
#endif
#if !defined (L2_RQSTS_DEMAND_DATA_RD_HIT)
    #define L2_RQSTS_DEMAND_DATA_RD_HIT                0xC124
#endif
#if !defined (L2_RQSTS_RFO_HIT)
    #define L2_RQSTS_RFO_HIT                           0xC224
#endif
#if !defined (L2_RQSTS_CODE_RD_HIT)
    #define L2_RQSTS_CODE_RD_HIT                       0xC424
#endif
#if !defined (L2_RQSTS_L2_PF_HIT)
    #define L2_RQSTS_L2_PF_HIT                         0xD024
#endif
#if !defined (L2_RQSTS_ALL_DEMAND_DATA_RD)
    #define L2_RQSTS_ALL_DEMAND_DATA_RD                0xE124
#endif
#if !defined (L2_RQSTS_ALL_RFO)
    #define L2_RQSTS_ALL_RFO                           0xE224
#endif
#if !defined (L2_RQSTS_ALL_CODE_RD)
    #define L2_RQSTS_ALL_CODE_RD                       0xE424
#endif
#if !defined (L2_RQSTS_ALL_DEMAND_REFERENCES)
    #define L2_RQSTS_ALL_DEMAND_REFERENCES             0xE724
#endif
#if !defined (L2_RQSTS_ALL_PF)
    #define L2_RQSTS_ALL_PF                            0xF824
#endif
#if !defined (L2_RQSTS_REFERENCES)
    #define L2_RQSTS_REFERENCES                        0xFF24
#endif
#if !defined (L2_DEMAND_RQSTS_WB_HIT)
    #define L2_DEMAND_RQSTS_WB_HIT                     0x5027
#endif
#if !defined (LONGEST_LAT_CACHE_MISS)
    #define LONGEST_LAT_CACHE_MISS                     0x412E
#endif
#if !defined (LONGEST_LAT_CACHE_REFERENCE)
    #define LONGEST_LAT_CACHE_REFERENCE                0x4F2E
#endif
#if !defined (CPU_CLK_UNHALTED_THREAD_P)
#define CPU_CLK_UNHALTED_THREAD_P                      0x003C
    #endif
#if !defined (CPU_CLK_UNHALTED_THREAD_P_ANY)
    #define CPU_CLK_UNHALTED_THREAD_P_ANY              0x003C
#endif
#if !defined (CPU_CLK_UNHALTED_THREAD_REF_XCLK)
    #define CPU_CLK_UNHALTED_THREAD_REF_XCLK           0x013C
#endif
#if !defined (CPU_CLK_UNHALTED_THREAD_REF_XCLK_ANY)
    #define CPU_CLK_UNHALTED_THREAD_REF_XCLK_ANY       0x013C
#endif
#if !defined (CPU_CLK_UNHALTED_REF_XCLK)
    #define CPU_CLK_UNHALTED_REF_XCLK                  0x013C
#endif
#if !defined (CPU_CLK_UNHALTED_REF_XCLOCK_ANY)
    #define CPU_CLK_UNHALTED_REF_XCLOCK_ANY            0x13C
#endif
#if !defined (CPU_CLK_THREAD_UNHALTED_ONE_THREAD_ACTIVE)
    #define CPU_CLK_THREAD_UNHALTED_ONE_THREAD_ACTIVE  0x023C
#endif
#if !defined (L1D_MISS_PENDING)
    #define L1D_MISS_PENDING                           0x0148
#endif
#if !defined (L1D_PEND_MISS_PENDING_CYCLES)
    #define L1D_PEND_MISS_PENDING_CYCLES               0x0148
#endif
#if !defined (L1D_PEND_MISS_PENDING_CYCLES_ANY)
    #define L1D_PEND_MISS_PENDING_CYCLES_ANY           0x0148
#endif
#if !defined (L1D_PEND_MISS_REQUEST_FB_FULL)       
    #define L1D_PEND_MISS_REQUEST_FB_FULL              0x0248
#endif
#if !defined (L1D_PEND_MISS_FB_FULL)
    #define L1D_PEND_MISS_FB_FULL                      0x0248
#endif
#if !defined (DTLB_STORE_MISSES_MISS_CAUSES_A_WALK)
    #define DTLB_STORE_MISSES_MISS_CAUSES_A_WALK       0x0149
#endif
#if !defined (DTLB_STORE_MISSES_WALK_COMPLETED_4K)
    #define DTLB_STORE_MISSES_WALK_COMPLETED_4K        0x0249
#endif
#if !defined (DTLB_STORE_MISSES_WALK_COMPLETED_2M_4M)
    #define DTLB_STORE_MISSES_WALK_COMPLETED_2M_4M     0x0449
#endif
#if !defined (DTLB_STORE_MISSES_WALK_COMPLETED_1G)
    #define DTLB_STORE_MISSES_WALK__COMPLETED_1G       0x0849
#endif
#if !defined (DTLB_STORE_MISSES_WALK_COMPLETED)
    #define DTLB_STORE_MISSES_WALK_COMPLETED           0x0E49
#endif
#if !defined (DTLB_STORE_MISS_WALK_DURATION)
    #define DTLB_STORE_MISS_WALK_DURATION              0x1049
#endif
#if !defined (DTLB_STORE_MISSES_STLB_HIT_4K)
    #define DTLB_STORE_MISSES_STLB_HIT_4K              0x2049
#endif
#if !defined (DTLB_STORE_MISSES_STLB_HIT_2M)
    #define DTLB_STORE_MISSES_STLB_HIT_2M              0x4049
#endif
#if !defined (DTLB_STORE_MISSES_STLB_HIT)
    #define DTLB_STORE_MISSES_STLB_HIT                 0x6049
#endif
#if !defined (DTLB_STORE_MISSES_PDE_CACHE_MISS)
    #define DTLB_STORE_MISSES_PDE_CACHE_MISS           0x8049
#endif
#if !defined (LOAD_HIT_PRE_SW_PF)
    #define LOAD_HIT_PRE_SW_PF                         0x014C
#endif
#if !defined (LOAD_HIT_PRE_HW_PF)
    #define LOAD_HIT_PRE_HW_PF                         0x024C
#endif
#if !defined (EPT_WALK_CYCLES)
    #define EPT_WALK_CYCLES                            0x104F
#endif
#if !defined (L1D_REPLACEMENT) 
    #define L1D_REPLACEMENT                            0x0151
#endif
#if !defined (TX_MEM_ABORT_CONFLICT)
    #define TX_MEM_ABORT_CONFLICT                      0x0154
#endif
#if !defined (TX_MEM_ABORT_CAPACITY_WRITE)
    #define TX_MEM_ABORT_CAPACITY_WRITE                0x0254
#endif
#if !defined (TX_MEM_ABORT_HLE_STORE_TO_ELIDED_LOCK)
    #define TX_MEM_ABORT_HLE_STORE_TO_ELIDED_LOCK      0x0454
#endif
#if !defined (TX_MEM_ABORT_HLE_ELISION_BUFFER_NOT_EMPTY)
    #define  TX_MEM_ABORT_HLE_ELISION_BUFFER_NOT_EMPTY 0x0845
#endif
#if !defined (TX_MEM_ABORT_HLE_ELISION_BUFFER_MISMATCH)
    #define TX_MEM_ABORT_HLE_ELISION_BUFFER_MISMATCH   0x1054
#endif
#if !defined (TX_MEM_ABORT_HLE_ELISION_BUFFER_UNSUPPORTED_ALIGNMENT)
    #define  TX_MEM_ABORT_HLE_ELISION_BUFFER_UNSUPPORTED_ALIGNMENT 0x2054
#endif
#if !defined (TX_MEM_HLE_ELISION_BUFFER_FULL)
    #define TX_MEM_HLE_ELISION_BUFFER_FULL             0x4054
#endif
#if !defined ( MOVE_ELIMINATION_INT_ELIMINATED)
    #define  MOVE_ELIMINATION_INT_ELIMINATED           0x0158
#endif
#if !defined (MOVE_ELIMINATION_SIMD_ELIMINATED)
    #define MOVE_ELIMINATION_SIMD_ELIMINATED           0x0258
#endif
#if !defined (MOVE_ELIMINATION_INT_NOT_ELIMINATED)
    #define MOVE_ELIMINATION_INT_NOT_ELIMINATED        0x0458
#endif
#if !defined (MOVE_ELIMINATION_SIMD_NOT_ELIMINATED)
    #define MOVE_ELIMINATION_SIMD_NOT_ELIMINATED       0x0858
#endif
#if !defined (CPL_CYCLES_RING0)
    #define CPL_CYCLES_RING0                           0x015C
#endif
#if !defined (CPL_CYCLES_RING0_TRANS)
    #define  CPL_CYCLES_RING0_TRANS                    0x015C
#endif
#if !defined (CPL_CYCLES_RING123)
    #define CPL_CYCLES_RING123                         0x025C
#endif
#if !defined (TX_EXEC_MISC1)
    #define TX_EXEC_MISC1                              0x015D
#endif
#if !defined (TX_EXEC_MISC2)
    #define TX_EXEC_MISC2                              0x025D
#endif
#if !defined (TX_EXEC_MISC3)
    #define  TX_EXEC_MISC3                             0x045D
#endif
#if !defined (TX_EXEC_MISC4)
    #define TX_EXEC_MISC4                              0x085D
#endif
#if !defined (TX_EXEC_MISC5)
    #define TX_EXEC_MISC5                              0x105D
#endif
#if !defined (RS_EVENTS_EMPTY_CYCLES)
    #define RS_EVENTS_EMPTY_CYCLES                     0x015E
#endif
#if !defined ( RS_EVENTS_EMPTY_END)
    #define  RS_EVENTS_EMPTY_END                       0x015E
#endif
#if !defined (OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD)
   #define OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD 0x0160
#endif
#if !defined (OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_DATA_RD)
    #define  OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_DATA_RD 0x0160
#endif
#if !defined (OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD_GE_6)
    #define OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD_GE_6 0x0160
#endif
#if !defined (OFFCORE_REQUESTS_OUTSTANDING_DEMAND_CODE_R)
    #define OFFCORE_REQUESTS_OUTSTANDING_DEMAND_CODE_R 0x0260
#endif
#if !defined (OFFCORE_REQUESTS_OUTSTANDING_DEMAND_RFO)
    #define OFFCORE_REQUESTS_OUTSTANDING_DEMAND_RFO    0x0460
#endif
#if !defined (OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD)
    #define  OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD  0x0860
#endif
#if !defined (LOCK_CYCLES_SPLIT_LOCK_UC_LOCK_DURATION)
    #define LOCK_CYCLES_SPLIT_LOCK_UC_LOCK_DURATION    0x0263
#endif
#if !defined (IDQ_EMPTY)
    #define IDQ_EMPTY                                  0x0279
#endif
#if !defined (IDQ_MITE_UOPS)
    #define IDQ_MITE_UOPS                              0x0479
#endif
#if !defined (IDQ_MITE_CYCLES)
    #define IDQ_MITE_CYCLES                            0x0479
#endif
#if !defined (IDQ_DSB_UOPS)
    #define IDQ_DSB_UOPS                               0x0879
#endif
#if !defined (IDQ_MS_DSB_UOPS)
    #define IDQ_MS_DSB_UOPS                            0x1079
#endif

namespace gms {

       namespace hardware {

#if defined __linux           
    __attribute__((align(64)))
            struct HaswellUarchEvents_t {
              //   ! Based on haswell_core_v28.json
              //   ! rNNNN -- as used by the perf -- UMASK+EVENT
	      // --------------------------------------------------------------
                   //   ! "EventName": "INST_RETIRED.ANY",
                   //   !  "BriefDescription": "Instructions retired from execution."
	          const char INST_RETIRED_ANY[] = "r100";
	          // !  "EventName": "CPU_CLK_UNHALTED.THREAD",
                  // !  "BriefDescription": "Core cycles when the thread is not in halt state.",	  
                  const char CPU_CLK_UNHALTED_THREAD[] = "r200";
		  //    "EventName": "CPU_CLK_UNHALTED.THREAD_ANY",
		  //    "BriefDescription": "Core cycles when at least one thread on the physical core is not in halt state.",
		  const char CPU_CLK_UNHALTED_THREAD_ANY[] = "r200";
		  //    "EventName": "CPU_CLK_UNHALTED.REF_TSC",
		  //    "BriefDescription": "Reference cycles when the core is not in halt state.",
		  const char CPU_CLK_UNHALTED_REF_TSC[] = "r300";
		  //     "EventName": "LD_BLOCKS.STORE_FORWARD",
                  //  !  "BriefDescription": "loads blocked by overlapping with store buffer that cannot be forwarded",
		  const char LD_BLOCKS_STORE_FORWARD[] = "r203";
		  //  !   "EventName": "LD_BLOCKS.NO_SR",
                  //  !   "BriefDescription": "The number of times that split load operations
                  //  !   are temporarily blocked because all resources for handling the split accesses are in use",
		  const char LD_BLOCKS_NO_SR[] = "r803";
		  //  !    "EventName": "MISALIGN_MEM_REF.LOADS",
                  //  !    "BriefDescription": "Speculative cache line split load uops dispatched to L1 cache",
		  const char MISALIGN_MEM_REF_LOADS[] = "r105";
		  // !   "EventName": "MISALIGN_MEM_REF.STORES",
                  // !   "BriefDescription": "Speculative cache line split STA uops dispatched to L1 cache",
		  const char MISALIGN_MEM_REF_STORES[] = "r205";
		  // !   "EventName": "LD_BLOCKS_PARTIAL.ADDRESS_ALIAS",
                  // !   "BriefDescription": "False dependencies in MOB due to partial compare on address.",
                  // !   "PublicDescription": "Aliasing occurs when a load is issued after a store and their memory addresses are offset by 4K.
                  // !   This event counts the number of loads that aliased with a precedin  !
                  // !   g store, resulting in an extended address check in the pipeline which can have a performance impact."
		  const char LD_BLOCKS_PARTIAL_ADDRESS_ALIAS[] = "r107";
		  //  !  "EventName": "DTLB_LOAD_MISSES.MISS_CAUSES_A_WALK",
                  //  !  "BriefDescription": "Load misses in all DTLB levels that cause page walks",
		  const char DTLB_BLOCKS_PARTIAL_ADDRESS_ALIAS[] = "r108";
		  //  !  "EventName": "DTLB_LOAD_MISSES.WALK_COMPLETED_4K",
                  //  !  "BriefDescription": "Demand load Miss in all translation lookaside buffer (TLB)
		  //     levels causes a page walk that completes (4K)."
		  const char DTLB_LOAD_MISSES_WALK_COMPLETED_4K[] = "r208";
		  //   !   "EventName": "DTLB_LOAD_MISSES.WALK_COMPLETED_2M_4M",
                  //   !   "BriefDescription": "Demand load Miss in all translation
                  //   !   lookaside buffer (TLB) levels causes a page walk that completes (2M/4M)."
		  const char DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M[] = "r408";
		  //   !    "EventName": "DTLB_LOAD_MISSES.WALK_COMPLETED_1G",
                  //   !    "BriefDescription": "Load miss in all TLB levels causes a page walk that completes. (1G)"
		  const char DTLB_LOAD_MISSES_WALK_COMPLETED_1G[] = "r808";
		  //   !    "EventName": "DTLB_LOAD_MISSES.WALK_COMPLETED",
                  //   !    "BriefDescription": "Demand load Miss in all translation
                  //   !    lookaside buffer (TLB) levels causes a page walk that completes of any page size."
		  const char DTLB_LOAD_MISSES_WALK_COMPLETED[] = "re08";
		  //   !    "EventName": "DTLB_LOAD_MISSES.WALK_DURATION",
                  //   !    "BriefDescription": "Cycles when PMH is busy with page walks",
                  //   !    "PublicDescription": "This event counts cycles when the
                  //   !    page miss handler (PMH) is servicing page walks caused by DTLB load misses.",
		  const char DTLB_LOAD_MISSES_WALK_DURATION[] = "r1008";
		  //   !     "EventName": "DTLB_LOAD_MISSES.STLB_HIT_4K",
                  //   !     "BriefDescription": "Load misses that miss the  DTLB and hit the STLB (4K)",
                  //   !     "PublicDescription": "This event counts load operations from a
                  //   !     4K page that miss the first DTLB level but hit the second and do not cause page walks.",
		  const char DTLB_LOAD_MISSES_STLB_HIT_4K[] = "r2008";
		  //   !     "EventName": "DTLB_LOAD_MISSES.STLB_HIT_2M",
                  //   !     "BriefDescription": "Load misses that miss the  DTLB and hit the STLB (2M)",
                  //   !     "PublicDescription": "This event counts load operations from
                  //   !     a 2M page that miss the first DTLB level but hit the second and do not cause page walks.",
		  const char DTLB_LOAD_MISSES_STLB_HIT_2M[] = "r4008";
		  //   !     "EventName": "DTLB_LOAD_MISSES.STLB_HIT",
                  //   !     "BriefDescription": "Load operations that miss the first DTLB level but hit the second and do not cause page walks"
		  const char DTLB_LOAD_MISSES_STLB_HIT[] = "r6008";
		  //   !      "EventName": "DTLB_LOAD_MISSES.PDE_CACHE_MISS",
                  //   !      "BriefDescription": "DTLB demand load misses with low part of linear-to-physical address translation missed"
		  const char DTLB_LOAD_MISSES_PDE_CACHE_MISS[] = "r8008";
		  //   !      "EventName": "INT_MISC.RECOVERY_CYCLES",
                  //   !  "BriefDescription": "Core cycles the allocator was stalled due to recovery from earlier clear event for this thread (e.g. misprediction or memory nuke)",
                  //   !  "PublicDescription": "This event counts the number of cycles spent waiting for a recovery after an event such as a processor nuke, JEClear, assist, hle/rtm abort etc
		  const char INT_MISC_RECOVERY_CYCLES[] = "r30D";
		  //   !     "EventName": "UOPS_ISSUED.ANY",
                  //   !    "BriefDescription": "Uops that Resource Allocation Table (RAT) issues to Reservation Station (RS)",
                  //   !    "PublicDescription": "This event counts the number of uops issued by the
                  //   !    Front-end of the pipeline to the Back-end. This event is counted at the
                  //   !    allocation stage and will count both retired and non  !   !-retired uops.",
		  const char UOPS_ISSUED_ANY[] = "r100E";
		  //   !   "EventName": "UOPS_ISSUED.STALL_CYCLES",
                  //   !   "BriefDescription": "Cycles when Resource Allocation Table (RAT) does not issue Uops to Reservation Station (RS) for the thread.
		  const char UOPS_ISSUED_STALL_CYCLES[] = "r100E";
		  //   !    "EventName": "UOPS_ISSUED.CORE_STALL_CYCLES",
                  //   !    "BriefDescription": "Cycles when Resource Allocation Table (RAT) does not issue Uops to Reservation Station (RS) for all threads."
		  const char UOPS_ISSUED_CORE_STALL_CYCLES[] = "r100E";
		  //   !   "EventName": "UOPS_ISSUED.FLAGS_MERGE",
                  //   !   "BriefDescription": "Number of flags-merge uops being allocated. Such uops considered perf sensitive; added by GSR u-arch."
		  const char UOPS_ISSUED_FLAGS_MERGE[] = "r100E";
		  //   !    "EventName": "UOPS_ISSUED.SLOW_LEA",
                  //   !    "BriefDescription": "Number of slow LEA uops being allocated. A uop is generally considered SlowLea if it has 3 sources
		  const char UOPS_ISSUED_SLOW_LEA[] = "r200E";
		  //   !    "EventName": "UOPS_ISSUED.SINGLE_MUL",
                  //   !    "BriefDescription": "Number of Multiply packed/scalar single precision uops allocated"
		  const char UOPS_ISSUED_SINGLE_MUL[] = "r400E";
		  //   !    "EventName": "ARITH.DIVIDER_UOPS",
                  //   !    "BriefDescription": "Any uop executed by the Divider. (This includes all divide uops, sqrt, ...)"
		  const char ARITH_DIVIDER_UOPS[] = "r0214";
		  //   !     "EventName": "L2_RQSTS.DEMAND_DATA_RD_MISS",
                  //   !    "BriefDescription": "Demand Data Read miss L2, no rejects",
		  const char L2_RQSTS_DEMAND_DATA_RD_MISS[] = "r2124";
		  //   !    "EventName": "L2_RQSTS.RFO_MISS",
                  //   !  "BriefDescription": "RFO requests that miss L2 cache",
                  //   !  "PublicDescription": "Counts the number of store RFO requests that miss the L2 cache."
		  const char L2_RQSTS_RFO_MISS[] = "r2224";
		  //   !  "EventName": "L2_RQSTS.CODE_RD_MISS",
                  //   !  "BriefDescription": "L2 cache misses when fetching instructions",
		  const char L2_RQSTS_CODE_RD_MISS[] = "r2424";
		  //   !   "EventName": "L2_RQSTS.ALL_DEMAND_MISS",
                  //   !   "BriefDescription": "Demand requests that miss L2 cache",
		  const char L2_RQSTS_ALL_DEMAND_MISS[] = "r2724";
		  //   !   "EventName": "L2_RQSTS.L2_PF_MISS",
                  //   !  "BriefDescription": "L2 prefetch requests that miss L2 cache",
                  //   !  "PublicDescription": "Counts all L2 HW prefetcher requests that missed L2.
		  const char L2_RQSTS_L2_PF_MISS[] = "r3024";
		  //   !   "EventName": "L2_RQSTS.MISS",
                  //   !   "BriefDescription": "All requests that miss L2 cache"
		  const char L2_RQSTS_MISS[] = "r3F24";
		  //   !   "EventName": "L2_RQSTS.DEMAND_DATA_RD_HIT",
                  //   !  "BriefDescription": "Demand Data Read requests that hit L2 cache",
                  //   !  "PublicDescription": "Counts the number of demand Data Read requests, initiated by load instructions, that hit L2 cache"
		  const char L2_RQSTS_DEMAND_DATA_RD_HIT[] = "rC124";
		  //   !   "EventName": "L2_RQSTS.RFO_HIT",
                  //   !  "BriefDescription": "RFO requests that hit L2 cache",
                  //   !  "PublicDescription": "Counts the number of store RFO requests that hit the L2 cache."
		  const char L2_RQSTS_RFO_HIT[] = "rC224";
		  //   !  "EventName": "L2_RQSTS.CODE_RD_HIT",
                  //   !  "BriefDescription": "L2 cache hits when fetching instructions, code reads.",
		  const char L2_RQSTS_CODE_RD_HIT[] = "rC424";
		  //   !   "EventName": "L2_RQSTS.L2_PF_HIT",
                  //   !   "BriefDescription": "L2 prefetch requests that hit L2 cache",
		  const char L2_RQSTS_L2_PF_HIT[] = "rD024";
		  //    !    "EventName": "L2_RQSTS.ALL_DEMAND_DATA_RD",
                  //    !  "BriefDescription": "Demand Data Read requests",
                  //    !  "PublicDescription": "Counts any demand and L1 HW prefetch data load requests to L2."
		  const char L2_RQSTS_ALL_DEMAND_DATA_RD[] = "rE124";
		  //    !   "EventName": "L2_RQSTS.ALL_RFO",
                  //    !  "BriefDescription": "RFO requests to L2 cache"
		  const char L2_RQSTS_ALL_RFO[] = "rE224";
		  //    !   "EventName": "L2_RQSTS.ALL_CODE_RD",
                  //    !   "BriefDescription": "L2 code requests
		  const char L2_RQSTS_ALL_CODE_RD[] = "rE424";
		  //    !   "EventName": "L2_RQSTS.ALL_DEMAND_REFERENCES",
                  //    !   "BriefDescription": "Demand requests to L2 cache",
		  const char L2_RQSTS_ALL_DEMAND_REFERENCES[] = "r724";
		  //    !   "EventName": "L2_RQSTS.ALL_PF",
                  //    !   "BriefDescription": "Requests from L2 hardware prefetchers",
		  const char L2_RQSTS_ALL_PF[] = "r824";
		  //    !     "EventName": "L2_RQSTS.REFERENCES",
                  //    !   "BriefDescription": "All L2 requests"
		  const char L2_RQSTS_REFERENCES[] = "rFF24";
		  //    !  "EventName": "L2_DEMAND_RQSTS.WB_HIT",
                  //    !   "BriefDescription": "Not rejected writebacks that hit L2 cache"
		  const char L2_DEMAND_RQSTS_WB_HIT[] = "r5027";
		  //    !    "EventName": "LONGEST_LAT_CACHE.MISS",
                  //    !  "BriefDescription": "Core-originated cacheable demand requests missed L3"
		  const char LONGEST_LAT_CACHE_MISS[] = "r412E";
		  //    !  "EventName": "LONGEST_LAT_CACHE.REFERENCE",
                  //    !  "BriefDescription": "Core-originated cacheable demand requests that refer to L3"
		  const char LONGEST_LAT_CACHE_REFERENCE[] = "r4F2E";
		  //    !  "EventName": "CPU_CLK_UNHALTED.THREAD_P",
                  //    !  "BriefDescription": "Thread cycles when thread is not in halt state",
                  //    !  "PublicDescription": "Counts the number of thread cycles while the thread is not in a halt state. 
                  //    !   The thread enters the halt state when it is running the HLT instruction.
                  //    !   The core frequency may change from time to time due to power or thermal throttling.",
		  const char CPU_CLK_UNHALTED_THREAD_P[] = "r003C";
		  //    !   "EventName": "CPU_CLK_THREAD_UNHALTED.REF_XCLK",
                  //    !  "BriefDescription": "Reference cycles when the thread is unhalted (counts at 100 MHz rate)",
                  //    !  "PublicDescription": "Increments at the frequency of XCLK (100 MHz) when not halted."
		  const char CPU_CLK_THREAD_UNHALTED_REF_XCLK[] = "r013C";
		  //     !    "EventName": "CPU_CLK_THREAD_UNHALTED.REF_XCLK_ANY",
                  //     !  "BriefDescription": "Reference cycles when the at least one thread on the physical core is unhalted (counts at 100 MHz rate)"
		  const char CPU_CLK_THREAD_UNHALTED_REF_XCLK_ANY[] = "r013c";
		  //     !    "EventName": "CPU_CLK_UNHALTED.REF_XCLK",
                  //     !   "BriefDescription": "Reference cycles when the thread is unhalted (counts at 100 MHz rate)"
		  const char CPU_CLK_UNHALTED_REF_XCLK[] = "r013C";
		  //     !   "EventName": "CPU_CLK_THREAD_UNHALTED.ONE_THREAD_ACTIVE",
                  //     !   "BriefDescription": "Count XClk pulses when this thread is unhalted and the other thread is halted."
                  const char  CPU_CLK_THREAD_ONE_THREAD_ACTIVE[] = "r023C";
                  //     !    "EventName": "CPU_CLK_UNHALTED.ONE_THREAD_ACTIVE",
                  //     !   "BriefDescription": "Count XClk pulses when this thread is unhalted and the other thread is halted.
                  const char  CPU_CLK_UNHALTED_ONE_THREAD_ACTIVE[] = "r023C";
                  //     !   "EventName": "L1D_PEND_MISS.PENDING",
                  //     !  "BriefDescription": "L1D miss oustandings duration in cycles",
                  //     !  "PublicDescription": "Increments the number of outstanding L1D misses every cycle. Set Cmask = 1 and Edge =1 to count occurrences."
                  const char  L1D_MISS_PENDING[]  = "r0148";
                  //     !   "EventName": "L1D_PEND_MISS.PENDING_CYCLES",
                  //     !   "BriefDescription": "Cycles with L1D load Misses outstanding.",
                  const char  L1_PEND_MISS_PENDING_CYCLES[]  = "r0148";
                  //     !    "EventName": "L1D_PEND_MISS.PENDING_CYCLES_ANY",
                  //     !    "BriefDescription": "Cycles with L1D load Misses outstanding from any thread on physical core.
                  const char  L1_PEND_MISS_PENDING_CYCLES_ANY[]   = "r0148";
                  //     !   "EventName": "L1D_PEND_MISS.REQUEST_FB_FULL",
                  //     !  "BriefDescription": "Number of times a request needed a FB entry but there was no entry available for it.
                  //     !   That is the FB unavailability was dominant reason for blocking the request. A request includes cacheable/uncacheable demands that is load, store or SW prefetch. HWP are e.
                  const char  L1D_PEND_MISS_REQUEST_FB_FULL[]  = "r0248";
                  //     !   "EventName": "L1D_PEND_MISS.FB_FULL",
                  //     !   "BriefDescription": "Cycles a demand request was blocked due to Fill Buffers inavailability."
                  const char  L1D_PEND_MISS_FB_FULL[]  = "r0248";
                  //     !    "EventName": "DTLB_STORE_MISSES.MISS_CAUSES_A_WALK",
                  //     !  "BriefDescription": "Store misses in all DTLB levels that cause page walks"
                  const char  DTLB_STORE_MISES_MISS_CAUSES_A_WALK[]    = "r0149";
                  //     !   "EventName": "DTLB_STORE_MISSES.WALK_COMPLETED_4K",
                  //     !  "BriefDescription": "Store miss in all TLB levels causes a page walk that completes. (4K)"
                  const char  DTLB_STORE_MISSES_WALK_COMPLETED_4K[]    = "r0249";
                  //     !   "EventName": "DTLB_STORE_MISSES.WALK_COMPLETED_2M_4M",
                  //     !  "BriefDescription": "Store misses in all DTLB levels that cause completed page walks (2M/4M)",
                  const char  DTLB_STORE_MISSES_WALK_COMPLETED_2M_4M[] = "r0449";
                  //     !  "EventName": "DTLB_STORE_MISSES.WALK_COMPLETED_1G",
                  //     !  "BriefDescription": "Store misses in all DTLB levels that cause completed page walks. (1G)
                  const char  DTLB_STORE_MISSES_WALK_COMPLETED_1G[]    = "r0849";
                  //     !   "EventName": "DTLB_STORE_MISSES.WALK_COMPLETED",
                  //     !   "BriefDescription": "Store misses in all DTLB levels that cause completed page walks",
                  const char  DTLB_STORE_MISSES_WALK_COMPLETED[]       = "r0e49";
                  //     !  "EventName": "DTLB_STORE_MISSES.WALK_DURATION",
                  //     !  "BriefDescription": "Cycles when PMH is busy with page walks",
                  //     !  "PublicDescription": "This event counts cycles when the  page miss handler (PMH) is servicing page walks caused by DTLB store misses.
                  const char  DTLB_STORE_MISSES_WALK_DURATION[]        = "r1049";
                  //     !   "EventName": "DTLB_STORE_MISSES.STLB_HIT_4K",
                  //     !  "BriefDescription": "Store misses that miss the  DTLB and hit the STLB (4K)",
                  const char  DTLB_STORE_MISSES_STLB_HIT_4K[]          = "r2049";
                  //     !  "EventName": "DTLB_STORE_MISSES.STLB_HIT_2M",
                  //     !  "BriefDescription": "Store misses that miss the  DTLB and hit the STLB (2M)",
                  //     !  "PublicDescription": "This event counts store operations from
		  //        a 2M page that miss the first DTLB level but hit the second and do not cause page walks.",
                  const char  DTLB_STORE_MISSES_STLB_HIT_2M[]          = "r4049";
		  //     !   "EventName": "DTLB_STORE_MISSES.STLB_HIT",
                  //     !  "BriefDescription": "Store operations that miss the first TLB level but hit the second and do not cause page walks",
                  const char  DTLB_STORE_MISSES_STLB_HIT[]             = "r6049";
                  //     !   "EventName": "DTLB_STORE_MISSES.PDE_CACHE_MISS",
                  //     !  "BriefDescription": "DTLB store misses with low part of linear-to-physical address translation missed"
                  const char  DTLB_STORE_MISSES_PDE_CACHE_MISS[]       = "r8049";
                  //     !   "EventName": "LOAD_HIT_PRE.SW_PF",
                  //     !  "BriefDescription": "Not software-prefetch load dispatches that hit FB allocated for software prefetch"
                  const char  LOAD_HIT_PRE_SW_PF[]                     = "r014C";
                  //     !   "EventName": "LOAD_HIT_PRE.HW_PF",
                  //     !  "BriefDescription": "Not software-prefetch load dispatches that hit FB allocated for hardware prefetch"
                  const char  LOAD_HIT_PRE_HW_PF[]                     = "r024C";
                  //     !   "EventName": "EPT.WALK_CYCLES",
                  //     !   "BriefDescription": "Cycle count for an Extended Page table walk.
                  const char  EPT_WALK_CYCLES[]                        = "r104F";
                  //     !  "EventName": "L1D.REPLACEMENT",
                  //     !  "BriefDescription": "L1D data line replacements",
                  const char  L1D_REPLACEMENT[]                        = "r0151";
                  //     !   "EventName": "TX_MEM.ABORT_CONFLICT",
                  //     !  "BriefDescription": "Number of times a transactional abort was signaled due to a data conflict on a transactionally accessed address.
                  const char  TX_MEM_ABORT_CONFLICT[]                  = "r0154";
                  //     !  "EventName": "TX_MEM.ABORT_CAPACITY_WRITE",
                  //     !  "BriefDescription": "Number of times a transactional abort was signaled due to a data capacity limitation for transactional writes.
                  const char  TX_MEM_ABORT_CAPACITY_WRITE[]            = "r0254";
                  //     !   "EventName": "TX_MEM.ABORT_HLE_STORE_TO_ELIDED_LOCK",
                  //     !   "BriefDescription": "Number of times a HLE transactional region aborted due to a non XRELEASE
                  //     !    prefixed instruction writing to an elided lock in the elision buffer.
                  const char  TX_MEM_ABORT_HLE_STORE_TO_ELIDED_LOCK[]  = "r0454";
                  //     !     "EventName": "TX_MEM.ABORT_HLE_ELISION_BUFFER_NOT_EMPTY",
                  //     !   "BriefDescription": "Number of times an HLE transactional execution aborted due to NoAllocatedElisionBuffer being non-zero
                  const char  TX_MEM_ABORT_HLE_ELISION_BUFFER_NOT_EMPTY[]  = "r0854";
                  //     !    "EventName": "TX_MEM.ABORT_HLE_ELISION_BUFFER_MISMATCH",
                  //     !    "BriefDescription": "Number of times an HLE transactional execution aborted due to
                  //     !    XRELEASE lock not satisfying the address and value requirements in the elision buffer.
                  const char  TX_MEM_ABORT_HLE_ELISION_BUFFER_MISMATCH[]  = "r1054";
                  //     !    "EventName": "TX_MEM.ABORT_HLE_ELISION_BUFFER_UNSUPPORTED_ALIGNMENT",
                  //     !    "BriefDescription": "Number of times an HLE transactional execution aborted due to an unsupported read alignment from the elision buffer.
                  const char  TX_MEM_ABORT_HLE_ELISION_BUFFER_UNSUPPORTED_ALIGNMENT[] = "r2054";
                  //     !     "EventName": "TX_MEM.HLE_ELISION_BUFFER_FULL",
                  //     !     "BriefDescription": "Number of times HLE lock could not be elided due to ElisionBufferAvailable being zero.
                  const char  TX_MEM_HLE_ELISION_BUFFER_FULL[]            = "r4054";
		  //     !      "EventName": "MOVE_ELIMINATION.INT_ELIMINATED",
                  //     !      "BriefDescription": "Number of integer Move Elimination candidate uops that were eliminated
                  const char  MOVE_ELIMINATION_INT_ELIMINATED[]           = "r0158";
                  //     !    "EventName": "MOVE_ELIMINATION.SIMD_ELIMINATED",
                  //     !    "BriefDescription": "Number of SIMD Move Elimination candidate uops that were eliminated."
                  const char  MOVE_ELIMINATION_SIMD_ELIMINATED[]          = "r0258";
                  //     !     "EventName": "MOVE_ELIMINATION.INT_NOT_ELIMINATED",
                  //     !     "BriefDescription": "Number of integer Move Elimination candidate uops that were not eliminated.
                  const char  MOVE_ELIMINATION_INT_NOT_ELIMINATED[]       = "r0458";
                  //     !    "EventName": "MOVE_ELIMINATION.SIMD_NOT_ELIMINATED",
                  //     !   "BriefDescription": "Number of SIMD Move Elimination candidate uops that were not eliminated
                  const char  MOVE_ELIMINATION_SIMD_NOT_ELIMINATED[]      = "r0858";
                  //     !    "EventName": "CPL_CYCLES.RING0",
                  //     !   "BriefDescription": "Unhalted core cycles when the thread is in ring 0
                  const char  CPL_CYCLES_RING0[]                          = "r015C";
                  //     !  "EventName": "CPL_CYCLES.RING0_TRANS",
                  //     !  "BriefDescription": "Number of intervals between processor halts while thread is in ring 0
                  const char  CPL_CYCLES_RING0_TRANS[]                    = "r015C";
                  //     !    "EventName": "CPL_CYCLES.RING123",
                  //     !   "BriefDescription": "Unhalted core cycles when thread is in rings 1, 2, or 3
                 const  char CPL_CYCLES_RING123[]                        = "r025C";
                  //     !  "EventName": "TX_EXEC.MISC1",
                  //     !  "BriefDescription": "Counts the number of times a class of instructions that may cause a transactional abort was executed." 
                  //     !  Since this is the count of execution, it may not always cause a transactional abort.
                 const  char  TX_EXEC_MISC1[]                             = "r015D";
                  //     !  "EventName": "TX_EXEC.MISC2",
                  //     !  "BriefDescription": "Counts the number of times a class of instructions (e.g., vzeroupper) that may cause a transactional abort was executed inside a transactional region.
                 const  char  TX_EXEC_MISC2[]                             = "r025D";
                  //     !   "EventName": "TX_EXEC.MISC3",
                  //     !   "BriefDescription": "Counts the number of times an instruction execution caused the transactional nest count supported to be exceeded.
                 const char   TX_EXEC_MISC3[]                             = "r045D";
                  //     !    "EventName": "TX_EXEC.MISC4",
                  //     !  "BriefDescription": "Counts the number of times a XBEGIN instruction was executed inside an HLE transactional region."
                 const char   TX_EXEC_MISC4[]                             = "r085D";
                  //     !   "EventName": "TX_EXEC.MISC5",
                  //     !   "BriefDescription": "Counts the number of times an HLE XACQUIRE instruction was executed inside an RTM transactional region.
                 const char   TX_EXEC_MISC5[]                             = "r105D";
                  //     !    "EventName": "RS_EVENTS.EMPTY_CYCLES",
                  //     !  "BriefDescription": "Cycles when Reservation Station (RS) is empty for the thread",
                  //     !  "PublicDescription": "This event counts cycles when the Reservation Station ( RS )
		  //     is empty for the thread. The RS is a structure that buffers allocated micro-ops from the Front-end.
                  //     !  If there are many cycles when the RS is empty, it may represent an underflow of instructions delivered from the Front-end.
                 const char   RS_EVENTS_EMPTY_CYCLES                      = "r015E";
                  //     !  "EventName": "RS_EVENTS.EMPTY_END",
                  //     !   "BriefDescription": "Counts end of periods where the Reservation Station (RS) was empty. Could be useful to precisely locate Frontend Latency Bound issues.
                 const char   RS_EVENTS_EMPTY_END[]                       = "r015E";
		 //      !   "EventName": "OFFCORE_REQUESTS_OUTSTANDING.DEMAND_DATA_RD",
                 //      !  "BriefDescription": "Offcore outstanding Demand Data Read transactions in uncore queue.",
                 //      !  "PublicDescription": "Offcore outstanding demand data read transactions in SQ to uncore. Set Cmask=1 to count cycles.
                 const char  OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD[] = "r0160";
                 //      !   "EventName": "OFFCORE_REQUESTS_OUTSTANDING.CYCLES_WITH_DEMAND_DATA_RD",
                 //      !  "BriefDescription": "Cycles when offcore outstanding Demand Data Read transactions are present in SuperQueue (SQ), queue to uncore."
                 const char  OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_DATA_RD[] = "r0160";
                 //      !   "EventName": "OFFCORE_REQUESTS_OUTSTANDING.DEMAND_DATA_RD_GE_6",
                 //      !   "BriefDescription": "Cycles with at least 6 offcore outstanding Demand Data Read transactions in uncore queue.
                 const char  OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD_GE_6[]        = "r0160";
                 //      !   "EventName": "OFFCORE_REQUESTS_OUTSTANDING.DEMAND_CODE_RD",
                 //      !  "BriefDescription": "Offcore outstanding code reads transactions in SuperQueue (SQ), queue to uncore, every cycle",
                 //      !  "PublicDescription": "Offcore outstanding Demand code Read transactions in SQ to uncore. Set Cmask=1 to count cycles.
                 const char  OFFCORE_REQUESTS_OUTSTANDING_DEMAND_CODE_R[]              = "r0260";
                 //      !   "EventName": "OFFCORE_REQUESTS_OUTSTANDING.DEMAND_RFO",
                 //      !  "BriefDescription": "Offcore outstanding RFO store transactions in SuperQueue (SQ), queue to uncore",
                 //      !  "PublicDescription": "Offcore outstanding RFO store transactions in SQ to uncore. Set Cmask=1 to count cycles."
                 const char  OFFCORE_REQUESTS_OUTSTANDING_DEMAND_RFO[]                = "r0460";
                 //      !   "EventName": "OFFCORE_REQUESTS_OUTSTANDING.ALL_DATA_RD",
                 //      !  "BriefDescription": "Offcore outstanding cacheable Core Data Read transactions in SuperQueue (SQ), queue to uncore",
                 //      !  "PublicDescription": "Offcore outstanding cacheable data read transactions in SQ to uncore. Set Cmask=1 to count cycles."
                 const char  OFFCORE_REQUESTS_OUTSTANDING_ALL_DATA_RD[]               = "r0860";
                 //      !   "EventName": "OFFCORE_REQUESTS_OUTSTANDING.CYCLES_WITH_DATA_RD",
                 //      !  "BriefDescription": "Cycles when offcore outstanding cacheable Core Data Read transactions are present in SuperQueue (SQ), queue to uncore."
                 const char  OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DATA_RD[]       = "r0860";
                 //      !   "EventName": "LOCK_CYCLES.SPLIT_LOCK_UC_LOCK_DURATION",
                 //      !  "BriefDescription": "Cycles when L1 and L2 are locked due to UC or split lock",
                 //      !  "PublicDescription": "Cycles in which the L1D and L2 are locked, due to a UC lock or split lock.
                 const char  LOCK_CYCLES_SPLIT_LOCK_UC_LOCK_DURATION[]                = "r0163";
                 //      !   "EventName": "LOCK_CYCLES.CACHE_LOCK_DURATION",
                 //      !  "BriefDescription": "Cycles when L1D is locked",
                 //      !  "PublicDescription": "Cycles in which the L1D is locked.
                 const char  LOCK_CYCLES_CACHE_LOCK_DURATION[]                        = "r0263";
                 //      !   "EventName": "IDQ.EMPTY",
                 //      !   "BriefDescription": "Instruction Decode Queue (IDQ) empty cycles"
                 const char   IDQ_EMPTY[]                                              = "r0279";
                 //      !    "EventName": "IDQ.MITE_UOPS",
                 //      !  "BriefDescription": "Uops delivered to Instruction Decode Queue (IDQ) from MITE path",
                 //      !  "PublicDescription": "Increment each cycle # of uops delivered to IDQ from MITE path. Set Cmask = 1 to count cycles.
                 const char  IDQ_MITE_UOPS[]                                          = "r0479";
		 //      !   "EventName": "IDQ.MITE_CYCLES",
                 //      !  "BriefDescription": "Cycles when uops are being delivered to Instruction Decode Queue (IDQ) from MITE path.
                 const char  IDQ_MITE_CYCLES[]                                        = "r0479";
                 //      !   "EventName": "IDQ.DSB_UOPS",
                 //      !  "BriefDescription": "Uops delivered to Instruction Decode Queue (IDQ) from the Decode Stream Buffer (DSB) path",
                 //      !  "PublicDescription": "Increment each cycle. # of uops delivered to IDQ from DSB path. Set Cmask = 1 to count cycles.
                 const char  IDQ_DSB_UOPS[]                                           = "r0879";
                 //      !   "EventName": "IDQ.MS_DSB_UOPS",
                 //      !  "BriefDescription": "Uops initiated by Decode Stream Buffer (DSB) that are being delivered to Instruction Decode Queue (IDQ) while Microcode Sequenser (MS) is busy",
                 //      !  "PublicDescription": "Increment each cycle # of uops delivered to IDQ when MS_busy by DSB. Set Cmask = 1 to count cycles. Add Edge=1 to count # of delivery
                const  char  IDQ_MS_DSB_UOPS[]                                        = "r1079";
                 //      !   "EventName": "IDQ.MS_DSB_CYCLES",
                 //      !  "BriefDescription": "Cycles when uops initiated by Decode Stream Buffer (DSB) are being delivered to Instruction Decode Queue (IDQ) while Microcode Sequenser (MS) is busy
                const  char  IDQ_MS_DSB_CYCLES[]                                      = "r1079";
                //       !    "EventName": "IDQ.MS_DSB_OCCUR",
                //       !   "BriefDescription": "Deliveries to Instruction Decode Queue (IDQ) initiated by Decode Stream Buffer (DSB) while Microcode Sequenser (MS) is busy
                const  char  IDQ_MS_DBS_OCCUR[]                                       = "r1079";
                //       !    "EventName": "IDQ.ALL_DSB_CYCLES_4_UOPS",
                //       !  "BriefDescription": "Cycles Decode Stream Buffer (DSB) is delivering 4 Uops",
                //       !  "PublicDescription": "Counts cycles DSB is delivered four uops. Set Cmask = 4.
                const  char  IDQ_ALL_DBS_CYCLES_4_UOPS[]                              = "r1879";
                //       !   "EventName": "IDQ.ALL_DSB_CYCLES_ANY_UOPS",
                //       !  "BriefDescription": "Cycles Decode Stream Buffer (DSB) is delivering any Uop",
                //       !  "PublicDescription": "Counts cycles DSB is delivered at least one uops. Set Cmask = 1.",
                const char   IDQ_ALL_DBS_CYCLES_ANY_UOPS[]                            = "r1879";
                //       !   "EventName": "IDQ.MS_MITE_UOPS",
                //       !  "BriefDescription": "Uops initiated by MITE and delivered to Instruction Decode Queue (IDQ) while Microcode Sequenser (MS) is busy",
                //       !  "PublicDescription": "Increment each cycle # of uops delivered to IDQ when MS_busy by MITE. Set Cmask = 1 to count cycles."
                const char   IDQ_MS_MITE_UOPS[]                                       = "r2079";
                //       !    "EventName": "IDQ.ALL_MITE_CYCLES_4_UOPS",
                //       !  "BriefDescription": "Cycles MITE is delivering 4 Uops",
                //       !  "PublicDescription": "Counts cycles MITE is delivered four uops. Set Cmask = 4.
                const char   IDQ_ALL_MITE_CYCLES_4_UOPS[]                             = "r2479";
	        //       !   "EventName": "IDQ.ALL_MITE_CYCLES_ANY_UOPS",
                //       !  "BriefDescription": "Cycles MITE is delivering any Uop",
                //       !  "PublicDescription": "Counts cycles MITE is delivered at least one uop. Set Cmask = 1
                const char   IDQ_ALL_MITE_CYCLES_ANY_UOPS[]                           = "r2479";
                //       !  "EventName": "IDQ.MS_UOPS",
                //       !  "BriefDescription": "Uops delivered to Instruction Decode Queue (IDQ) while Microcode Sequenser (MS) is busy",
                //       !   "PublicDescription": "This event counts uops delivered by the Front-end with the assistance of the microcode sequencer.
                //       !  Microcode assists are used for complex instructions or scenarios that can't    
                //       be handled by the standard decoder.  Using other instructions, if possible, will usually improve performance."
                const char   IDQ_MS_UOPS[]                                            = "r3079";
                //       !  "EventName": "IDQ.MS_CYCLES",
                //       !  "BriefDescription": "Cycles when uops are being delivered to Instruction Decode Queue (IDQ) while Microcode Sequenser (MS) is busy",
                //       !  "PublicDescription": "This event counts cycles during which the microcode sequencer assisted the Front-end in delivering uops.
                //       !  Microcode assists are used for complex instructions or scenarios that !
		//       can't be handled by the standard decoder.  Using other instructions, if possible, will usually improve performance.",
                //       !   "Counter": "0,1,2,3",
               const char    IDQ_MS_CYCLES[]                                          = "r3079";
                //       !  "EventName": "IDQ.MS_SWITCHES",
                //       !  "BriefDescription": "Number of switches from DSB (Decode Stream Buffer) or MITE (legacy decode pipeline) to the Microcode Sequencer.",
                //       !  "PublicDescription": "Number of switches from DSB (Decode Stream Buffer) or MITE (legacy decode pipeline) to the Microcode Sequencer."
               const char    IDQ_MS_SWITCHES[]                                        = "r3079";
                //       !  "EventName": "IDQ.MITE_ALL_UOPS",
                //       !   "BriefDescription": "Uops delivered to Instruction Decode Queue (IDQ) from MITE path",
                //       !  "PublicDescription": "Number of uops delivered to IDQ from any path.",
               const char    IDQ_MITE_ALL_UOPS[]                                      = "r3c79";
                //       !  "EventName": "ICACHE.HIT",
                //       !  "BriefDescription": "Number of Instruction Cache, Streaming Buffer and Victim Cache Reads. both cacheable and noncacheable, including UC fetches.
               const char    ICACHE_HIT[]                                             = "r0180";
                //       !   "EventName": "ICACHE.MISSES",
                //       !  "BriefDescription": "Number of Instruction Cache, Streaming Buffer and Victim Cache Misses. Includes Uncacheable accesses.",
                //       !  "PublicDescription": "This event counts Instruction Cache (ICACHE) misses.
               const char    ICACHE_MISSES[]                                          = "r0280";
                //       !  "EventName": "ICACHE.IFETCH_STALL",
                //       !   "BriefDescription": "Cycles where a code fetch is stalled due to L1 instruction-cache miss.",
                //       !   "PublicDescription": "Cycles where a code fetch is stalled due to L1 instruction-cache miss."
               const char    ICACHE_IFETCH_STALL[]                                    = "r0480";
                //       !   "EventName": "ICACHE.IFDATA_STALL",
                //       !  "BriefDescription": "Cycles where a code fetch is stalled due to L1 instruction-cache miss.
               const char    ICACHE_IFDATA_STALL[]                                    = "r0480";
               //        !  "EventName": "ITLB_MISSES.MISS_CAUSES_A_WALK",
               //        !  "BriefDescription": "Misses at all ITLB levels that cause page walks",
               //        !  "PublicDescription": "Misses in ITLB that causes a page walk of any page size.
               const char    ITLB_MISSES_MISS_CAUSES_A_WALK[]                         = "r0185";
               //        !   "EventName": "ITLB_MISSES.WALK_COMPLETED_4K",
               //        !  "BriefDescription": "Code miss in all TLB levels causes a page walk that completes. (4K)"
               const char    ITLB_MISSES_WALK_COMPLETED_4K[]                          = "r0285";
               //        !   "EventName": "ITLB_MISSES.WALK_COMPLETED_2M_4M",
               //        !   "BriefDescription": "Code miss in all TLB levels causes a page walk that completes. (2M/4M)",
               const char    ITLB_MISSES_WALK_COMPLETED_2M_4M[]                       = "r0485";
               //        !   "EventName": "ITLB_MISSES.WALK_COMPLETED_1G",
               //        !   "BriefDescription": "Store miss in all TLB levels causes a page walk that completes. (1G)
               const char    ITLB_MISSES_WALK_COMPLETED_1G[]                          = "r0885";
               //        !   "EventName": "ITLB_MISSES.WALK_COMPLETED",
               //        !   "BriefDescription": "Misses in all ITLB levels that cause completed page walks
               const char    ITLB_MISSES_WALK_COMPLETED[]                             = "r0e85";
               //        !    "EventName": "ITLB_MISSES.WALK_DURATION",
               //        !   "BriefDescription": "Cycles when PMH is busy with page walks
               const char    ITLB_MISSES_WALK_DURATION[]                              = "r1085";
               //        !    "EventName": "ITLB_MISSES.STLB_HIT_4K",
               //        !    "BriefDescription": "Core misses that miss the  DTLB and hit the STLB (4K)"
               const char    ITLB_MISSES_STLB_HIT_4K[]                                = "r2085";
               //        !     "EventName": "ITLB_MISSES.STLB_HIT_2M",
               //        !   "BriefDescription": "Code misses that miss the  DTLB and hit the STLB (2M)"
               const char    ITLB_MISSES_STLB_HIT_2M[]                                = "r4085";
               //        !   "EventName": "ITLB_MISSES.STLB_HIT",
               //        !  "BriefDescription": "Operations that miss the first ITLB level but hit the second and do not cause any page walks"
               const char    ITLB_MISSES_STLB_HIT[]                                   = "r6085";
               //        !   "EventName": "ILD_STALL.LCP",
               //        !  "BriefDescription": "Stalls caused by changing prefix length of the instruction.",
               //        !  "PublicDescription": "This event counts cycles where the decoder is stalled on an instruction with a length changing prefix (LCP).
               const char    ILD_STALL_LCP[]                                          = "r0187";
               //        !   "EventName": "ILD_STALL.IQ_FULL",
               //        !  "BriefDescription": "Stall cycles because IQ is full",
               const char    ILD_STALL_IQ_FULL[]                                      = "r0487";
               //        ! "EventName": "BR_INST_EXEC.NONTAKEN_CONDITIONAL",
               //        !  "BriefDescription": "Not taken macro-conditional branches."
               const char    BR_INST_EXEC_NONTAKEN_CONDITIONAL[]                      = "r4188";
               //        !   "EventName": "BR_INST_EXEC.TAKEN_CONDITIONAL",
               //        !   "BriefDescription": "Taken speculative and retired macro-conditional branches.
               const char    BR_INST_EXEC_TAKEN_CONDITIONAL[]                         = "r8188";
               //        !  "EventName": "BR_INST_EXEC.TAKEN_DIRECT_JUMP",
               //        !  "BriefDescription": "Taken speculative and retired macro-conditional branch instructions excluding calls and indirects.
               const char    BR_INST_EXEC_TAKEN_DIRECT_JUMP[]                         = "r8288";
               //        !   "EventName": "BR_INST_EXEC.TAKEN_INDIRECT_JUMP_NON_CALL_RET",
               //        !  "BriefDescription": "Taken speculative and retired indirect branches excluding calls and returns."
               const char    BR_INST_EXEC_TAKEN_INDIRECT_JUMP_NON_CALL_RET[]          = "r8488";
               //        ! "EventName": "BR_INST_EXEC.TAKEN_INDIRECT_NEAR_RETURN",
               //        !  "BriefDescription": "Taken speculative and retired indirect branches with return mnemonic."
               const char    BR_INST_EXEC_TAKEN_INDIRECT_NEAR_RETURN[]                = "r8888";
               //        !  "EventName": "BR_INST_EXEC.TAKEN_DIRECT_NEAR_CALL",
               //        !  "BriefDescription": "Taken speculative and retired direct near calls."
               const char    BR_INST_EXEC_TAKEN_DIRECT_NEAR_CALL[]                    = "r9088";
               //        !  "EventName": "BR_INST_EXEC.TAKEN_INDIRECT_NEAR_CALL",
               //        !   "BriefDescription": "Taken speculative and retired indirect calls.
               const char    BR_INST_EXEC_TAKEN_INDIRECT_NEAR_CALL[]                  = "rA088";
               //        !  "EventName": "BR_INST_EXEC_ALL_CONDITIONAL",
               //        !  "BriefDescription": "Speculative and retired macro-conditional branches.",
               const char    BR_INST_EXEC_ALL_CONDITIONAL[]                           = "rC188";
               //        !   "EventName": "BR_INST_EXEC_ALL_DIRECT_JMP",
               //        !  "BriefDescription": "Speculative and retired macro-unconditional branches excluding calls and indirects."
               const char    BR_INST_EXEC_ALL_DIRECT_JMP[]                            = "rC288";
               //        !  "EventName": "BR_INST_EXEC_ALL_INDIRECT_JUMP_NON_CALL_RET",
               //        !  "BriefDescription": "Speculative and retired indirect branches excluding calls and returns.
               const char    BR_INST_EXEC_ALL_INDIRECT_JUMP_NON_CALL_RET[]            = "rC488";
               //        !    "EventName": "BR_INST_EXEC_ALL_INDIRECT_NEAR_RETURN",
               //        !  "BriefDescription": "Speculative and retired indirect return branches."
               const char    BR_INST_EXEC_ALL_INDIRECT_NEAR_RETURN[]                  = "rC888";
               //        !  "EventName": "BR_INST_EXEC.ALL_DIRECT_NEAR_CALL",
               //        !  "BriefDescription": "Speculative and retired direct near calls."
               const char    BR_INST_EXEC_ALL_DIRECT_NEAR_CALL[]                      = "rD088";
               //        !  "EventName": "BR_INST_EXEC.ALL_BRANCHES",
               //        !   "BriefDescription": "Speculative and retired  branches",
               const char    BR_INST_EXEC_ALL_BRANCHES[]                              = "rFF88";
               //        !  "EventName": "BR_MISP_EXEC.NONTAKEN_CONDITIONAL",
               //        !  "BriefDescription": "Not taken speculative and retired mispredicted macro conditional branches.
               const char    BR_MISP_EXEC_NONTAKEN_CONDITIONAL[]                      = "r4189";
               //        !   "EventName": "BR_MISP_EXEC.TAKEN_CONDITIONAL",
               //        !   "BriefDescription": "Taken speculative and retired mispredicted macro conditional branches."
               const char    BR_MISP_EXEC_TAKEN_CONDITIONAL[]                         = "r4189";
               //        !    "EventName": "BR_MISP_EXEC.TAKEN_INDIRECT_JUMP_NON_CALL_RET",
               //        !   "BriefDescription": "Taken speculative and retired mispredicted indirect branches excluding calls and returns.
               const char    BR_MISP_EXEC_TAKEN_INDIRECT_JUMP_NON_CALL_RET[]          = "r8489";
               //        !   "EventName": "BR_MISP_EXEC.TAKEN_RETURN_NEAR",
               //        !   "BriefDescription": "Taken speculative and retired mispredicted indirect branches with return mnemonic.
               const char    BR_MISP_EXEC_TAKEN_RETURN_NEAR[]                         = "r8889";
               //        ! "EventName": "BR_MISP_EXEC.TAKEN_INDIRECT_NEAR_CALL",
               //        !  "BriefDescription": "Taken speculative and retired mispredicted indirect calls.
               const char    BR_MISP_EXEC_TAKEN_INDIRECT_NEAR_CALL[]                  = "rA089";
               //        !  "EventName": "BR_MISP_EXEC.ALL_CONDITIONAL",
               //        !  "BriefDescription": "Speculative and retired mispredicted macro conditional branches.
               const char    BR_MISP_EXEC_ALL_CONDITIONAL[]                           = "rC189";
               //        !  "EventName": "BR_MISP_EXEC.ALL_INDIRECT_JUMP_NON_CALL_RET",
               //        !  "BriefDescription": "Mispredicted indirect branches excluding calls and returns.
               const char    BR_MISP_EXEC_ALL_INDIRECT_JUMP_NON_CALL_RET[]            = "rC489";
               //        !   "EventName": "BR_MISP_EXEC.ALL_BRANCHES",
               //        !   "BriefDescription": "Speculative and retired mispredicted macro conditional branches"
               const char    BR_MISP_EXEC_ALL_BRANCHES[]                              = "rFF89";
               //        !   "EventName": "IDQ_UOPS_NOT_DELIVERED.CORE",
               //        !  "BriefDescription": "Uops not delivered to Resource Allocation Table (RAT) per thread when backend of the machine is not stalled",
               //        !  "PublicDescription": "This event count the number of undelivered (unallocated) uops 
               //        from the Front-end to the Resource Allocation Table (RAT) while the Back-end of the processor is not stalled.
               //        !  The  !Front-end can allocate up to 4 uops per cycle so this event can increment 0-4 times per
	       //        cycle depending on the number of unallocated uops. This event is counted on a per-core basis."
               const char    IDQ_UOPS_NOT_DELIVERED_CORE[]                            = "r019C";
              /* !   "EventName": "IDQ_UOPS_NOT_DELIVERED.CYCLES_0_UOPS_DELIV.CORE",
               !  "BriefDescription": "Cycles per thread when 4 or more uops are not delivered to Resource Allocation Table (RAT) when backend of the machine is not stalled",
               !  "PublicDescription": "This event counts the number cycles during which the Front-end allocated exactly zero uops
               !  to the Resource Allocation Table (RAT) while the Back-end of the processor is not stalled.  This event is counted on a per-core basis.",*/
               const char    IDQ_UOPS_NOT_DELIVERED_CYCLES_0_UOPS_DELIV_CORE[]        = "r019C";
              /* !     "EventName": "IDQ_UOPS_NOT_DELIVERED.CYCLES_LE_1_UOP_DELIV.CORE",
               !     "BriefDescription": "Cycles per thread when 3 or more uops are not delivered to
               !     Resource Allocation Table (RAT) when backend of the machine is not stalled.*/
               const char    IDQ_UOPS_NOT_DELIVERED_CYCLES_LE_1_UOP_DELIV_CORE[]      = "r019C";
               /*!   "EventName": "IDQ_UOPS_NOT_DELIVERED.CYCLES_LE_2_UOP_DELIV.CORE",
               !  "BriefDescription": "Cycles with less than 2 uops delivered by the front end."*/
               const char    IDQ_UOPS_NOT_DELIVERED_CYCLES_LE_2_UOP_DELIV_CORE[]      = "r019C";
               /*!   "EventName": "IDQ_UOPS_NOT_DELIVERED.CYCLES_LE_3_UOP_DELIV.CORE",
               !   "BriefDescription": "Cycles with less than 3 uops delivered by the front end.*/
               const char    IDQ_UOPS_NOT_DELIVERED_CYCLES_LE_3_UOP_DELIV_CORE[]      = "r019C";
               /*!   "EventName": "IDQ_UOPS_NOT_DELIVERED.CYCLES_FE_WAS_OK",
               !   "BriefDescription": "Counts cycles FE delivered 4 uops or Resource Allocation Table (RAT) was stalling FE.*/
              const char     IDQ_UOPS_NOT_DELIVERED_CYCLES_FE_WAS_OK[]                = "r019C";
              /* !    "EventName": "UOPS_EXECUTED_PORT.PORT_0",
               !    "BriefDescription": "Cycles per thread when uops are executed in port 0*/
              const char     UOPS_EXECUTED_PORT_PORT_0[]                              = "r01A1";
              /* !    "EventName": "UOPS_EXECUTED_PORT.PORT_0_CORE",
               !    "BriefDescription": "Cycles per core when uops are executed in port 0.*/
              const char     UOPS_EXECUTED_PORT_PORT_0_CORE[]                         = "r01A1";
              /* !     "EventName": "UOPS_DISPATCHED_PORT.PORT_0",
               !     "BriefDescription": "Cycles per thread when uops are executed in port 0.*/
              const char     UOPS_DISPATCHED_PORT_PORT_0[]                            = "r01A1";
              /* !      "EventName": "UOPS_EXECUTED_PORT.PORT_1",
               !     "BriefDescription": "Cycles per thread when uops are executed in port 1"*/
              const char     UOPS_EXECUTED_PORT_PORT_1[]                              = "r02A1";
              /*!     "EventName": "UOPS_EXECUTED_PORT.PORT_1_CORE",
              !     "BriefDescription": "Cycles per core when uops are executed in port 1*/
              const char     UOPS_EXECUTED_PORT_PORT_1_CORE[]                         = "r02A1";
              /*!      "EventName": "UOPS_DISPATCHED_PORT.PORT_1",
              !      "BriefDescription": "Cycles per thread when uops are executed in port 1*/
              const char     UOPS_DISPATCHED_PORT_PORT_1[]                            = "r02A1";
             /* !      "EventName": "UOPS_EXECUTED_PORT.PORT_2",
              !      "BriefDescription": "Cycles per thread when uops are executed in port 2*/
              const char     UOPS_EXECUTED_PORT_PORT_2[]                              = "r04A1";
              /*!       "EventName": "UOPS_EXECUTED_PORT.PORT_2_CORE",
              !       "BriefDescription": "Cycles per core when uops are dispatched to port 2.*/
              const char     UOPS_EXECUTED_PORT_PORT_2_CORE[]                         = "r04A1";
              /*!        "EventName": "UOPS_DISPATCHED_PORT.PORT_2",
              !        "BriefDescription": "Cycles per thread when uops are executed in port 2*/
              const char     UOPS_DISPATCHED_PORT_PORT_2[]                            = "r04A1";
              /*!      "EventName": "UOPS_EXECUTED_PORT.PORT_3",
              !      "BriefDescription": "Cycles per thread when uops are executed in port 3"*/
              const char     UOPS_EXECUTED_PORT_PORT_3[]                              = "r08A1";
              /*!       "EventName": "UOPS_EXECUTED_PORT.PORT_3_CORE",
              !       "BriefDescription": "Cycles per core when uops are dispatched to port 3*/
              const char     UOPS_EXECUTED_PORT_PORT_3_CORE[]                         = "r08A1";
              /*!      "EventName": "UOPS_DISPATCHED_PORT.PORT_3",
              !      "BriefDescription": "Cycles per thread when uops are executed in port 3.*/
              const  char     UOPS_DISPATCHED_PORT_PORT_3[]                            = "r08A1";
              /*!       "EventName": "UOPS_EXECUTED_PORT.PORT_4",
              !       "BriefDescription": "Cycles per thread when uops are executed in port 4"*/
              const  char     UOPS_EXECUTED_PORT_PORT_4[]                              = "r10A1";
              /*!        "EventName": "UOPS_EXECUTED_PORT.PORT_4_CORE",
              !         "BriefDescription": "Cycles per core when uops are executed in port 4.*/
              const  char    UOPS_EXECUTED_PORT_PORT_4_CORE[]                         = "r10A1";
              /*!         "EventName": "UOPS_DISPATCHED_PORT.PORT_4",
              !   "BriefDescription": "Cycles per thread when uops are executed in port 4.*/
              const  char    UOPS_DISPATCHED_PORT_PORT_4[]                            = "r10A1";
              /*!    "EventName": "UOPS_EXECUTED_PORT.PORT_5",
              !    "BriefDescription": "Cycles per thread when uops are executed in port 5"*/
              const  char     UOPS_EXECUTED_PORT_PORT_5[]                              = "r20A1";
              /*!      "EventName": "UOPS_EXECUTED_PORT.PORT_5_CORE",
              !      "BriefDescription": "Cycles per core when uops are executed in port 5."*/
              const  char     UOPS_EXECUTED_PORT_PORT_5_CORE[]                         = "r20A1";
              /*!       "EventName": "UOPS_DISPATCHED_PORT.PORT_5",
              !       "BriefDescription": "Cycles per thread when uops are executed in port 5.*/
              const  char     UOPS_DISPATCHED_PORT_PORT_5[]                            = "r20A1";
              /*!       "EventName": "UOPS_EXECUTED_PORT.PORT_6",
              !       "BriefDescription": "Cycles per thread when uops are executed in port 6*/
              const  char      UOPS_EXECUTED_PORT_PORT_6[]                              = "r40A1";
              /*!    "EventName": "UOPS_EXECUTED_PORT.PORT_6_CORE",
              !    "BriefDescription": "Cycles per core when uops are executed in port 6.*/
              const  char      UOPS_EXECUTED_PORT_PORT_6_CORE[]                         = "r40A1";
              /*!   "EventName": "UOPS_DISPATCHED_PORT.PORT_6",
              !    "BriefDescription": "Cycles per thread when uops are executed in port 6.*/
              const  char      UOPS_DISPATCHED_PORT_PORT_6[]                            = "r40A1";
              /*!     "EventName": "UOPS_EXECUTED_PORT.PORT_7",
              !  "BriefDescription": "Cycles per thread when uops are executed in port 7"*/
              const  char      UOPS_EXECUTED_PORT_PORT_7[]                              = "r80A1";
              /*!    "EventName": "UOPS_EXECUTED_PORT.PORT_7_CORE",
              !    "BriefDescription": "Cycles per core when uops are dispatched to port 7."*/
              const  char      UOPS_EXECUTED_PORT_PORT_7_CORE[]                         = "r80A1";
              /*!     "EventName": "UOPS_DISPATCHED_PORT.PORT_7",
              !    "BriefDescription": "Cycles per thread when uops are executed in port 7.*/
              const  char      UOPS_DISPATCHED_PORT_PORT_7[]                            = "r80A1";
              /*!   "EventName": "RESOURCE_STALLS.ANY",
              !   "BriefDescription": "Resource-related stall cycles"*/
              const  char      RESOURCES_STALLS_ANY[]                                   = "r01A2";
              /*!     "EventName": "RESOURCE_STALLS.RS",
              !  "BriefDescription": "Cycles stalled due to no eligible RS entry available."*/
              const  char      RESOURCE_STALLS_RS[]                                     = "r04A2";
              /*!   "EventName": "RESOURCE_STALLS.SB",
              !   "BriefDescription": "Cycles stalled due to no store buffers available. (not including draining form sync).",
              !   "PublicDescription": "This event counts cycles during which no instructions were allocated because no Store Buffers (SB) were available.*/
              const  char       RESOURCE_STALLS_SB[]                                     = "r08A2";
              /*!    "EventName": "RESOURCE_STALLS.ROB",
              !  "BriefDescription": "Cycles stalled due to re-order buffer full.",*/
              const  char       RESOURCE_STALLS_ROB[]                                    = "r10A2";
              /*!    "EventName": "CYCLE_ACTIVITY.CYCLES_L2_PENDING",
              !  "BriefDescription": "Cycles with pending L2 cache miss loads.",
              !  "PublicDescription": "Cycles with pending L2 miss loads. Set Cmask=2 to count cycle."*/
             const   char      CYCLE_ACTIVITIES_CYCLES_L2_PENDING[]                     = "r01A3";
             /* !   "EventName": "CYCLE_ACTIVITY.CYCLES_LDM_PENDING",
              !  "BriefDescription": "Cycles with pending memory loads.",
              !  "PublicDescription": "Cycles with pending memory loads. Set Cmask=2 to count cycle.*/
             const   char      CYCLES_ACTIVITIES_CYCLES_LDM_PENDING[]                   = "r02A3";
              /*!    "EventName": "CYCLE_ACTIVITY.CYCLES_NO_EXECUTE",
              !  "BriefDescription": "This event increments by 1 for every cycle where there was no execute for this thread.",
              !   "PublicDescription": "This event counts cycles during which no instructions were executed in the execution stage of the pipeline.*/
             const   char      CYCLE_ACTIVITY_CYCLES_NO_EXECUTE[]                       = "r04A3";
              /*!    "EventName": "CYCLE_ACTIVITY.STALLS_L2_PENDING",
              !    "BriefDescription": "Execution stalls due to L2 cache misses."*/
             const   char      CYCLE_ACTIVITY_STALLS_L2_PENDING[]                       = "r05A3";
              /*!     "EventName": "CYCLE_ACTIVITY.STALLS_LDM_PENDING",
              !  "BriefDescription": "Execution stalls due to memory subsystem.",
              !  "PublicDescription": "This event counts cycles during which no instructions were executed in the execution stage of
              !   the pipeline and there were memory instructions pending (waiting for data).*/
             const   char      CYCLE_ACTIVITY_STALLS_LDM_PENDING[]                      = "r06A3";
              /*!    "EventName": "CYCLE_ACTIVITY.CYCLES_L1D_PENDING",
              !   "BriefDescription": "Cycles with pending L1 cache miss loads.",
              !  "PublicDescription": "Cycles with pending L1 data cache miss loads. Set Cmask=8 to count cycle."*/
             const   char      CYCLE_ACTIVITY_CYCLES_L1D_PENDING[]                      = "r08A3";
              /*!    "EventName": "CYCLE_ACTIVITY.STALLS_L1D_PENDING",
              !  "BriefDescription": "Execution stalls due to L1 data cache misses",
              !  "PublicDescription": "Execution stalls due to L1 data cache miss loads. Set Cmask=0CH.*/
             const   char      CYCLE_ACTIVITY_STALLS_L1D_PENDING[]                      = "r0CA3";
              /*!    "EventName": "LSD.UOPS",
              !   "BriefDescription": "Number of Uops delivered by the LSD.",*/
             const   char      LSD_UOPS[]                                               = "r01A8";
              /*!     "EventName": "LSD.CYCLES_ACTIVE",
              !   "BriefDescription": "Cycles Uops delivered by the LSD, but didn't come from the decoder.",*/
            const    char      LSD_CYCLES_ACTIVE[]                                      = "r01A8";
              /*!    "EventName": "LSD.CYCLES_4_UOPS",
              !   "BriefDescription": "Cycles 4 Uops delivered by the LSD, but didn't come from the decoder.*/
            const    char      LSD_CYCLES_4_UOPS[]                                      = "r01A8";
              /*!    "EventName": "DSB2MITE_SWITCHES.PENALTY_CYCLES",
              !    "BriefDescription": "Decode Stream Buffer (DSB)-to-MITE switch true penalty cycles."*/
            const    char      DSB2MITE_SWITCHES_PENALTY_CYCLES[]                       = "r02AB";
              /*!     "EventName": "ITLB.ITLB_FLUSH",
              !     "BriefDescription": "Flushing of the Instruction TLB (ITLB) pages, includes 4k/2M/4M pages.*/
            const    char      ITLB_ITLB_FLUSH[]                                        = "r01AE";
              /*!      "EventName": "OFFCORE_REQUESTS.DEMAND_DATA_RD",
              !      "BriefDescription": "Demand Data Read requests sent to uncore",*/
            const    char      OFFCORE_REQUESTS_DEMAND_DATA_RD[]                        = "r01B0";
              /*!       "EventName": "OFFCORE_REQUESTS.DEMAND_CODE_RD",
              !       "BriefDescription": "Cacheable and noncachaeble code read requests",*/
            const    char      OFFCORE_REQUESTS_DEMAND_CODE_RD[]                        = "r02B0";
              /*!        "EventName": "OFFCORE_REQUESTS.DEMAND_RFO",
              !  "BriefDescription": "Demand RFO requests including regular RFOs, locks, ItoM",
              !  "PublicDescription": "Demand RFO read requests sent to uncore, including regular RFOs, locks, ItoM.",*/
            const    char      OFFCORE_REQUESTS_DEMAND_RFO[]                            = "r04B0";
              /*!    "EventName": "OFFCORE_REQUESTS.ALL_DATA_RD",
              !  "BriefDescription": "Demand and prefetch data reads",
              !  "PublicDescription": "Data read requests sent to uncore (demand and prefetch).*/
            const     char      OFFCORE_REQUESTS_ALL_DATA_RD[]                           = "r08B0";
              /*!   "EventName": "UOPS_EXECUTED.STALL_CYCLES",
              !  "BriefDescription": "Counts number of cycles no uops were dispatched to be executed on this thread."*/
            const     char      UOPS_EXECUTED_STALL_CYCLES[]                             = "r01B1";
              /*!   "EventName": "UOPS_EXECUTED.CYCLES_GE_1_UOP_EXEC",
              !  "BriefDescription": "Cycles where at least 1 uop was executed per-thread",
              !  "PublicDescription": "This events counts the cycles where at least one uop was executed. It is counted per thread.",*/
            const     char      UOPS_EXECUTED_CYCLES_GE_1_UOP_EXEC[]                     = "r01B1";
              /*!  "EventName": "UOPS_EXECUTED.CYCLES_GE_2_UOPS_EXEC",
              !  "BriefDescription": "Cycles where at least 2 uops were executed per-thread",
              !  "PublicDescription": "This events counts the cycles where at least two uop were executed. It is counted per thread.",*/
            const     char      UOPS_EXECUTED_CYCLES_GE_2_UOP_EXEC[]                     = "r01B1";
              /*!   "EventName": "UOPS_EXECUTED.CYCLES_GE_3_UOPS_EXEC",
              !  "BriefDescription": "Cycles where at least 3 uops were executed per-thread",
              !  "PublicDescription": "This events counts the cycles where at least three uop were executed. It is counted per thread.",*/
            const     char      UOPS_EXECUTED_CYCLES_GE_3_UOPS_EXEC[]                     = "r01B1";
              /*!   "EventName": "UOPS_EXECUTED.CYCLES_GE_4_UOPS_EXEC",
              !  "BriefDescription": "Cycles where at least 4 uops were executed per-thread.",
              !  "PublicDescription": "Cycles where at least 4 uops were executed per-thread.",*/
            const     char      UOPS_EXECUTED_CYCLES_GE_4_UOPS_EXEC[]                     = "r01B1";
              /*!   "EventName": "UOPS_EXECUTED.CORE",
              !  "BriefDescription": "Number of uops executed on the core.",
              !  "PublicDescription": "Counts total number of uops to be executed per-core each cycle.",*/
            const    char       UOPS_EXECUTED_CORE[]                                      = "r02B1";
              /*!   "EventName": "UOPS_EXECUTED.CORE_CYCLES_GE_1",
              !   "BriefDescription": "Cycles at least 1 micro-op is executed from any thread on physical core.",*/
            const    char      UOPS_EXECUTED_CORE_CYCLES_GE_1[]                          = "r02B1";
              /*!    "EventName": "UOPS_EXECUTED.CORE_CYCLES_GE_2",
              !    "BriefDescription": "Cycles at least 2 micro-op is executed from any thread on physical core.",*/
            const    char     UOPS_EXECUTED_CORE_CYCLES_GE_2[]                          = "r02B1";
              /*!     "EventName": "UOPS_EXECUTED.CORE_CYCLES_GE_3",
              !   "BriefDescription": "Cycles at least 3 micro-op is executed from any thread on physical core.",*/
            const    char     UOPS_EXECUTED_CORE_CYCLES_GE_3[]                         = "r02B1";
              /*!    "EventName": "UOPS_EXECUTED.CORE_CYCLES_GE_4",
              !   "BriefDescription": "Cycles at least 4 micro-op is executed from any thread on physical core.",*/
            const    char     UOPS_EXECUTED_CORE_CYCLES_GE_4[]                        = "r02B1";
              /*!    "EventName": "UOPS_EXECUTED.CORE_CYCLES_NONE",
              !    "BriefDescription": "Cycles with no micro-ops executed from any thread on physical core.",*/
            const    char     UOPS_EXECUTED_CORE_CYCLES_NONE[]                         = "r02B1";
              /*!    "EventName": "OFFCORE_REQUESTS_BUFFER.SQ_FULL",
              !   "BriefDescription": "Offcore requests buffer cannot take more entries for this thread core.",*/
            const    char     OFFCORE_REQUESTS_BUFFER_SQ_FULL[]                        = "r01b2";
              /*!   "EventName": "OFFCORE_RESPONSE",
              !   "BriefDescription": "Offcore response can be programmed only with a specific pair of event select and counter MSR,
              !and with specific event codes and predefine mask bit value in a dedicated MSR to specify attributes of the offcore transaction.",*/
           const     char     OFFCORE_RESPONSE[]                                       = "r0B7";
             /* !  "EventName": "PAGE_WALKER_LOADS.DTLB_L1",
             ! "BriefDescription": "Number of DTLB page walker hits in the L1+FB"*/
           const     char     PAGE_WALKER_LOADS_DTLB_L1[]                              = "r11BC";
             /*  !   "EventName": "PAGE_WALKER_LOADS.DTLB_L2",
               ! "BriefDescription": "Number of DTLB page walker hits in the L2",*/
           const     char     PAGE_WALKER_LOADS_DTLB_L2[]                              = "r12BC";
             /*  !    "EventName": "PAGE_WALKER_LOADS.DTLB_L3",
               !  BriefDescription": "Number of DTLB page walker hits in the L3 + XSNP",*/
           const     char     PAGE_WALKER_LOADS_DTLB_L3[]                              = "r14BC";
             /*  !   "EventName": "PAGE_WALKER_LOADS.DTLB_MEMORY",
               ! "BriefDescription": "Number of DTLB page walker hits in Memory",*/
           const     char     PAGE_WALKER_LOADS_DTLB_MEMORY[]                          = "r18BC";
             /*  !  "EventName": "PAGE_WALKER_LOADS.ITLB_L1",
               !  "BriefDescription": "Number of ITLB page walker hits in the L1+FB",*/
           const    char      PAGE_WALKER_LOADS_ITLB_L1[]                              = "r21BC";
             /*  !    "EventName": "PAGE_WALKER_LOADS.ITLB_L2",
               !  "BriefDescription": "Number of ITLB page walker hits in the L2",*/
           const    char      PAGE_WALKER_LOADS_ITLB_L2[]                             = "r22BC";
             /*  !    "EventName": "PAGE_WALKER_LOADS.ITLB_L3",
               !  "BriefDescription": "Number of ITLB page walker hits in the L3 + XSNP"*/
           const    char      PAGE_WALKER_LOADS_ITLB_L3[]                              = "r24BC";
             /*  !  "EventName": "PAGE_WALKER_LOADS.ITLB_MEMORY",
               !  "BriefDescription": "Number of ITLB page walker hits in Memory",*/
           const    char      PAGE_WALKER_LOADS_ITLB_MEMORY[]                          = "r28BC";
            /*   !  "EventName": "PAGE_WALKER_LOADS.EPT_DTLB_L1",
               !  "BriefDescription": "Counts the number of Extended Page Table walks from
               !  the DTLB that hit in the L1 and FB."*/
           const    char      PAGE_WALKER_LOADS_EPT_DTLB_L1[]                          = "r41BC";
            /*   !   "EventName": "PAGE_WALKER_LOADS.EPT_DTLB_L2",
               !  "BriefDescription": "Counts the number of Extended Page Table walks
               !  from the DTLB that hit in the L2.*/
           const    char      PAGE_WALKER_LOADS_EPT_DTLB_L2[]                          = "r42BC";
            /*   !   "EventName": "PAGE_WALKER_LOADS.EPT_DTLB_L3",
               !   "BriefDescription": "Counts the number of Extended Page Table walks
               !   from the DTLB that hit in the L3.",*/
           const    char      PAGE_WALKER_LOADS_EPT_DTLB_L3[]                          = "r44BC";
             /*  !   "EventName": "PAGE_WALKER_LOADS.EPT_DTLB_MEMORY",
               !  "BriefDescription": "Counts the number of Extended Page Table walks
               !   from the DTLB that hit in memory.*/
           const    char      PAGE_WALKER_LOADS_EPT_DTLB_MEMORY[]                      = "r48BC";
              /* !    "EventName": "PAGE_WALKER_LOADS.EPT_ITLB_L1",
               !  BriefDescription": "Counts the number of Extended Page Table
               !waks from the ITLB that hit in the L1 and FB."*/
           const    char      PAGE_WALKER_LOADS_EPT_ITLB_L1[]                          = "r81BC";
              /* ! "EventName": "PAGE_WALKER_LOADS.EPT_ITLB_L2",
               !  "BriefDescription": "Counts the number of Extended Page Table walks
               !  from the ITLB that hit in the L2."*/
          const     char      PAGE_WALKER_LOADS_EPT_ITLB_L2[]                          = "r82BC";
              /* !  "EventName": "PAGE_WALKER_LOADS.EPT_ITLB_L3",
               !  "BriefDescription": "Counts the number of Extended Page Table walks
               !  from the ITLB that hit in the L2.",*/
          const     char      PAGE_WALKER_LOADS_EPT_ITLB_L3[]                          = "r84BC";
              /* !   "EventName": "PAGE_WALKER_LOADS.EPT_ITLB_MEMORY",
               !  "BriefDescription": "Counts the number of Extended Page Table walks
               !   from the ITLB that hit in memory."*/
          const      char      PAGE_WALKER_LOADS_EPT_ITLB_MEMORY[]                      = "r88BC";
              /* !   "EventName": "TLB_FLUSH.DTLB_THREAD",
               !  "BriefDescription": "DTLB flush attempts of the thread-specific entries",*/
          const      char      TLB_FLUSH_DTLB_THREAD[]                                  = "r01BD";
             /*  !     "EventName": "TLB_FLUSH.STLB_ANY",
               !  "BriefDescription": "STLB flush attempts",*/
          const      char      TLB_FLUSH_STLB_ANY[]                                     = "r20BD";
              /* !  "EventName": "INST_RETIRED.ANY_P",
               !  "BriefDescription": "Number of instructions retired.
               ! General Counter   - architectural event",*/
          const      char      INST_RETIRED_ANY_P[]                                    = "r00C0";
             /*  !   "EventName": "INST_RETIRED.PREC_DIST",
               !  "BriefDescription": "Precise instruction retired event with HW
               !  to reduce effect of PEBS shadow in IP distribution",*/
          const      char      INT_RETIRED_PREC_DIST[]                                 = "r01C0";
            /*   !   "EventName": "INST_RETIRED.X87",
               !  "BriefDescription": "FP operations retired. X87 FP operations that have no exceptions: Counts also flows that have several X87 or flows that use X87 uops in the exception handling.",
               !  "PublicDescription": "This is a precise version (that is, uses PEBS) of the event that counts FP operations retired.
               !  For X87 FP operations that have no exceptions counting also includes flows that have several X87, or flows that use X87 uops in the exception handling.",*/
          const      char      INST_RETIRED_X87[]                                      = "r02C0";
             /*  !   "EventName": "OTHER_ASSISTS.AVX_TO_SSE",
               !   "BriefDescription": "Number of transitions from AVX-256 to legacy SSE when penalty applicable"*/
          const      char      OTHER_ASSISTS_AVX_TO_SSE[]                              = "r08C1";
             /*  !    "EventName": "OTHER_ASSISTS.SSE_TO_AVX",
               !  "BriefDescription": "Number of transitions from legacy SSE to AVX-256 when penalty applicable ",*/
          const      char      OTHER_ASSISTS_SSE_TO_AVX[]                              = "r10C1";
              /* !    "EventName": "OTHER_ASSISTS.ANY_WB_ASSIST"*/
          const      char      OTHER_ASSISTS_ANY_WB_ASSIST[]                           = "r40C1";
             /*  !   "EventName": "UOPS_RETIRED.ALL",
               !  "BriefDescription": "Actually retired uops.*/
          const      char      UOPS_RETIRED_ALL[]                                      = "r01C2";
             /*  !  "EventName": "UOPS_RETIRED.STALL_CYCLES",
               !   "BriefDescription": "Cycles no executable uops retired",*/
          const       char     UOPS_RETIRED_STALL_CYCLES[]                             = "r01C2";
             /*  !    "EventName": "UOPS_RETIRED.TOTAL_CYCLES",
               !  "BriefDescription": "Number of cycles using always true condition applied to  PEBS uops retired event."*/
          const       char     UOPS_RETIRED_TOTAL_CYCLES[]                             = "r01C2";
              /* !   "EventName": "UOPS_RETIRED.CORE_STALL_CYCLES",
               !  "BriefDescription": "Cycles no executable uops retired on core",*/
          const       char     UOPS_RETIRED_CORE_STALL_CYLES[]                         = "r01C2";
              /* !   "EventName": "UOPS_RETIRED.RETIRE_SLOTS",
               !  "BriefDescription": "Retirement slots used.*/
          const       char     UOPS_RETIRED_RETIRE_SLOTS[]                             = "r02C2";
              /* !   "EventName": "MACHINE_CLEARS.CYCLES",
               !  "BriefDescription": "Cycles there was a Nuke. Account for both thread-specific and All Thread Nukes."*/
          const       char     MACHINE_CLEARS_CYCLES[]                                 = "r01C3";
             /*  !   "EventName": "MACHINE_CLEARS.COUNT",
               !   "BriefDescription": "Number of machine clears (nukes) of any type."*/
          const       char     MACHINE_CLEARS_COUNT[]                                  = "r01C3";
             /*  !    "EventName": "MACHINE_CLEARS.MEMORY_ORDERING",
               !  "BriefDescription": "Counts the number of machine clears due to memory order conflicts.",
               !  "PublicDescription": "This event counts the number of memory ordering machine clears detected.
               !  Memory ordering machine clears can result from memory address aliasing or snoops from another hardware thread or core to data inflight in the pipeline.
               !  Machine clears can have a significant performance impact if they are happening frequently.",*/
          const       char     MACHINE_CLEARS_MEMORY_ORDERING[]                        = "r02C3";
             /*  !   "EventName": "MACHINE_CLEARS.SMC",
               !   "BriefDescription": "Self-modifying code (SMC) detected.*/
          const       char      MACHINE_CLEARS_SMC[]                                    = "r04C3";
             /*  !    "EventName": "MACHINE_CLEARS.MASKMOV",
               !  "BriefDescription": "This event counts the number of executed Intel AVX masked load operations that refer to an
               !   illegal address range with the mask bits set to 0.*/
          const       char      MACHINE_CLEARS_MASKMOV[]                                = "r20C3";
             /*  !    "EventName": "BR_INST_RETIRED.ALL_BRANCHES",
               !   "BriefDescription": "All (macro) branch instructions retired.",*/
          const       char      BR_INST_RETIRED_ALL_BRANCHES[]                          = "r00C4";
              /* !    "EventName": "BR_INST_RETIRED.CONDITIONAL",
               !   "BriefDescription": "Conditional branch instructions retired.",*/
          const       char      BR_INST_RETIRED_CONDITIONAL[]                           = "r01C4";
              /* !    "EventName": "BR_INST_RETIRED.NEAR_CALL",
               !   "BriefDescription": "Direct and indirect near call instructions retired.*/
          const       char      BR_INST_RETIRED_NEAR_CALL[]                             = "r02C4";
              /* !   "EventName": "BR_INST_RETIRED.NEAR_CALL_R3",
               !  "BriefDescription": "Direct and indirect macro near call instructions retired (captured in ring 3)*/
          const       char      BR_INST_RETIRED_NEAR_CALL_R3[]                          = "r02C4";
              /* !  "EventName": "BR_INST_RETIRED.ALL_BRANCHES_PEBS",
               !   "BriefDescription": "All (macro) branch instructions retired.",*/
          const       char      BR_INST_RETIRED_ALL_BRANCHES_PEBS[]                     = "r04C4";
            /*   !   "EventName": "BR_INST_RETIRED.NEAR_RETURN",
               !  "BriefDescription": "Return instructions retired.*/
          const       char      BR_INST_RETIRED_NEAR_RETURN[]                           = "r08C4";
             /*  !    "EventName": "BR_INST_RETIRED.NOT_TAKEN",
               !    "BriefDescription": "Counts all not taken macro branch instructions retired.",*/
          const       char      BR_INST_RETIRED_NOT_TAKEN[]                             = "r10C4";
             /*  !     "EventName": "BR_INST_RETIRED.NEAR_TAKEN",
               !    "BriefDescription": "Taken branch instructions retired.",*/
          const       char      BR_INST_RETIRED_NEAR_TAKEN[]                            = "r20C4";
              /* !   "EventName": "BR_INST_RETIRED.FAR_BRANCH",
               !  "BriefDescription": "Counts the number of far branch instructions retired.",*/
          const        char     BR_INST_RETIRED_FAR_BRANCH[]                            = "r40C4";
             /*  !   "EventName": "BR_MISP_RETIRED.ALL_BRANCHES",
               !  "BriefDescription": "All mispredicted macro branch instructions retired.*/
          const        char     BR_MISP_RETIRED_ALL_BRANCHES[]                          = "r00C5";
             /*  !     "EventName": "BR_MISP_RETIRED.CONDITIONAL",
               !  "BriefDescription": "Mispredicted conditional branch instructions retired.*/
          const        char      BR_MISP_RETIRED_CONDITIONAL[]                           = "r01C5";
             /*  !    "EventName": "BR_MISP_RETIRED.ALL_BRANCHES_PEBS",
               !  "BriefDescription": "Mispredicted macro branch instructions retired. ",*/
          const        char      BR_MISP_RETIRED_ALL_BRANCHES_PEBS[]                     = "r04C5";
             /*  !   "EventName": "BR_MISP_RETIRED.NEAR_TAKEN",
               !  "BriefDescription": "number of near branch instructions retired that were mispredicted and taken.",*/
          const        char      BR_MISP_RETIRED_NEAR_TAKEN[]                            = "r20C5";
             /*  !   "EventName": "AVX_INSTS.ALL",
               !  "BriefDescription": "Approximate counts of AVX & AVX2 256-bit instructions, including non-arithmetic instructions, loads, and stores.  
               ! May count non-AVX instructions that employ 256-bit operations,
	       !including (but not necessarily limited to) rep string instructions that use 256-bit loads and stores for optimized performance,
               ! XSAVE* and XRSTOR*, and operations that transition the x87 FPU data registers between x87 and MMX.",*/
          const        char  AVX_INSTS_ALL[]                                         = "r07C6";
                                                    
              /* !     "EventName": "HLE_RETIRED.START",
               !  "BriefDescription": "Number of times an HLE execution started.",*/
          const        char   HLE_RETIRED_START[]                                     = "r01C8";
              /* !   "EventName": "HLE_RETIRED.COMMIT",
               !  "BriefDescription": "Number of times an HLE execution successfully committed.",*/
          const        char   HLE_RETIRED_COMMIT[]                                    = "r02C8";
             /*  !  "EventName": "HLE_RETIRED.ABORTED",
               !  "BriefDescription": "Number of times an HLE execution aborted due to any reasons (multiple categories may count as one)."*/
          const        char HLE_RETIRED_ABORTED[]                                   = "r04C8";
             /*  !   "EventName": "HLE_RETIRED.ABORTED_MISC1",
               !  "BriefDescription": "Number of times an HLE execution aborted due to various memory events (e.g., read/write capacity and conflicts)*/
          const        char  HLE_RETIRED_ABORTED_MISC1[]                             = "r08C8";
             /*  !  "EventName": "HLE_RETIRED.ABORTED_MISC2",
               !  "BriefDescription": "Number of times an HLE execution aborted due to uncommon conditions."*/
          const        char  HLE_RETIRED_ABORTED_MISC2[]                             = "r10C8";
             /*  !   "EventName": "HLE_RETIRED.ABORTED_MISC3",
               !  "BriefDescription": "Number of times an HLE execution aborted due to HLE-unfriendly instructions.",*/
          const        char   HLE_RETIRED_ABORTED_MISC3[]                             = "r20C8";
             /*  !    "EventName": "HLE_RETIRED.ABORTED_MISC4",
               !  "BriefDescription": "Number of times an HLE execution aborted due to incompatible memory type*/
          const        char    HLE_RETIRED_ABORTED_MISC4[]                             = "r40C8";
             /*  !  "EventName": "FP_ASSIST.X87_OUTPUT",
               !  "BriefDescription": "output - Numeric Overflow, Numeric Underflow, Inexact Result ",*/
          const        char    FP_ASSIST_X87_OUTPUT[]                                  = "r02CA";
             /*  !   "EventName": "FP_ASSIST.X87_INPUT",
               !  "BriefDescription": "input - Invalid Operation, Denormal Operand, SNaN Operand ",*/
          const       char     FP_ASSIST_X87_INPUT[]                                   = "r04CA";
             /*  !    "EventName": "FP_ASSIST.SIMD_OUTPUT",
               !  "BriefDescription": "SSE* FP micro-code assist when output value is invalid.*/
          const       char     FP_ASSIST_SIMD_OUTPUT[]                                 = "r08CA";
              /* !    "EventName": "FP_ASSIST.SIMD_OUTPUT",
               !  "BriefDescription": "SSE* FP micro-code assist when output value is invalid. ",*/
          const       char     FP_ASSIST_SIMD_OUTPUT[]                                 = "r08CA";
             /*  !    "EventName": "FP_ASSIST.SIMD_INPUT",
               !   "BriefDescription": "Any input SSE* FP Assist ",*/
          const       char     FP_ASSIST_SIMD_INPUT[]                                  = "r10CA";
             /*  !    "EventName": "FP_ASSIST.ANY",
               !  "BriefDescription": "Counts any FP_ASSIST umask was incrementing ",*/
          const       char     FP_ASSIST_ANY[]                                         = "r1ECA";
            /*   !     "EventName": "ROB_MISC_EVENTS.LBR_INSERTS",
               !   "BriefDescription": "Count cases of saving new LBR",*/
          const       char     ROB_MISC_EVENTS_LBR_INSERTS[]                           = "r20CC";
              /* !    "EventName": "MEM_TRANS_RETIRED.LOAD_LATENCY_GT_4",
               !  "BriefDescription": "Randomly selected loads with latency value being above 4*/
          const       char     MEM_TRANS_RETIRED_LOAD_LATENCY_GT_4[]                   = "r01CD";
	   };
#else
#error Include this file only in Linux GMS (Guided-Missile-Simulation-Modeling) projects!!
#endif
      }

}





#endif /*__GMS_HASWELL_MARCH_EVENTS_H__*/
