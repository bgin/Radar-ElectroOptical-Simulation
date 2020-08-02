
#ifndef __GMS_PMC_EVENTS_H__
#define __GMS_PMC_EVENTS_H_
















static const char * const LD_BLOCKS[2] = {
    "ld_blocks.store_forward",
    "ld_blocks.no_sr"
};

static const char * const MISALIGN_MEM_REF[2] = {
    "misalign_mem_ref.loads",
    "misalign_mem_ref.stores"
};

static const char * const  LD_BLOCKS_PARTIAL[1] = {
    "ld_blocks_partial.address_alias"
};

static const char * const DTLB_LOAD_MISSES_1[4] = {
    "dtlb_load_misses.misses_causes_a_walk",
    "dtlb_load_misses.walk_completed_4k",
    "dtlb_load_misses.walk_completed_2m_4m",
    "dtlb_load_misses.walk_completed"
   
};

static const char * const DTLB_LOAD_MISSES_2[5] = {
     "dtlb_load_misses.walk_duration",
     "dtlb_load_misses.stlb_hit_4k",
     "dtlb_load_misses.stlb_hit_2m",
     "dtlb_load_misses.stlb_hit",
     "dtlb_load_misses.pde_cache_miss"
};

static const char * const INT_MISC[1] = {
     "int_misc.recovery_cycles"
};

static const char * const UOPS_ISSUED[4] = {
     "uops_issued.any",
     "uops_issued.flags_merge",
     "uops_issued.slow_lea",
     "uops_issued.single_mul"
};

static const char * const L2_RQSTS_1[4] = {
     "l2_rqsts.demand_data_rd_miss",
     "l2_rqsts.demand_data_rd_hit",
     "l2_rqsts.all_demand_data_rd",
     "l2_rqsts.crfo_hit"
};

static const char * const L2_RQSTS_2[4] = {
     "l2_rqsts.rfo_miss",
     "l2_rqsts.all_rfo",
     "l2_rqsts.code_rd_hit",
     "l2_rqsts.code_rd_miss"
};

static const char * const L2_RQSTS_3[4] = {
     "l2_rqsts.all_demand_miss",
     "l2_rqsts.all_demand_references",
     "l2_rqsts.all_code_rd",
     "l2_rqsts.l2_pf_hit"
};

static const char * const L2_RQSTS_4[4] = {
     "l2_rqsts.l2_pf_hit",
     "l2_rqsts.l2_pf_miss",
     "l2_rqsts.all_pf",
     "l2_rqsts.miss"
};

static const char * const L2_RQSTS_5[1] = {
     "l2_rqsts.references"
};

static const char * const L2_DEMAND_RQSTS[1] = {
     "l2_demand_rqsts.wb_hit"
};

static const char * const LONGEST_LAT_CACHE[2] = {
      "longest_lat_cache.reference",
      "longest_lat_cache.miss"
};

static const char * const CPU_CLK_UNHALTED[2] = {
      "cpu_clk_unhalted.core_clk",
      "cpu_clk_unhalted.ref_xclk"
};

static const char * const L1_PEND_MISS[1] = {
      "l1_pend_miss.pending"
};

static const char * const DTLB_STORE_MISSES_1[4] = {
       "dtlb_store_misses.miss_causes_a_walk",
       "dtlb_store_misses.walk_completed_4k",
       "dtlb_store_misses.walk_completed_2m_4m",
       "dtlb_store_misses.walk_completed"
};

static const char * const DTLB_STORE_MISSES_2[4] = {
       "dtlb_store_misses.walk_duration",
       "dtlb_store_misses.stlb_hit_4k",
       "dtlb_store_misses.stlb_hit_2m",
       "dtlb_store_misses.stlb_hit"
};

static const char * const DTLB_STORE_MISSES_3[1] = {
       "dtlb_store_misses.pde_cache_misses"
};

static const char * const LOAD_HIT_PRE[2] = {
       "load_hit_pre.sw_pf",
       "load_hit_pre.hw_pf"
};

static const char * const L1D[1] = {
       "l1d.replacement"
};

static const char * const TX_MEM_1[4] = {
       "tx_mem.abort_conflict",
       "tx_mem.abort_capacity_write",
       "tx_mem.abort_hle_store_to_elided_lock",
       "tx_mem.abort_hle_elision_buffer_not_empty"
       
};


static const char * const TX_MEM_2[3] = {
       "tx_mem.abort_hle_elision_buffer_mismatch",
       "tx_mem.abort_hle_elision_buffer_unsupported_alignment",
       "tx_mem.hle_elision_buffer_full"
};

static const char * const MOVE_ELIMINATION[4] = {
       "move_elimination.int_not_eliminated",
       "move_elimination.simd_not_eliminated",
       "move_elimination.int_eliminated",
       "move_elimination.simd_eliminated"
};

static const char * const CPL_CYCLES[2] = {
       "cpl_cycles.ring0",
       "cpl_cycles.ring123"
};

static const char * const TX_EXEC_1[4] = {
       "tx_exec.misc1",
       "tx_exec.misc2",
       "tx_exec.misc3",
       "tx_exec.misc4"
};

static const char * const TX_EXEC_2[1] = {
       "tx_exec.misc5"
};

static const char * const RS_EVENTS[1] = {
       "rs_events.empty_cycles"
};

static const char * const OFFCORE_REQUESTS_OUTSTANDING[4] = {
       "offcore_requests_outstanding.demand_data_rd",
       "offcore_requests_outstanding.demand_code_rd",
       "offcore_requests_outstanding.demand_rfo",
       "offcore_requests_outstanding.all_data_rd"
};

static const char * const LOCK_CYCLES[2] = {
       "lock_cycles.split_lock_uc_lock_duration",
       "lock_cycles.cache_lock_duration"
};

static const char * const IDQ_1[4] = {
       "idq.empty",
       "idq.mite_uops",
       "idq.dsb_uops",
       "idq.ms_dsb_uops"
};

static const char * const IDQ_2[4] = {
       "idq.ms_mite_uops",
       "idq.ms_uops",
       "idq.all_dsb_cycles_any_uops",
       "idq.all_dsb_cycles_4_uops"
};

static const char * const IDQ_3[3] = {
       "idq.all_mite_cycles_any_uops",
       "idq.all_mite_cycles_4_uops",
       "idq.mite_all_uops"
};

static const char * const ICACHE[1] = {
       "icache.hit"
};

static const char * const ITLB_MISSES_1[4] = {
       "itlb_misses.miss_causes_a_walk",
       "itlb_misses.walk_completed_4k",
       "itlb_misses.walk_completed_2m_4m",
       "itlb_misses.walk_completed"
};

static const char * const ITLB_MISSES_2[3] = {
       "itlb_misses.stlb_hit_4k",
       "itlb_misses.stlb_hit_2m",
       "itlb_misses.stlb_hit"
};

static const char * const ILD_STALL[2] = {
       "ild_stall.lcp",
       "ild_stall.iq_full"
};

static const char * const BR_INST_EXEC_1[4] = {
       "br_inst_exec.cond",
       "br_inst_exec.direct_jmp",
       "br_inst_exec.indirect_jmp_non_call_ret",
       "br_inst_exec.return_near"
};

static const char * const BR_INST_EXEC_2[4] = {
       "br_inst_exec.direct_near_call",
       "br_inst_exec.indirect_near_call",
       "br_inst_exec.nontaken",
       "br_inst_exec.taken"
};

static const char * const BR_INST_EXEC_3[1] = {
       "br_inst_exec.all_branches"
};

static const char * const BR_MISP_EXEC_1[4] = {
        "br_misp_exec.cond",
	"br_misp_exec.indirect_jmp_non_call_ret",
	"br_misp_exec.return_near",
	"br_misp_exec.direct_near_call"
};

static const char * const BR_MISP_EXEC_2[4] = {
        "br_misp_exec.indirect_near_call",
	"br_misp_exec.nontaken",
	"br_misp_exec.taken",
	"br_misp_exec.all_branches"
};

static const char * const IDQ_UOPS_NOT_DELIVERED[1] = {
        "idq_uops_not_delivered.core"
};

static const char * const UOPS_EXECUTED_PORT_1[4] = {
        "uops_executed_port.port_0",
	"uops_executed_port.port_1",
	"uops_executed_port.port_2",
	"uops_executed_port.port_3"
};

static const char * const UOPS_EXECUTED_PORT_2[4] = {
        "uops_executed_port.port_4",
	"uops_executed_port.port_5",
	"uops_executed_port.port_6",
	"uops_executed_port.port_7"
};

static const char * const RESOURCE_STALLS[4] = {
        "resource_stalls.any",
	"resource_stalls.rs",
	"resource_stalls.sb",
	"resource_stalls.rob"
};

static const char * const CYCLE_ACTIVITY_1[4] = {
        "cycle_activity.cycles_l2_pending",
	"cycle_activity.cycles_ldm_pending",
	"cycle_activity.stalls_l2_pending",
	"cycle_activity.stalls_l1d_pending"
};

static const char * const CYCLE_ACTIVITY_2[1] = {
        "cycle_activity.cycles_l1d_pending"
};

static const char * const LSD[1] = {
        "lsd.uops"
};

static const char * const ITLB[1] = {
        "itlb.itlb_flush"
};

static const char * const OFFCORE_REQUESTS[4] = {
        "offcore_requests.demand_data_rd",
	"offcore_requests.demand_core_rd",
	"offcore_requests.demand_rfo",
	"offcore_requests.all_data_rd"
};

static const char * const UOPS_EXECUTED[1] = {
        "uops_executed.core"
};

static const char * const PAGE_WALKER_LOADS_1[4] = {
        "page_walker_loads.dtlb_l1",
	"page_walker_loads.itlb_l1",
	"page_walker_loads.dtlb_l2",
	"page_walker_loads.itlb_l2"
};

static const char * const PAGE_WALKER_LOADS_2[4] = {
        "page_walker_loads.dtlb_l3",
	"page_walker_loads.itlb_l3",
	"page_walker_loads.dtlb_memory",
	"page_walker_loads.itlb_memory"
};

static const char * const TLB_FLUSH[2] = {
        "tlb_flush.dtlb_thread",
	"tlb_flush.stlb_any"
};

static const char * const INST_RETIRED[2] = {
        "inst_retired.any",
	"inst_retired.prec_dist"
};

static const char * const OTHER_ASSIST[3] = {
        "other_assist.avx_to_sse",
	"other_assist.sse_to_avx",
	"other_assist.any_wb_assist"
};

static const char * const UOPS_RETIRED[2] = {
        "uops_retired.all",
	"uops_retired.retired_slots"
};

static const char * const MACHINE_CLEARS[3] = {
        "machine_clears.memory_ordering",
	"machine_clears.smc",
	"machine_clears.maskmov"
};

static const char * const BR_INST_RETIRED_1[4] = {
        "br_inst_retired.all_branches",
	"br_inst_retired.conditional",
	"br_inst_retired.near_call",
	"br_inst_retired.all_branches_pebs"
};

static const char * const BR_INST_RETIRED_2[4] = {
        "br_instr_retired.near_return",
	"br_instr_retired.not_taken",
	"br_instr_retired.near_taken",
	"br_instr_retired.far_branch"
};

static const char * const BR_MISP_RETIRED[4] = {
        "br_misp_retired.all_branches",
	"br_misp_retired.conditional",
	"br_misp_retired.all_branches_pebs",
	"br_misp_retired.near_taken"
};

static const char * const HLE_RETIRED_1[4] = {
        "hle_retired.start",
	"hle_retired.commit",
	"hle_retired.aborted",
	"hle_retired.aborted_mem"
};

static const char * const HLE_RETIRED_2[4] = {
         "hle_retired.aborted_timer",
	 "hle_retired.aborted_unfriendly",
	 "hle_retired.aborted_memtype",
	 "hle_retired.aborted_events"
};

static const char * const RTM_RETIRED_1[4] = {
        "rtm_retired.start",
	"rtm_retired.commit",
	"rtm_retired.aborted",
	"rtm_retired.aborted_mem"
};

static const char * const RTM_RETIRED_2[4] = {
         "rtm_retired.aborted_timer",
	 "rtm_retired.aborted_unfriendly",
	 "rtm_retired.aborted_memtype",
	 "rtm_retired.aborted_events"
};

static const char * const FP_ASSIST_1[4] = {
         "fp_assist.x87_output",
	 "fp_assist.x87_input",
	 "fp_assist.simd_output",
	 "fp_assist.simd_input"
};

static const char * const FP_ASSIST_2[1] = {
         "fp_assist.any"
};

static const char * const ROB_MISC_EVENTS[1] = {
         "rob_misc_events.lbr_inserts"
};

static const char * const MEM_UOPS_RETIRED[4] = {
         "mem_uops_retired.stlb_miss_load",
	 "mem_uops_retired.stlb_miss_stores",
	 "mem_uops_retired.lock_loads",
	 "mem_uops_retired.split_loads"
};

static const char * const MEM_UOPS_RETIRED_2[3] = {
         "mem_uops_retired.split_stores",
	 "mem_uops_retired.all_loads",
	 "mem_uops_retired.all_stores"
};

static const char * const MEM_LOAD_UOPS_RETIRED_1[4] = {
         "mem_load_uops_retired.l1_hit",
	 "mem_load_uops_retired.l2_hit",
	 "mem_load_uops_retired.l3_hit",
	 "mem_load_uops_retired.l1_miss"
};

static const char * const MEM_LOAD_UOPS_RETIRED_2[3] = {
         "mem_load_uops_retired.l2_miss",
	 "mem_load_uops_retired.l3_miss",
	 "mem_load_uops_retired.hit_lfb"
};

static const char * const MEM_LOAD_UOPS_L3_HIT_RETIRED[4] = {
         "mem_load_uops_l3_hit_retired.xsnp_miss",
	 "mem_load_uops_l3_hit_retired.xsnp_hit",
	 "mem_load_uops_l3_hit_retired.xsnp_hitm",
	 "mem_load_uops_l3_hit_retired.xsnp_none"
};


#endif /*__GMS_PMC_EVENTS_H__*/
