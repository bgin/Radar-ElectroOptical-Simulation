
#ifndef __GMS_PMC_EVENTS_H__
#define __GMS_PMC_EVENTS_H_

// PMC sampling schedules

// By 4-element batch

static const char * const EVENTS_BATCH1[4] = {
  "uops_issued.any",
  "uops_issued.any<1",
  "uops_issued.any>=1",
  "uops_issued.any>=2"
};

static const char * const EVENTS_BATCH2[4] = {
  "uops_issued.any>=3",
  "uops_issued.any>=4",
  "uops_issued.any>=5",
  "uops_issued.any>=6"
};

static const char * const EVENTS_BATCH3[4] = {
  "uops_executed_port.port_0",
  "uops_executed_port.port_1",
  "uops_executed_port.port_2",
  "uops_executed_port.port_3"
};

static const char * const EVENTS_BATCH4[4] = {
  "uops_executed_port.port_4",
  "uops_executed_port.port_5",
  "uops_executed_port.port_6",
  "uops_executed_port.port_7"
};

static const char * const EVENTS_BATCH5[4] = {
  "resource_stalls.any",
  "resource_stalls.rs",
  "resource_stalls.sb",
  "resource_stalls.rob"
};

static const char * const EVENTS_BATCH6[4] = {
  "uops_retired.all",
  "uops_retired.all<1",
  "uops_retired.all>=1",
  "uops_retired.all>=2"
};

static const char * const EVENTS_BATCH7[4] = {
  "uops_retired.all>=3",
  "uops_retired.all>=4",
  "uops_retired.all>=5",
  "uops_retired.all>=6"
};

static const char * const EVENTS_BATCH8[4] = {
  "inst_retired.any_p",
  "inst_retired.any_p<1",
  "inst_retired.any_p>=1",
  "inst_retired.any_p>=2"
};

static const char * const EVENTS_BATCH9[4] = {
  "inst_retired.any_p>=3",
  "inst_retired.any_p>=4",
  "inst_retired.any_p>=5",
  "inst_retired.any_p>=6"
};

static const char * const EVENTS_BATCH10[4] = {
  "idq_uops_not_delivered.core",
  "idq_uops_not_delivered.core<1",
  "idq_uops_not_delivered.core>=1",
  "idq_uops_not_delivered.core>=2"
};

static const char * const EVENTS_BATCH11[4] = {
  "idq_uops_not_delivered.core>=3",
  "idq_uops_not_delivered.core>=4",
  "rs_events.empty_cycles",
  "idq.empty"
};

static const char * const EVENTS_BATCH12[4] = {
  "idq.mite_uops",
  "idq.dsb_uops",
  "idq.ms_dsb_uops",
  "idq.ms_mite_uops"
};

static const char * const  EVENTS_BATCH13[4] = {
  "idq.mite_all_uops",
  "idq.mite_all_uops<1",
  "idq.mite_all_uops>=1",
  "idq.mite_all_uops>=2"
};

static const char * const EVENTS_BATCH14[4] = {
   "idq.mite_all_uops>=3",
   "idq.mite_all_uops>=4",
   "move_elimination.int_not_eliminated",
   "move_elimination.simd_not_eliminated"
};

static const char * const EVENTS_BATCH15[4] = {
   "lsd.uops",
   "lsd.uops<1",
   "lsd.uops>=1",
   "lsd.uops>=2" 
};

static const char * const EVENTS_BATCH16[4] = {
    "lsd.uops>=3",
    "lsd.uops>=4",
    "ild_stall.lcp",
    "ild_stall.iq_full"
};

static const char * const EVENTS_BATCH17[4] = {
    "br_inst_exec.all_branches",
    "br_inst_exec.0x81",
    "br_inst_exec.0x82",
    "icache.misses"
};

static const char * const EVENTS_BATCH18[4] = {
     "br_misp_exec.all_branches",
     "br_misp_exec.0x81",
     "br_misp_exec.0x82",
     "fp_assist.any"
};

static const char * const EVENTS_BATCH19[4] = {
     "cpu_clk_unhalted.core_clk",
     "cpu_clk_unhalted.ref_xclk",
     "baclears.any",
     "idq.ms_uops"
};

static const char * const EVENTS_BATCH20[4] = {
     "ld_blocks.store_forward",
     "ld_blocks.no_sr",
     "misalign_mem_reference.loads",
     "misalign_mem_reference.stores"
};

static const char * const EVENTS_BATCH21[4] = {
     "ld_blocks_partial_address_alias",
     "dtlb_load_misses.miss_causes_a_walk",
     "dtlb_load_misses.walk_completed_4k",
     "dtlb_load_misses.walk_completed_2m_4m"
};

static const char * const EVENTS_BATCH22[4] = {
    "dtlb_load_misses.walk_completed",
    "dtlb_load_misses.stlb_hit_4k",
    "dtlb_load_misses.stlb_hit_2m",
    "dtlb_load_misses.stlb_hit"
};

static const char * const EVENTS_BATCH23[4] = {
    "dtlb_load_misses.pde_cache_miss",
    "int_misc.recovery_cycles",
    "l2_rqsts.rfo_miss",
    "l2_rqsts.rfo_hit"
};

static const char * const EVENTS_BATCH24[4] = {
    "l2_rqsts.all_rfo",
    "l2_rqsts.demand_data_rd_miss",
    "l2_rqsts.demand_data_rd_hit",
    "l2_rqsts.all_demand_data_rd"
};

static const char * const EVENTS_BATCH25[4] = {
    "l2_rqsts.code_rd_hit",
    "l2_rqsts.code_rd_miss",
    "l2_rqsts.all_demand_miss",
    "l2_rqsts.all_demand_references"
};

static const char * const EVENTS_BATCH26[4] = {
    "l2_rqsts.all_code_rd",
    "l2_rqsts.l2_pf_hit",
    "l2_rqsts.l2_pf_miss",
    "l2_rqsts.all_pf"
};

static const char * const EVENTS_BATCH27[4] = {
     "l2_rqsts.miss",
     "l2_rqsts.references",
     "l2_demand_rqsts.wb_hit",
     "llc.miss"
};

static const char * const EVENTS_BATCH28[4] = {
     "llc.hit",
     "cpu_clk_unhalted.core_clk",
     "cpu_clk_unhalted.ref_xclk",
     "l1_pend_miss.pending"
};

static const char * const EVENTS_BATCH29[4] = {
      "dtlb_store_misses.miss_causes_a_walk",
      "dtlb_store_misses.walk_completed_4k",
      "dtlb_store_misses.walk_completed_2m_4m",
      "dtlb_store_misses.walk_completed"
};

static const char * const EVENTS_BATCH30[4] = {
      "dtlb_store_misses.walk_duration",
      "dtlb_store_misses.stlb_hit_4k",
      "dtlb_store_misses.stlb_hit_2m",
      "dtlb_store_misses.stlb_hit"
};

static const char * const EVENTS_BATCH31[4] = {
      "dtlb_store_misses.pde_cache_miss",
      "load_hit_pre.sw_pf",
      "load_hit_pre.hw_pf",
      "l1d.replacement"
};

static const char * const EVENTS_BATCH32[4] = {
      "tx_mem.abort_conflict",
      "tx_mem.abort_capacity_write",
      "tx_mem.abort_hle_store_to_elided_lock",
      "tx_mem.abort_hle_elision_buffer_not_empty"
};

static const char * const EVENTS_BATCH33[4] = {
      "tx_mem.abort_hle_elision_buffer_mismatch",
      "tx_mem.abort_hle_elision_buffer_unsupported_alignment",
      "tx_mem.hle_elision_buffer_full",
      "move_elimination.int_eliminated"
};

static const char * const EVENTS_BATCH34[4] = {
      "move_elimination.simd_eliminated",
      "cpl_cycles.ring0",
      "cpl_cycles.ring123",
      "tx_exec.misc1"
};

static const char * const EVENTS_BATCH35[4] = {
      "tx_exec.misc2",
      "tx_exec.misc3",
      "tx_exec.misc4",
      "tx_exec.misc5"
};

static const char * const EVENTS_BATCH36[4] = {
      "offcore_requests_outstanding.demand_data_rd",
      "offcore_requests_outstanding.demand_code_rd",
      "offcore_requests_outstanding.demand_rfo",
      "offcore_requests_outstanding.all_data_rd"
};

static const char * const EVENTS_BATCH37[4] = {
      "lock_cycles.split_lock_uc_lock_duration",
      "lock_cycles.cache_lock_duration",
      
};

static const char * const EVENTS_BATCH38[4] = {
      "idq.ms_uops",
      "idq.all_dsb_cycles_any_uops",
      "idq.all_dsb_cycles_any_uops<1",
      "idq.all_dsb_cycles_any_uops>=1"
};

static const char * const EVENTS_BATCH39[4] = {
      "idq.all_dsb_cycles_any_uops>=2",
      "idq.all_dsb_cycles_any_uops>=3",
      "idq.all_dsb_cycles_any_uops>=4",
      "idq.all_dsb_cycles_any_uops>=5"
};

static const char * const EVENTS_BATCH40[4] = {
      "idq.all_dsb_cycles_any_uops>=6",
      "idq.all_dsb_cycles_4_uops",
      "idq.all_dsb_cycles_4_uops<1",
      "idq.all_dsb_cycles_4_uops>=1"
}

static const char * const EVENTS_BATCH41[4] = {
      "idq.all_mite_cycles_any_uops",
      "idq.all_mite_cycles_any_uops<1",
      "idq.all_mite_cycles_any_uops>=1",
      "idq.all_mite_cycles_any_uops>=2"
};

static const char * const EVENTS_BATCH42[4] = {
      "idq.all_mite_cycles_any_uops>=3",
      "idq.all_mite_cycles_any_uops>=4",
      "idq.all_mite_cycles_any_uops>=5",
      "idq.all_mite_cycles_any_uops>=6"
};

static const char * const EVENTS_BATCH43[4] = {
      "idq.all_mite_cycles_4_uops"
      "idq.all_mite_cycles_4_uops<1",
      "idq.all_mite_cycles_4_uops>=1",
      "idq.all_mite_cycles_4_uops>=2",
      "idq.all_mite_cycles_4_uops>=3"
     
};

static const char * const EVENTS_BATCH44[4] = {
       "idq.all_mite_cycles_4_uops>=4",
       "idq.all_mite_cycles_4_uops>=5",
       "idq.all_mite_cycles_4_uops>=6",
       "idq.mite_all_uops"
};

static const char * const EVENTS_BATCH45[4] = {
       "idq.mite_all_uops<1",
       "idq.mite_all_uops>=1",
       "idq.mite_all_uops>=2",
       "idq.mite_all_uops>=3"
};

static const char * const EVENTS_BATCH46[4] = {
       "idq.mite_all_uops>=4",
       "idq.mite_all_uops>=5",
       "idq.mite_all_uops>=6",
       "icache.misses"
};

static const char * const EVENTS_BATCH47[4] = {
       "itlb_misses.miss_causes_a_walk",
       "itlb_misses.miss_causes_a_walk<1",
       "itlb_misses.miss_causes_a_walk>=1",
       "itlb_misses.miss_causes_a_walk>=2"
};

static const char * const EVENTS_BATCH48[4] = {
       "itlb_misses.miss_causes_a_walk>=3",
       "itlb_misses.miss_causes_a_walk>=4",
       "itlb_misses.miss_causes_a_walk>=5",
       "itlb_misses.miss_causes_a_walk>=6"
};

static const char * const EVENTS_BATCH49[4] = {
      "itlb_misses.walk_completed_4k",
      "itlb_misses.walk_completed_4k<1",
      "itlb_misses.walk_completed_4k>=1",
      "itlb_misses.walk_completed_4k>=2"
};

static const char * const EVENTS_BATCH50[4] = {
      "itlb_misses.walk_completed_4k>=3",
      "itlb_misses.walk_completed_4k>=4",
      "itlb_misses.walk_completed_4k>=5",
      "itlb_misses.walk_completed_4k>=6"
};

static const char * const EVENTS_BATCH51[4] = {
      "itlb_misses.walk_completed_2m_4m",
      "itlb_misses.walk_completed_2m_4m<1",
      "itlb_misses.walk_completed_2m_4m>=1",
      "itlb_misses.walk_completed_2m_4m>=2"
};

static const char * const EVENTS_BATCH52[4] = {
      "itlb_misses.walk_completed_2m_4m>=3",
      "itlb_misses.walk_completed_2m_4m>=4",
      "itlb_misses.walk_completed_2m_4m>=5",
      "itlb_misses.walk_completed_2m_4m>=6"
};

static const char * const EVENTS_BATCH53[4] = {
      "itlb_misses.walk_completed",
      "itlb_misses.walk_completed<1",
      "itlb_misses.walk_completed>=1",
      "itlb_misses.walk_completed>=2"
};

static const char * const EVENTS_BATCH54[4] = {
      "itlb_misses.walk_completed>=3",
      "itlb_misses.walk_completed>=4",
      "itlb_misses.walk_completed>=5",
      "itlb_misses.walk_completed>=6"
};

static const char * const EVENTS_BATCH55[4] = {
      "itlb_misses.walk_duration",
      "itlb_misses.walk_duration<1",
      "itlb_misses.walk_duration>=1",
      "itlb_misses.walk_duration>=2"
};

static const char * const EVENTS_BATCH56[4] = {
      "itlb_misses.walk_duration>=3",
      "itlb_misses.walk_duration>=4",
      "itlb_misses.walk_duration>=5",
      "itlb_misses.walk_duration>=6"
};

static const char * const EVENTS_BATCH57[4] = {
       "itlb_misses.stlb_hit_4k",
       "itlb_misses.stlb_hit_4k<1",
       "itlb_misses.stlb_hit_4k>=1",
       "itlb_misses.stlb_hit_4k>=2"
};

static const char * const EVENTS_BATCH58[4] = {
       "itlb_misses.stlb_hit_4k>=3",
       "itlb_misses.stlb_hit_4k>=4",
       "itlb_misses.stlb_hit_4k>=5",
       "itlb_misses.stlb_hit_4k>=6"
};

static const char * const EVENTS_BATCH59[4] = {
       "itlb_misses.stlb_hit_2m",
       "itlb_misses.stlb_hit_2m<1",
       "itlb_misses.stlb_hit_2m>=1",
       "itlb_misses.stlb_hit_2m>=2"
};

static const char * const EVENTS_BATCH60[4] = {
       "itlb_misses.stlb_hit_2m>=3",
       "itlb_misses.stlb_hit_2m>=4",
       "itlb_misses.stlb_hit_2m>=5",
       "itlb_misses.stlb_hit_2m>=6"
};

static const char * const EVENTS_BATCH61[4] = {
       "itlb_misses.stlb_hit",
       "itlb_misses.stlb_hit<1",
       "itlb_misses.stlb_hit>=1",
       "itlb_misses.stlb_hit>=2"
};

static const char * const EVENTS_BATCH62[4] = {
       "itlb_misses.stlb_hit>=3",
       "itlb_misses.stlb_hit>=4",
       "itlb_misses.stlb_hit>=5",
       "itlb_misses.stlb_hit>=6"
};

static const char * const EVENTS_BATCH63[4] = {
       "branch_inst_exec.cond",
       "branch_inst_exec.direct_jump",
       "branch_inst_exec.indirect_jump_non_call_ret",
       "branch_inst_exec.return_near"
};

static const char * const EVENTS_BATCH64[4] = {
       "branch_inst_exec.direct_near_call",
       "branch_inst_exec.indirect_near_call",
       "branch_inst_exec.nontaken",
       "branch_inst_exec.all_branches"
};




static const char * const ALL_EVENTS[] = {


 
  // By 4-elements batch
  "uops_issued.any",
  "uops_issued.any<1",
  "uops_issued.any>=1",
  "uops_issued.any>=2",
  "uops_issued.any>=3",
  "uops_issued.any>=4",
  "uops_issued.any>=5",
  "uops_issued.any>=6",
  "uops_executed_port.port_0",
  "uops_executed_port.port_1",
  "uops_executed_port.port_2",
  "uops_executed_port.port_3",
   "uops_executed_port.port_4",
  "uops_executed_port.port_5",
  "uops_executed_port.port_6",
  "uops_executed_port.port_7",
  "resource_stalls.any",
  "resource_stalls.rs",
  "resource_stalls.sb",
  "resource_stalls.rob",
  "uops_retired.all",
  "uops_retired.all<1",
  "uops_retired.all>=1",
  "uops_retired.all>=2",
  "uops_retired.all>=3",
  "uops_retired.all>=4",
  "uops_retired.all>=5",
  "uops_retired.all>=6",
  "inst_retired.any_p",
  "inst_retired.any_p<1",
  "inst_retired.any_p>=1",
  "inst_retired.any_p>=2",
  "inst_retired.any_p>=3",
  "inst_retired.any_p>=4",
  "inst_retired.any_p>=5",
  "inst_retired.any_p>=6",
  "idq_uops_not_delivered.core",
  "idq_uops_not_delivered.core<1",
  "idq_uops_not_delivered.core>=1",
  "idq_uops_not_delivered.core>=2",
  "idq_uops_not_delivered.core>=3",
  "idq_uops_not_delivered.core>=4",
  "rs_events.empty",
  "idq.empty",
   "idq.mite_uops",
  "idq.dsb_uops",
  "idq.ms_dsb_uops",
  "idq.ms_mite_uops",
  "idq.mite_all_uops",
  "idq.mite_all_uops<1",
  "idq.mite_all_uops>=1",
  "idq.mite_all_uops>=2",
  "idq.mite_all_uops>=3",
  "idq.mite_all_uops>=4",
  "move_elimination.int_not_eliminated",
  "move_elimination.simd_not_eliminated",
  "lsd.uops",
  "lsd.uops<1",
  "lsd.uops>=1",
  "lsd.uops>=2",
  "lsd.uops>=3",
  "lsd.uops>=4",
  "ild_stall.lcp",
  "ild_stall.iq_full",
  "br_inst_exec.all_branches",
  "br_inst_exec.0x81",
  "br_inst_exec.0x82",
  "icache.misses",
  "br_misp_exec.all_branches",
  "br_misp_exec.0x81",
  "br_misp_exec.0x82",
  "fp_assist.any",
  "cpu_clk_unhalted.core_clk",
  "cpu_clk_unhalted.ref_xclk",
  "baclears.any",
  "idq.ms_uops",
  "ld_blocks.store_forward",
  "ld_blocks.no_sr",
  "misalign_mem_reference.loads",
  "misalign_mem_reference.stores",
  "ld_blocks_partial_address_alias",
  "dtlb_load_misses.miss_causes_a_walk",
  "dtlb_load_misses.walk_completed_4k",
  "dtlb_load_misses.walk_completed_2m_4m",
  "dtlb_load_misses.walk_completed",
  "dtlb_load_misses.stlb_hit_4k",
  "dtlb_load_misses.stlb_hit_2m",
  "dtlb_load_misses.stlb_hit",
  "dtlb_load_misses.pde_cache_miss",
  "int_misc.recovery_cycles",
  "l2_rqsts.rfo_miss",
  "l2_rqsts.rfo_hit",
  "l2_rqsts.all_rfo",
  "l2_rqsts.demand_data_rd_miss",
  "l2_rqsts.demand_data_rd_hit",
  "l2_rqsts.all_demand_data_rd",
  "l2_rqsts.code_rd_hit",
  "l2_rqsts.code_rd_miss",
  "l2_rqsts.all_demand_miss",
  "l2_rqsts.all_demand_references",
  "l2_rqsts.all_code_rd",
  "l2_rqsts.l2_pf_hit",
  "l2_rqsts.l2_pf_miss",
  "l2_rqsts.all_pf",
  "l2_rqsts.miss",
  "l2_rqsts.references",
  "l2_demand_rqsts.wb_hit",
  "llc.miss",
  "llc.hit",
  "cpu_clk_unhalted.core_clk",
  "cpu_clk_unhalted.ref_xclk",
  "l1_pend_miss.pending",
  "dtlb_store_misses.miss_causes_a_walk",
  "dtlb_store_misses.walk_completed_4k",
  "dtlb_store_misses.walk_completed_2m_4m",
  "dtlb_store_misses.walk_completed",
  "dtlb_store_misses.walk_duration",
  "dtlb_store_misses.stlb_hit_4k",
  "dtlb_store_misses.stlb_hit_2m",
  "dtlb_store_misses.stlb_hit",
  "dtlb_store_misses.pde_cache_miss",
  "load_hit_pre.sw_pf",
  "load_hit_pre.hw_pf",
  "l1d.replacement",
  
};

static const int ALL_EVENTS_COUNT = sizeof(ALL_EVENTS)/sizeof(*ALL_EVENTS);


#endif /*__GMS_PMC_EVENTS_H__*/
