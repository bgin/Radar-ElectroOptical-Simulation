
#ifndef __GMS_PMC_SCHEDULES_H__
#define __GMS_PMC_SCHEDULES_H_

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
  "rs_events.empty",
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
  "dtlb_load_misses.stlb_hit"
};

static const int ALL_EVENTS_COUNT = sizeof(ALL_EVENTS)/sizeof(*ALL_EVENTS);


#endif /*__GMS_PMC_SCHEDULES_H__*/
