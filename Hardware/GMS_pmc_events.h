
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



 
  
 


#endif /*__GMS_PMC_EVENTS_H__*/
