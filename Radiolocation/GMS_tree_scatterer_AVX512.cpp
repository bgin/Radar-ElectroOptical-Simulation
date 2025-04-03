

#include <omp.h>
#if defined __ICC || defined __INTEL_COMPILER
#include <svrng.h>
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#include <random>
#include <ctime>
#include <cstdlib>
#endif
//
#include "GMS_tree_scatterer_AVX512.h"
//
#include "GMS_malloc.h"
#include "GMS_tree_scatterer_common.h"
#include "GMS_indices.h"
#include "GMS_common.h"

#if !defined(TREE_SCATTERER_AVX512_HW_PMC_PROLOGE_BODY)
      #define TREE_SCATTERER_AVX512_HW_PMC_PROLOGE_BODY             \
        do {                                                        \
               uint64_t prog_counters_start[4] = {};                \
               uint64_t prog_counters_end[4]   = {};                \
               uint64_t tsc_start,tsc_stop;                         \
               uint64_t act_cyc_start,act_cyc_stop;                 \
               uint64_t ref_cyc_start,ref_cyc_stop;                 \
               uint64_t volatile dummy1,dummy2,dummy3,dummy4;       \
               int32_t core_counter_width;                          \
               double utilization,nom_ghz,avg_ghz;                  \
               dummy1 = rdtsc();                                    \
               dummy2 = rdtsc();                                    \
               dummy3 = rdpmc(0);                                   \
               core_counter_width = get_core_counter_width();       \
               for(int32_t i = 0; i != 4; ++i) {                    \
                   prog_counters_start[i] = rdpmc(i);               \
               }                                                    \
               act_cyc_start = rdpmc_actual_cycles();               \
               ref_cyc_start = rdpmc_reference_cycles();            \
               tsc_start = rdtsc();                                
	} while(0)
#endif

#if !defined(TREE_SCATTERER_AVX512_HW_PMC_EPILOGE_BODY)
    #define TREE_SCATTERER_AVX512_HW_PMC_EPILOGE_BODY                                                                                                                         \
       do {                                                                                                                                                                   \
                for(int32_t i = 0; i != 4; ++i) {                                                                                                                             \
                     prog_counters_stop[i] = rdpmc(i);                                                                                                                        \                                           }                                                                                                                                                             \
                act_cyc_stop = rdpmc_actual_cycles();                                                                                                                         \
                ref_cyc_stop = rdpmc_reference_cycles();                                                                                                                      \
                tsc_stop = rdtsc();                                                                                                                                           \
	        dummy4 = get_core_counter_width();                                                                                                                            \
	        nom_ghz = get_TSC_frequency()/1.0e9;                                                                                                                          \
	        utilization = (double)(ref_cyc_end-ref_cyc_start)/(double)(tsc_end-tsc_start-rdtscp_latency);                                                                 \
	        avg_ghz = (double)(act_cyc_end-act_cyc_start)/(double)(tsc_end-tsc_start-rdtscp_latency)*nom_ghz;                                                             \
	        syslog(LOG_INFO,"%-10s:\n", __PRETTY_FUNCTION__);                                                                                                             \
	        syslog(LOG_INFO, "*************** Hardware Counters -- Dump Begin **************");                                                                           \
	        syslog(LOG_INFO,"Core utilization                      : %f\n",utilization );                                                                                 \
	        syslog(LOG_INFO,"Core average frequency                : %f\n",avg_ghz);                                                                                      \
	        syslog(LOG_INFO,"Reference cycles                      : %20lld\n",ref_cyc_end-ref_cyc_start);                                                                \
	        syslog(LOG_INFO,"Actual cycles                         : %20lld\n",act_cyc_stop-act_cyc_start);                                                               \
	        syslog(LOG_INFO,"L2 cycles  at entry                   : %20lld\n",cycles_lvl2);                                                                              \
	        syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event1                    , corrected_pmc_delta(prog_counters_stop[0],prog_counters_start[0],core_counter_width);      \
	        syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event2                    , corrected_pmc_delta(prog_counters_stop[1],prog_counters_start[1],core_counter_width);      \
	        syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event3                    , corrected_pmc_delta(prog_counters_stop[2],prog_counters_start[2],core_counter_width);      \
	        syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event4                    , corrected_pmc_delta(prog_counters_stop[3],prog_counters_start[3],core_counter_width);      \    
	        syslog(LOG_INFO, "*************** Hardware Counters -- Dump End   **************");
       } while(0)
#endif
      
gms::math::
TreeScattererAVX512::TreeScattererAVX512() {

     m_tsc.nleaves             = -1;
     m_tsc.nbranches           = -1;
     m_tsc.nsteps              = -1;
     m_tsc.ordinal             = -1;
     m_tsc.trunk_param_npts    = -1;
     m_tsc.leaves_param_npts   = -1;
     m_tsc.branches_param_npts = -1;
     m_tsc.tree_height         = -1.0f;
     m_tsc.trunk_height        = -1.0f;
     m_tsc.trunk_radius        = -1.0f;
     m_tsc.crown_height        = -1.0f;
     m_tsc.crown_area          = -1.0f;
     m_tsc.trunk_area          = -1.0f;
     m_tsc.tree_area           = -1.0f;
     m_tsc.tree_lat            = -1.0f;
     m_tsc.tree_lon            = -1.0f;
     m_tsc_elevation           = -1.0f;
     m_tsc.leaves_moist        = NULL;
     m_tsc.branches_moist      = NULL;
     m_tsc.trunk_xparam        = NULL;
     m_tsc.trunk_yparam        = NULL;
     m_tsc.trunk_zparam        = NULL;
     m_tsc.leaves_thick        = NULL;
     m_tsc.leaves_dens         = NULL;
     m_tsc.leaves_xparam       = NULL;
     m_tsc.leaves_yparam       = NULL;
     m_tsc.branches_thick      = NULL;
     m_tsc.branches_dens       = NULL;
     m_tsc.branches_xparam     = NULL;
     m_tsc.branches_yparam     = NULL;
     m_tsc.branches_zparam     = NULL;
     m_tsh.tree_dtheta         = 0.0f;
     m_tsh.tree_dphi           = 0.0f;
     m_tsh.tree_rcs            = 0.0f;
     m_tsh.crown_rcs           = 0.0f;
     m_tsh.trunk_rcs           = 0.0f;
     m_tsh.leaves_rcs          = NULL;
     m_tsh.leaves_reflect      = NULL;
     m_tsh.branches_rcs        = NULL;
     m_tsh.branches_reflect    = NULL;
     m_lp.l4x4phm              = NULL;
     m_lp.l2x2mp               = NULL;
     m_lp.l2x2mn               = NULL;
     m_lp.stokes4x4m           = NULL;
     m_lp.scat2x2m             = NULL;

}


#if !defined(GMS_TREE_SCATTERER_AVX512_COLD_ALLOC_CTOR)
    #define GMS_TREE_SCATTERER_AVX512_COLD_ALLOC_CTOR                                                                   \
      m_tsc.leaves_moist        = gms_eimalloca4(static_cast<size_t>(m_tsc.nleaves),64);                                \
      m_tsc.branches_moist      = gms_eimalloca4(static_cast<size_t>(m_tsc.nbranches),64);                              \          
      m_tsc.trunk_xparam        = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.trunk_param_npts),64);                 \
      m_tsc.trunk_yparam        = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.trunk_param_npts),64);                 \
      m_tsc.trunk_zparam        = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.trunk_param_npts),64);                 \
      m_tsc.leaves_thick        = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nleaves),64);                          \
      m_tsc.leaves_dens         = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nleaves),64);                          \
      m_tsc.leaves_xparam       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nleaves*m_tsc.leaves_param_npts),64);  \
      m_tsc.leaves_yparam       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nleaves*m_tsc.leaves_param_npts),64);  \
      m_tsc.branches_thick      = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches),64);                        \
      m_tsc.branches_dens       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches),64);                        \
      m_tsc.branches_xparam     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches*m_tsc.branches_param_npts),64);   \
      m_tsc.branches_yparam     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches*m_tsc.branches_param_npts),64);   \
      m_tsc.branches_zparam     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches*m_tsc.branches_param_npts),64); 
#endif

#if !defined(GMS_TREE_SCATTERER_AVX512_HOT_ALLOC_CTOR)
    #define GMS_TREE_SCATTERER_AVX512_HOT_ALLOC_CTOR                                                                      \
      m_tsh.leaves_rcs          = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_reflect      = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \ 
      m_tsh.branches_rcs        = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_reflect    = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         
    
#endif

#if !defined(GMS_TREE_SCATTERER_LEAVES_PHASE_CTOR_BODY)
    #define GMS_TREE_SCATTERER_LEAVES_PHASE_CTOR_BODY                                                 \
   
     m_lp.l4x4phm          = gms_efmalloca(static_cast<size_t>(4*4*4*m_tsc.nleaves),64);                     \
     m_lp.l2x2mp           = gms_cmplxr4_emalloca(static_cast<size_t>(2*2*m_tsc.nleaves),64);                  \
     m_lp.l2x2mn           = gms_cmplxr4_emalloca(static_cast<size_t>(2*2*m_tsc.nleaves),64);                  \
     m_lp.stokes4x4m       = gms_efmalloca(static_cast<size_t>(4*4*m_tsc.nleaves),64);                        \
     m_lp.scat2x2m         = gms_cmplxr4_efmalloca(static_cast<size_t>(2*2*m_tsc.nleaves),64);
#endif

gms::math::
TreeScattererAVX512::TreeScattererAVX512( const int32_t nleaves,
                                          const int32_t nbranches,
			                  const int32_t nsteps,
			                  const int32_t ordinal,
			                  const int32_t trunk_param_npts,
			                  const int32_t leaves_param_npts,
			                  const int32_t branches_param_npts,
			                  const float   tree_height,
			                  const float   trunk_height,
			                  const float   trunk_radius,
			                  const float   crown_height,
                                          const float   tree_lat,
			                  const float   tree_lon,
			                  const float   tree_elev) {
      using namespace gms::common;
      m_tsc.nleaves             = nleaves;
      m_tsc.nbranches           = nbranches;
      m_tsc.nsteps              = nsteps;
      m_tsc.ordinal             = rdinal;
      m_tsc.trunk_param_npts    = trunk_param_npts;
      m_tsc.leaves_param_npts   = leaves_param_npts;
      m_tsc.branches_param_npts = branches_param_npts;
      m_tsc.tree_height         = tree_height;
      m_tsc.trunk_height        = trunk_height;
      m_tsc.trunk_radius        = trunk_radius;
      m_tsc.crown_height        = crown_height;
      m_tsc.crown_area          = 0.0f;
      m_tsc.trunk_area          = 0.0f;
      m_tsc.tree_area           = 0.0f;
      m_tsc.tree_lat            = tree_lat;
      m_tsc.tree_lon            = tree_lon;
      m_tsc.tree_elevation      = tree_elev;
      GMS_TREE_SCATTERER_AVX512_COLD_ALLOC_CTOR
      m_tsh.tree_dtheta         = 0.0f;
      m_tsh.tree_dphi           = 0.0f;
      m_tsh.tree_rcs            = 0.0f;
      m_tsh.crown_rcs           = 0.0f;
      m_tsh.trunk_rcs           = 0.0f;
      GMS_TREE_SCATTERER_AVX512_HOT_ALLOC_CTOR
      GMS_TREE_SCATTERER_LEAVES_PHASE_CTOR_BODY   
}


gms::math::
TreeScattererAVX512::~TreeScattererAVX512() {

      _mm_free(m_tsc.leaves_moist);
     m_tsc.leaves_moist    = NULL;
     _mm_free(m_tsc.branches_moist);
     m_tsc.branches_moist  = NULL;
     _mm_free(m_tsc.trunk_xparam);
     m_tsc.trunk_xparam    = NULL;
     _mm_free(m_tsc.trunk_yparam);
     m_tsc.trunk_yparam    = NULL;
     _mm_free(m_tsc.trunk_yparam);
     m_tsc.trunk_yparam    = NULL;
     _mm_free(m_tsc.trunk_zparam);
     m_tsc.trunk_zparam    = NULL;
     _mm_free(m_tsc.leaves_thick);
     m_tsc.leaves_thick    = NULL;
     _mm_free(m_tsc.leaves_dens);
     m_tsc.leaves_dens     = NULL;
     _mm_free(m_tsc.leaves_xparam);
     m_tsc.leaves_xparam   = NULL;
     _mm_free(m_tsc.leaves_yparam);
     m_tsc.leaves_yparam   = NULL;
     _mm_free(m_tsc.branches_thick);
     m_tsc.branches_thick  = NULL;
     _mm_free(m_tsc.branches_dens);
     m_tsc.branches_dens   = NULL;
     _mm_free(m_tsc.branches_xparam);
     m_tsc.branches_xparam = NULL;
     _mm_free(m_tsc.branches_yparam);
     m_tsc.branches_yparam = NULL;
     _mm_free(m_tsc.branches_zparam);
     m_tsc.branches_zparam = NULL;
     _mm_free(m_tsh.leaves_rcs);
     m_tsh.leaves_rcs      = NULL;
     _mm_free(m_tsh.leaves_reflect);
     m_tsh.leaves_reflect  = NULL;
     _mm_free(m_tsh.branches_rcs);
     m_tsh.branches_rcs    = NULL;
     _mm_free(m_tsh.branches_reflect);
     m_tsh.branches_reflect = NULL;
     _mm_free(m_lp.l4x4phm);
     m_lp.l4x4phm           = NULL;
     _mm_free(m_lp.l2x2mp);
     m_lp.l2x2mp            = NULL;
     _mm_free(m_lp.l2x2mn);
     m_lp.l2x2mn            = NULL;
     _mm_free(m_lp.stokes4x4m);
     m_lp.stokes4x4m        = NULL;
     _mm_free(m_lp.scat2x2m);
     m_lp.scat2x2m          = NULL;
}

void
gms::math::
TreeScattererAVX512::
SetMistnessMask() {

   
    std::random_device rd;
     // Memory first touch here!!
    for(int32_t i = 0; i != m_tsc.nleaves-7; i += 8){
        m_tsc.leaves_moist[i+0] = 0;
	m_tsc.leaves_moist[i+1] = 0;
	m_tsc.leaves_moist[i+2] = 0;
	m_tsc.leaves_moist[i+3] = 0;
	m_tsc.leaves_moist[i+4] = 0;
	m_tsc.leaves_moist[i+5] = 0;
	m_tsc.leaves_moist[i+6] = 0;
	m_tsc.leaves_moist[i+7] = 0;
    }
    std::mt19937 rgen(rd());
    std::uniform_int_distribution<> distr(0,1);
    for(int32_t i = 0; i != m_tsc.nleaves-3; i += 4) {
        m_tsc.leaves_moist[i+0] = distr(rgen);
	m_tsc.leaves_moist[i+1] = distr(rgen);
	m_tsc.leaves_moist[i+2] = distr(rgen);
	m_tsc.leaves_moist[i+3] = distr(rgen);
    }
    for(int32_t i = 0; i != m_tsc.nbranches-7; i += 8) {
        m_tsc.branches_moist[i+0] = 0;
	m_tsc.branches_moist[i+1] = 0;
	m_tsc.branches_moist[i+2] = 0;
	m_tsc.branches_moist[i+3] = 0;
	m_tsc.branches_moist[i+4] = 0;
	m_tsc.branches_moist[i+5] = 0;
	m_tsc.branches_moist[i+6] = 0;
	m_tsc.branches_moist[i+7] = 0;
        
    }
    for(int32_t i = 0; i != m_tsc.nbranches-3; i += 4) {
        m_tsc.branches_moist[i+0] = distr(rgen);
	m_tsc.branches_moist[i+1] = distr(rgen);
	m_tsc.branches_moist[i+2] = distr(rgen);
	m_tsc.branches_moist[i+3] = distr(rgen);
    }
   
}

void
gms::math::TreeScattererAVX512::
ComputeTrunkParamEq_zmm16r4(const int32_t zpoints ) {
                          

       struct _T0_ {
        AVX512Vec16 vtheta0;
	AVX512Vec16 vtheta1;
	AVX512Vec16 vtheta2;
	AVX512Vec16 vtheta3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
        AVX512Vec16 vthinc0;
	AVX512Vec16 vthinc1;
	AVX512Vec16 vthinc2;
	AVX512Vec16 vthinc3;
     } __ATTR_ALIGN__(64) t1;

     struct _T2_ {
        AVX512Vec16 vhinc0;
	AVX512Vec16 vhinc1;
	AVX512Vec16 vhinc2;
	AVX512Vec16 vhinc3;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
        AVX512Vec16 vhinit0;
	AVX512Vec16 vhinit1;
	AVX512Vec16 vhinit2;
	AVX512Vec16 vhinit3;
     } __ATTR_ALIGN__(64) t3;

     AVX512Vec16 tmp1, tmp2;
     AVX512Vec16 vrad, vNPTS;
     // Locals first-touch
     t0.vtheta0 = ZERO;
     t0.vtheta1 = ZERO;
     t0.vtheta2 = ZERO;
     t0.vtheta3 = ZERO;
     t1.vthinc0 = ZERO;
     t1.vthinc1 = ZERO;
     t1.vthinc2 = ZERO;
     t1.vthinc3 = ZERO;
     t2.vhinc0  = ZERO;
     t2.vhinc1  = ZERO;
     t2.vhinc2  = ZERO;
     t2.vhinc3  = ZERO;
     t3.vhinit0 = ZERO;
     t3.vhinit1 = ZERO;
     t3.vhinit2 = ZERO;
     t3.vhinit3 = ZERO;
     vrad       = ZERO;
     vNPTS      = ZERO;
     tmp1       = ZERO;
     tmp2       = ZERO;
     vNPTS      = AVX512Vec16{static_cast<float>(m_tsc.trunk_param_npts)};
     tmp1       = twopi/vNPTS;
     t1.vthinc0 = tmp1;
     t1.vthinc0 = t1.vthinc0*VINC0;
     t1.vthinc1 = tmp1;
     t1.vthinc1 = t1.vthinc1*VINC1;
     t1.vthinc2 = tmp1;
     t1.vthinc2 = t1.vthinc2*VINC2;
     t1.vthinc3 = tmp1;
     t1.vthinc3 = t1.vthinc3*VINC3;
     vrad       = AVX512Vec16{m_tsc.trunk_radius};
     zpoints    = zpoints+m_tsc.trunk_param_npts;
     tmp2       = AVX512Vec16{height/static_cast<float>(zpoints)};
     t2.vhinc0  = tmp2;
     t2.vhinc0  = t2.vhinc0*VINC0;
     t2.vhinc1  = tmp2;
     t2.vhinc1  = t2.vhinc1*VINC1;
     t2.vhinc2  = tmp2;
     t2.vhinc2  = t2.vhinc2*VINC2;
     t2.vhinc3  = tmp2;
     t2.vhinc3  = t2.vhinc3*VINC3;
     // First memory touch.
     gms::common::avx512vec16_init_unroll8x(&m_tsc.trunk_xparam[0],
                                        static_cast<int64_t>(m_tsc.trunk_param_npts),
					ZERO);
     gms::common::avx512vec16_init_unroll8x(&m_tsc.trunk_yparam[0],
                                        static_cast<int64_t>(m_tsc.trunk_param_npts),
					ZERO);
     gms::common::avx512vec16_init_unroll8x(&m_tsc.trunk_zparam[0],
                                        static_cast<int64_t>(m_tsc.trunk_param_npts),
					ZERO);
#if defined __GNUC__ && !defined __INTEL_COMPILER
     
   #pragma omp simd aligned(m_tsc.trunk_xparam,m_tsc.trunk_yparam,m_tsc.trunk_zparam:64)
#elif defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(m_tsc.trunk_xparam,64);
     __assume_aligned(m_tsc.trunk_yparam,64);
     __assume_aligned(m_tsc.trunk_zparam,64);
   #pragma vector always
   #pragma vector vectorlength(16)
#pragma code_align(32)
#endif
       for(int32_t i = 0; i != tsc.m_trunk_param_npts-3; i += 4) {

	   t0.vtheta = t0.vtheta0+t1.vthinc0;
	   m_tsc.trunk_xparam[i+0] = vrad*cos(t0.vtheta0);
	   m_tsc.trunk_yparam[i+0] = vrad*sin(t0.vtheta0);
	   t0.vtheta1 = t0.vtheta1+t1.vthinc1;
	   m_tsc.trunk_xparam[i+1] = vrad*cos(t0.vtheta1);
	   m_tsc.trunk_yparam[i+1] = vrad*sin(t0.vtheta1);
	   t0.vtheta2 = t0.vtheta2+t1.vthinc2;
	   m_tsc.trunk_xparam[i+2] = vrad*cos(t0.vtheta2);
	   m_tsc.trunk_yparam[i+2] = vrad*sin(t0.vtheta2);
	   t0.vtheta3 = t0.vtheta3+t1.vthinc3;
	   m_tsc.trunk_xparam[i+3] = vrad*cos(t0.vtheta3);
	   m_tsc.trunk_yparam[i+3] = vrad*sin(t0.vtheta3);
	   t3.vhinit0 = t3.vhinit0+t2.vhinc0;
	   m_tsc.trunk_zparam[i+0] = t3.vhinit0;
	   t3.vhinit1 = t3.vhinit1+t2.vhinc1;
	   m_tsc.trunk_zparam[i+1] = t3.vhinit1;
	   t3.vhinit2 = t3.vhinit2+t2.vhinc2;
	   m_tsc.trunk_zparam[i+2] = t3.vhinit2;
	   t3.vhinit3 = t3.vhinit3+t2.vhinit3;
	   m_tsc.trunk_zparam[i+3] = t3.vhinit3;
       }
}

void
gms::math::TreeScattererAVX512::
SetThicknessDensAng_zmm16r4(const AVX512vEC16 * __restrict bradii) {

#if defined __ICC || defined __INTEL_COMPILER   
    svrng_float16_t vrand1,vrand2,vrand3,vrand4,vrand5,vrand6;
    svrng_engine_t engine;
    svrng_distribution_t uniform1,uniform2,uniform3,uniform4,
                         uniform5;
    uint32_t seed;
    int32_t result;
    int32_t err;
    // Memory first-touch
    gms::common::avx512vec16_init_unroll8x(&m_tsc.leaves_thick[0],
                                       static_cast<int64_t>(m_tsc.nleaves),
				       ZERO);
    gms::common::avx512vec16_init_unroll8x(&m_tsc.leaves_dens[0],
                                       static_cast<int64_t>(m_tsc.nleaves),
				       ZERO);
    result = _rdrand32_step(&seed)
    if(!result) seed = 1458963254U;
    engine = svrng_new_mt19937_engine(seed);
    uniform1 = svrng_new_uniform_distribution_float(0.1f,0.7f);
    uniform2 = svrng_new_uniform_distribution_float(0.1f,0.6f);
    err = 0;
    err = svrng_get_status();
    if(err != SVRNG_STATUS_OK) {
      svrng_delete_engine(engine);
      return;
    }
    __assume_aligned(m_tsc.leaves_thick,64);
    __assume_aligned(m_tsc.leaves_dens,64);
#pragma vector always
#pragma vectorlength(16)
#pragms code_align(32)
#endif
      for(int32_t i = 0; i != m_tsc.nleaves; ++i) {
         vrand1 = svrng_generate16_float(engine,uniform1);
         m_tsc.leaves_thick[i] = *(AVX512Vec16*)&vrand1;
	 vrand2 = svrng_generate16_float(engine,uniform2);
	 m_tsc.leaves_dens[i]  = *(AVX512Vec16*)&vrand2;
     }
      //  gms::common::avx512vec16_init_unroll8x(&m_tsc.leaves_incang[0],
      //                                      static_cast<int64_t>(2*m_tsc.nleaves),
//					    ZERO);
//    uniform3 = svrng_new_uniform_distribution_float(0.3f,0.7f);
//   for(int32_t i = 0; i != 1; ++i) {
//
//        __assume_aligned(m_tsc.leaves_incang,64);
//#pragma vector always
//#pragma vectorlength(16)
//#pragma code_align(32)
//#endif
//        for(int32_t j = 0; j != m_tsc.nleaves; ++j) {
//            vrand3 = svrng_generate16_float(engine,uniform3);
//	    m_tsc.leaves_incang[Ix2D(i,m_tsc.nleaves,j)] = *(AVX512Vec16*)&vrand3;
//	}
//     }
//     uniform4 = svrng_new_uniform_distribution_float(0.75f,1.5f);
//     for(int32_t i = 1; i != 2; ++i) {
//
//        __assume_aligned(m_tsc.leaves_incang,64);
//#pragma vector always
//#pragma vectorlength(16)
//#pragma code_align(32)
//#endif
//        for(int32_t j = 0; i != m_tsc.nleaves; ++j) {
//            vrand4 = svrng_generate16_float(engine,uniform4);
//	    m_tsc.leaves_incang[Ix2D(i,m_tsc.nleaves,j)] = *(AVX512Vec16*)&vrand4;
//	}
//     }
//     gms::common::avx512vec16_init_unroll8x(&m_tsc.branches_incang[0],
//                                            static_cast<int64_t>(2*m_tsc.nbranches),
//				            ZERO);
//     for(int32_t i = 0; i != 1; ++i) {
//
//        __assume_aligned(m_tsc.branches_incang,64);
//#pragma vector always
//#pragma vectorlength(16)
//#pragma code_align(32)
//#endif
//        for(int32_t j = 0; j != m_tsc.nbranches; ++j) {
//            vrand5 = svrng_generate16_float(engine,uniform4);
//	    m_tsc.branches_incang[Ix2D(i,m_tsc.nbranches,j)] = *(AVX512Vec16*)&vrand5;
//	}
//     }
//     uniform5 = svrng_new_uniform_distribution_float(0.75f,1.0f);
//     for(int32_t i = 1; i != 2; ++i) {
//
//       __assume_aligned(m_tsc.branches_incang,64);
//#pragma vector always
//#pragma vectorlength(16)
//#endif
//        for(int32_t j = 0; j != m_tsc.nbranches; ++j) {
//            vrand6 = svrng_generate16_float(engine,uniform5);
//	    m_tsc.branches_incang[Ix2D(i,m_tsc.nbranches,j)] = *(AVX512Vec16*)&vrand6;
//	}
//     }
    
#elif defined __GNUC__ && !defined __INTEL_COMPILER
      float * __restrict __ATTR_ALIGN__(64) plthick  = NULL;
      float * __restrict __ATTR_ALIGN__(64) pldense  = NULL;
      // float * __restrict __ATTR_ALIGN__(64) plincang = NULL;
      // float * __restrict __ATTR_ALIGN__(64) pbincang = NULL;
      const int32_t leaves_len = 16*m_tsc.nleaves;
      const int32_t branch_len = 16*m_tsc.nbranches;
     
      plthick  = gms::common::gms_efmalloca(static_cast<size_t>(leaves_len),64);
      pldense  = gms::common::gms_efmalloca(static_cast<size_t>(leaves_len),64);
    //  plincang = gms::common::gms_efmalloca(static_cast<size_t>(2*leaves_len),64);
    //  pbincang = gms::common::gms_efmalloca(static_cast<size_t>(2*branch_len),64);
      std::clock_t seed;
      seed = std::clock();
      auto srand1 = std::bind(std::uniform_real_distribution<float>(0.1f,0.7f),
                              std::mt19937(seed));
      auto srand2 = std::bind(std::uniform_real_distribution<float>(0.1f,0.6f),
                              std::mt19937(seed));
      // Will GCC vectorize a srand1(2) functor calls -- hmmm... probably not.
      avx256_init_unroll4x_ps(&plthick[0],
                              static_cast<int64_t>(leaves_len),
			      0.0f);
      avx256_init_unroll4x_ps(&pldense[0],
                              static_cast<int64_t>(leaves_len),
			      0.0f);
      for(int32_t i = 0; i != leaves_len; ++i) {
          float rf1 = srand1();
	  plthick[i] = rf1;
	  float rf2 = srand2();
	  pldense[i] = rf2;
      }
    //  avx256_init_unroll4x_ps(&plincang[0],
    //                          static_cast<int64_t>(2*16*m_tsc.nleaves),
//			      0.0f);
     
 //     auto srand3 = std::bind(std::uniform_real_distribution<float>(0.3f,0.7f),
 //                             std::mt19937(seed));
 //     for(int32_t i = 0; i != 1; ++i) {
 //         for(int32_t j = 0; j != leaves_len; ++j) {
 //             float rf = srand3();
//	      plincang[Ix2D(i,leaves_len,j)] = rf;
//	  }
 //     }
 //     auto srand4 = std::bind(std::uniform_real_distribution<float>(0.75f,1.5f),
 //                             std::mt19937(seed));
 //     for(int32_t i = 1; i != 2; ++i) {
  //        for(int32_t j = 0; j != leaves_len; ++j) {
  //            float rf = srand4();
//	      plincang[Ix2D(i,leaves_len,j)] = rf;
//	  }
  //    }
 //     avx256_init_unroll4x_ps(&pbincang[0],
 //                             static_cast<int64_t>(2*16*m_tsc.nbranches),
//			      0.0f);
 //     for(int32_t i = 0; i != 1; ++i) {
 //         for(int32_t j = 0; j != branch_len; ++j) {
 //             float rf = srand4();
//	      pbincang[Ix2D(i,branch_len,j)] = rf;
//	  }
 //     }
 //     auto srand5 = std::bind(std::uniform_real_distribution<float>(0.75f,1.0f),
 //                            std::mt19937(seed));
 //     for(int32_t i = 1; i != 2; ++i) {
 //         for(int32_t j = 0; j != branch_len; ++j) {
 //             float rf = srand5();
//	      pbincang[Ix2D(i,branch_len,j)] = rf;
//	  }
  //    }
      gms::common::avx512vec16_copy_from_r4(&m_tsc.leaves_thick[0],
                                           &plthick[0],
					   leaves_len);
      gms::common::avx512vec16_copy_from_r4(&m_tsc.leaves_dens[0],
                                           &pldense[0],
					   leaves_len);
  //    gms::common::avx512vec16_copy_from_r4(&m_tsc.mleaves_incang[0],
  //                                         &plincang[0],
//					   2*leaves_len);
 //     gms::common::avx512vec16_copy_from_r4(m_tsc.branches_incang[0],
 //                                          &pbincang[0],
//					   2*branch_len);
      _mm_free(plthick);
      plthick = NULL;
      _mm_free(pldense);
      pldense = NULL;
   //   _mm_free(plincang);
   //   plincang = NULL;
   //   _mm_free(pbincang);
   //   pbincang = NULL;
#endif // End __GNUC__ part
       // ! Density set to 0.0 (must find the exact data)
       //    ! Setting only the radii
       //    ! First touch
       avx512vec16_init_unroll8x(&m_tsc.branches_thick[0],
                             static_cast<int64_t>(m_tsc.nbranches),
			     ZERO);
       avx512vec16_init_unroll8x(&m_tsc.branches_dens[0],
                             static_cast<int64_t>(m_tsc.nbranches),
			     ZERO);
       avx512vec16_copy_unroll8x(&m_tsc.branches_thick[0],
                             &bradii[0],
			     static_cast<int64_t>(m_tsc.nbranches));
#if defined __ICC || defined __INTEL_COMPILER
     svrng_delete_engine(engine);
#endif
   		     
}

#if (SAMPLE_HW_PMC) == 1
    #include "GMS_fast_pmc_access.h"
    #include "GMS_pmc_events_names.h"
    #if (CHECK_L2_LICENSE) == 1
        #include "GMS_avx512_warmup_loops.h"
    #endif
    #include <string.h>
    #include <syslog.h>
#endif

void
gms::math::TreeScattererAVX512::
ComputeLeavesParamEq_zmm16r4( const AVX512Vec16 va,
                              const AVX512Vec16 vb) {
			     






     struct _T0_ {
       AVX512Vec16 vthinc0;
       AVX512Vec16 vthinc1;
       AVX512Vec16 vthinc2;
       AVX512Vec16 vthinc3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
       AVX512Vec16 vtheta0;
       AVX512Vec16 vtheta1;
       AVX512Vec16 vtheta2;
       AVX512Vec16 vtheta3;
     } __ATTR_ALIGN__(64) t1

     struct _T2_ {
       AVX512Vec16 vsqrt;
       AVX512Vec16 vsqrtarg;
       AVX512Vec16 vC;
       AVX512Vec16 tmp;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
       AVX512Vec16 tva;
       AVX512Vec16 tvb;
     } __ATTR_ALIGN__(64) t3;
     AVX512Vec16 vNPTS;
     const int64_t xyparam_len = static_cast<int64_t>(m_tsc.nleaves*m_tsc.leaves_param_npts);
     std::clock_t seedx,seedy;
     // Locals first memory-touch
     t0.vthinc0  = ZERO;
     t0.vthinc1  = ZERO;
     t0.vthinc2  = ZERO;
     t0.vthinc3  = ZERO;
     t1.vtheta0  = ZERO;
     t1.vtheta1  = ZERO;
     t1.vtheta2  = ZERO;
     t1.vtheta3  = ZERO;
     t2.vsqrt    = ZERO;
     t2.vsqrtarg = ZERO;
     t2.vC       = ZERO;
     t2.tmp      = ZERO;
     vNPTS       = ZERO;
     t3.tva      = va;
     t3.tvb      = vb;
     // Memory first touch
     gms::common::avx512vec16_init_unroll8x(&m_tsc.leaves_xparam[0],
                                        xyparam_len,
					ZERO);
     gms::common::avx512vec16_init_unroll8x(&m_tsc.leaves_yparam[0],
                                        xyparam_len,
     					ZERO);
     vNPTS = AVX512Vec16{static_cast<float>(m_tsc.leaves_param_npts)};
#if (SAMPLE_HW_PMC) == 1
     uint64_t cycles_lvl2 = 0xFFFFFFFFFFFFFFFFULL;
     constexpr int32_t rdtscp_latency = 42; //Should be empiracally confirmed.
   #if (CHECK_L2_LICENSE) == 1
         const int32_t niter = 1000000;
	 cycles_lvl2 = rdpmc(4); // assume first PMC is 0
         if(!cycles_lvl2) {
            __m512 a = _mm512_setr_ps(0.1f,0.2f,0.3f,0.4f,
	                              0.5f,0.6f,0.7f,0.8f,
				      0.9f,0.10f,0.11f,0.12f,
				      0.13f,0.14f,0.15f,0.16f);
	    __m512 b = _mm512_setr_ps(0.0f,0.11f,0.12f,0.13f,
	                              0.14f,0.15f,0.16f,0.17f,
				      0.18f,0.19f,0.20f,0.21f,
				      0.22f,0.23f,0.24f,0.25f);
	    __m512 c = _mm512_setr_ps(1.0f,2.0f,3.0f,4.0f,
	                              5.0f,6.0f,7.0f,8.0f,
				      9.0f,10.0f,11.0f,12.0f,
				      13.0f,14.0f,15.0f,16.0f);
	    __m512 volatile d = _mm512_setzero_ps();
	    d = avx512_warmup_loop2_ps(a,b,c,niter);
	 }
    #endif
      TREE_SCATTERER_AVX512_HW_PMC_PROLOGE_BODY
#endif
       for(int32_t i = 0; i != m_tsc.nleaves; ++i) {
           // loop over leaves
           seedx = std::clock();
	   auto rand_x = std::bind(std::uniform_real_distribution<float>(0.1f,1.0f),
	                           std::mt19937(seedx));
	   const float xtemp = rand_x();
	   //
	   t3.tva      = t3.tva + xtemp;
	   seedy = std::clock();
	   auto rand_y = std::bind(std::uniform_real_distribution<float>(0.1f,1.0f),
	                           std::mt19937(seedy));
	   const float ytemp = rand_y();
	   t3.tvb      = t3.tvb + ytemp;
	   t2.vsqrtarg = TWO*(t3.tva*t3.tva+t3.tvb*t3.tvb);
	   t2.vsqrt    = sqrt(t2.vsqrt);
	   t2.vC       = PI*t2.vsqrt;
	   t2.tmp      = t2.vC/vNPTS;
	   t0.vthinc0  = t2.tmp;
	   t0.vthinc0  += VINC0;
	   t0.vthinc1  = t2.tmp;
	   t0.vthinc1  += VINC1;
	   t0.vthinc2  = t2.tmp;
	   t0.vthinc2  += VINC2;
	   t0.vthinc3  = t2.tmp;
	   t0.vthinc3  += VINC3;
	   t1.vtheta0  = ZERO;
	   t1.vtheta1  = ZERO;
	   t1.vtheta2  = ZERO;
	   t1.vtheta3  = ZERO;
           __assume_aligned(m_tsc.leaves_xparam,64);
	   __assume_aligned(m_tsc.leaves_yparam,64);
#pragma vector always
#pragma vector vectorlength(16)
#pragma code_align(32)
	   for(int32_t j = 0; j != m_tsc.leaves_param_npts-3; j += 4) {
               // for each leaf -- do ...
	       t1.vtheta0  = t1.vtheta0+t0.vthinc0;
	       m_tsc.leaves_xparam[Ix2D(i,m_tsc.leaves_param_npts,j+0)] = t3.tva*cos(t1.vtheta0);
	       m_tsc.leaves_yparam[Ix2D(i,m_tsc.leaves_param_npts,j+0)] = t3.tvb*sin(t1.vtheta0);
	       t1.vtheta1  = t1.vtheta1+t0.vthinc1;
	       m_tsc.leaves_xparam[Ix2D(i,m_tsc.leaves_param_npts,j+1)] = t3.tva*cos(t1.vtheta1);
	       m_tsc.leaves_yparam[Ix2D(i,m_tsc.leaves_param_npts,j+1)] = t3.tvb*sin(t1.vtheta1);
	       t1.vtheta2  = t1.vtheta2+t0.vthinc2
	       m_tsc.leaves_xparam[Ix2D(i,m_tsc.leaves_param_npts,j+2)] = t3.tva*cos(t1.vtheta2);
	       m_tsc.leaves_yparam[Ix2D(i,m_tsc.leaves_param_npts,j+2)] = t3.tvb*sin(t1.vtheta2);
	       t1.vtheta3  = t1.vtheta3+t0.vthinc3;
	       m_tsc.leaves_xparam[Ix2D(i,m_tsc.leaves_param_npts,j+3)] = t3.tva*cos(t1.vtheta3);
	       m_tsc.leaves_yparam[Ix2D(i,m_tsc.leaves_param_npts,j+3)] = t3.tva*sin(t1.vtheta3);
	   }	   

     }
#if (SAMPLE_HW_PMC) == 1
          
	  TREE_SCATTERER_AVX512_HW_PMC_EPILOGE_BODY  
	     
#endif  
}      
    

   



bool
gms::math::TreeScattererAVX512::
ComputeBranchParamEq_zmm16r4( const int32_t nzpts) {
                             

			     



    
     struct _T0_ {
       AVX512Vec16 vtheta0;
       AVX512Vec16 vtheta1;
       AVX512Vec16 vtheta2;
       AVX512Vec16 vtheta3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
       AVX512Vec16 vthinc0;
       AVX512Vec16 vthinc1;
       AVX512Vec16 vthinc2;
       AVX512Vec16 vthinc3;
     } __ATTR_ALIGN__(64) t1;

     struct _T2_ {
       AVX512Vec16 vhinc0;
       AVX512Vec16 vhinc1;
       AVX512Vec16 vhinc2;
       AVX512Vec16 vhinc3;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
       AVX512Vec16 vhinit0;
       AVX512Vec16 vhinit1;
       AVX512Vec16 vhinit2;
       AVX512Vec16 vhinit3;
     } __ATTR_ALIGN__(64) t3;

     struct _T4_ {
       AVX512Vec16 vNPTS;
       AVX512Vec16 vNZPTS;
       AVX512Vec16 tmp1;
       AVX512Vec16 tmp2;
     } __ATTR_ALIGN__(64) t4;

     struct _T5_ {
       AVX512Vec16 tvrad;
       AVX512Vec16 tvz;
     } __ATTR_ALIGN__(64) t5;
     const AVX512Vec16 rScale{100.0f}; // unit of mm.
     const AVX512Vec16 hScale{300.0f}; // unit of cm.
     const int64_t xyznpts = static_cast<int64_t>(m_tsc.nbranches*m_tsc.branches_param_npts);
     std::clock_t seedr,seedz;
     //Locals first-touch
     t0.vtheta0 = ZERO;
     t0.vtheta1 = ZERO;
     t0.vtheta2 = ZERO;
     t0.vtheta3 = ZERO;
     t1.vthinc0 = ZERO;
     t1.vthinc1 = ZERO;
     t1.vthinc2 = ZERO;
     t1.vthinc3 = ZERO;
     t2.vhinc0  = ZERO;
     t2.vhinc1  = ZERO;
     t2.vhinc2  = ZERO;
     t2.vhinc3  = ZERO;
     t3.vhinit0 = ZERO;
     t3.vhinit1 = ZERO;
     t3.vhinit2 = ZERO;
     t3.vhinit3 = ZERO;
     t4.vNPTS   = ZERO;
     t4.vNZPTS  = ZERO;
     t4.tmp1    = ZERO;
     t4.tmp2    = ZERO;
     t4.vNPTS   = AVX512Vec16{static_cast<float>(m_tsc.branches_param_npts)};
     t4.vNZPTS  = AVX512Vec16{static_cast<float>(m_tsc.branches_param_npts+nzpts)};
     t4.tmp1    = PI/t4.vNPTS;
     t5.tvrad   = ZERO;
     t5.tvz     = ZERO;
        // Memory first touch
     gms::common::avx512vec16_init_unroll8x(&m_tsc.branches_xparam[0],
                                       xyznpts,
				       ZERO);
     t1.vthinc0 = t4.tmp1;
     t1.vthinc0 += VINC0;
     t1.vthinc1 = t4.tmp1;
     t1.vthinc1 += VINC1;
     gms::common::avx512vec16_init_unroll8x(&m_tsc.branches_yparam[0],
                                        xyznpts,
					ZERO);
     t1.vthinc2 = t4.tmp1;
     t1.vthinc2 += VINC2;
     t1.vthinc3 = t4.tmp1;
     t1.vthinc3 += VINC3;
     gms::common::avx512vec16_init_unroll8x(&m_tsc.branches_zparam[0],
                                        xyznpts,
					ZERO);
     t5.tvrad = vrad;
     t5.tvz   = vz;
#if (SAMPLE_HW_PMC) == 1
     uint64_t cycles_lvl2 = 0xFFFFFFFFFFFFFFFFULL;
     constexpr int32_t rdtscp_latency = 42; //Should be empiracally confirmed.
   #if (CHECK_L2_LICENSE) == 1
         const int32_t niter = 1000000;
	 cycles_lvl2 = rdpmc(4); // assume first PMC is 0
         if(!cycles_lvl2) {
            __m512 a = _mm512_setr_ps(0.1f,0.2f,0.3f,0.4f,
	                              0.5f,0.6f,0.7f,0.8f,
				      0.9f,0.10f,0.11f,0.12f,
				      0.13f,0.14f,0.15f,0.16f);
	    __m512 b = _mm512_setr_ps(0.0f,0.11f,0.12f,0.13f,
	                              0.14f,0.15f,0.16f,0.17f,
				      0.18f,0.19f,0.20f,0.21f,
				      0.22f,0.23f,0.24f,0.25f);
	    __m512 c = _mm512_setr_ps(1.0f,2.0f,3.0f,4.0f,
	                              5.0f,6.0f,7.0f,8.0f,
				      9.0f,10.0f,11.0f,12.0f,
				      13.0f,14.0f,15.0f,16.0f);
	    __m512 volatile d = _mm512_setzero_ps();
	    d = avx512_warmup_loop2_ps(a,b,c,niter);
	 }
    #endif
      TREE_SCATTERER_AVX512_HW_PMC_PROLOGE_BODY
#endif   
     for(int32_t i = 0; i != m_tsc.nbranches; ++i) {
         // Loop over branches -- do...
         seedr       = std::clock();
	 auto rand_r = std::bind(std::uniform_real_distribution<float>(0.1f,1.0f),
	                         std::mt19937(seedr));
	 const float rtemp = rand_r();
         t5.tvrad = rScale * rtemp;
	 seedz       = std::clock();
	 auto rand_z = std::bind(std::uniform_real_distribution<float>(0.1f,1.0f),
	                         std::mt19937(seedz));
	 const float ztemp = rand_z();
	

	 t5.tvz   = hScale * ztemp;
	 t4.tmp2  = t5.tvz/t4.vNZPTS;
	 t2.vhinc0 = t4.tmp2;
	 t2.vhinc0 += VINC0;
	 t2.vhinc1 = t4.tmp2;
	 t2.vhinc1 += VINC1;
	 t2.vhinc2 = t4.tmp2;
	 t2.vhinc2 += VINC2;
	 t2.vhinc3 = t4.tmp2;
	 t2.vhinc3 += VINC3;
	 t0.vtheta0 = ZERO;
	 t3.vhinit0 = ZERO;
	 t0.vtheta1 = ZERO;
	 t3.vhinit1 = ZERO;
	 t0.vtheta2 = ZERO;
	 t3.vhinit2 = ZERO;
	 t0.vtheta3 = ZERO;
	 t3.vhinit3 = ZERO;
	 __assume_aligned(m_tsc.branches_xparam,64);
	 __assume_aligned(m_tsc.branches_yparam,64);
	 __assume_aligned(m_tsc.branches_zparam,64);
#pragma vector always
#pragma vector vectorlength(16)
#pragma code_align(32)
           for(int32_t j = 0; j != m_tsc.branches_param_npts-3; j += 4) {
            // for every point -- do...
	       t0.vtheta0 = t0.vtheta0+t1.vthinc0;
	       m_tsc.branches_xparam[Ix2D(i,m_tsc.branches_param_npts,j+0)] = t5.tvrad*cos(t0.vtheta0);
	       m_tsc.branches_yparam[Ix2D(i,m_tsc.branches_param_npts,j+0)] = t5.tvrad*sin(t0.vtheta0);
	       t3.vhinit0 = t3.vhinit0+t2.vhinc0;
	       m_tsc.branches_zparam[Ix2D(i,m_tsc.branches_param_npts,j+0)] = t3.vhinit0;
	       t0.vtheta1 = t0.vtheta1+t1.vthinc1;
	       m_tsc.branches_xparam[Ix2D(i,m_tsc.branches_param_npts,j+1)] = t5.tvrad*cos(t0.vtheta1);
	       m_tsc.branches_yparam[Ix2D(i,m_tsc.branches_param_npts,j+1)] = t5.tvrad*sin(t0.vtheta1);
	       t3.vhinit1 = t3.vhinit1+t2.vhinc1;
	       m_tsc.branches_zparam[Ix2D(i,m_tsc.branches_param_npts,j+1)] = t3.vhinit1;
	       t0.vtheta2 = t0.vtheta2+t1.vthinc2;
	       m_tsc.branches_xparam[Ix2D(i,m_tsc.branches_param_npts,j+2)] = t5.tvrad*cos(t0.vtheta2);
	       m_tsc.branches_yparam[Ix2D(i,m_tsc.branches_param_npts,j+2)] = t5.tvrad*sin(t0.vtheta2);
	       t3.vhinit2 = t3.vhinit2+t2.vhinc2;
	       m_tsc.branches_zparam[Ix2D(i,m_tsc.branches_param_npts,j+2)] = t3.vhinit2;
	       t0.vtheta3 = t0.vtheta3+t1.vthinc3;
	       m_tsc.branches_xparam[Ix2D(i,m_tsc.branches_param_npts,j+3)] = t5.tvrad*cos(t0.vtheta3);
	       m_tsc.branches_yparam[Ix2D(i,m_tsc.branches_param_npts,j+3)] = t5.tvrad*sin(t0.vtheta3);
	       t3.vhinit3 = t3.vhinit3+t2.vhinc3;
	       m_tsc.branches_zparam[Ix2D(i,m_tsc.branches_param_npts,j+3)] = t3.vhinit3;
	     }

       }
#if (SAMPLE_HW_PMC) == 1
            TREE_SCATTERER_AVX512_HW_PMC_EPILOGE_BODY
#endif     
}


void
gms::math::TreeScattererAVX
::ComputeLeafPhaseMatrices(const float * __restrict __ATTR_ALIGN__(64) ldiam,
			   const float * __restrict __ATTR_ALIGN__(64) lthick,
			   const float * __restrict __ATTR_ALIGN__(64) lmg,
			   const float * __restrict __ATTR_ALIGN__(64) lrho,
			   const float * __restrict __ATTR_ALIGN__(64) ldens,
			   const std::complex<float> * __restrict __ATTR_ALIGN__(64) epsr,
			   const int32_t * __restrict __ATTR_ALIGN__(64) lorient,
		           const float theta,
			   const float ctheta,
			   const float stheta,
			   const float rad_freq,
			   const float rad_wv,
			   const float rad_k0) {
 
#if defined __ICC || defined __INTEL_COMPILER
         __assume_aligned(ldiam,64);
	 __assume_aligned(lthick,64);
	 __assume_aligned(lmg,64);
	 __assume_aligned(lrho,64);
	 __assume_aligned(ldens,64);
	 __assume_aligned(epsr,64);
	 __assume_aligned(lorient,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
         ldiam   = (const float*)__builtin_assume_aligned(ldiam,64);
	 lthick  = (const float*)__builtin_assume_aligned(lthick,64);
	 lmg     = (const float*)__builtin_assume_aligned(lmg,64);
	 lrho    = (const float*)__builtin_assume_aligned(lrho,64);
	 ldens   = (const float*)__builtin_assume_aligned(ldens,64);
	 epsr    = (const std::complex<float>*)__builtin_assume_aligned(epsr,64);
	 lorient = (const int32_t*)__builtin_assume_aligned(lorient,64);
#endif
         int32_t i,off64,off16,off4;
         off64 = 64;
	 off16 = 16;
	 off4  = 4;
         // First touch for result arrays
	 gms::common::avx256_init_unroll4x_ps(&m_lp.l4x4phm[0],
                                 static_cast<int64_t>(64*m_tsc.nleaves),
			         0.0f);
	 gms::common::avx256_init_unroll4x_ps(&m_lp.stokes4x4m[0],
                                 static_cast<int64_t>(16*m_tsc.nleaves),
			         0.0f);
	 gms::common::init_unroll8x_cmplxr4(&m_lp.l2x2mp[0],
					   static_cast<int64_t>(4*m_tsc.nleaves),
					   {0.0f,0.0f});
	 gms::common::init_unroll8x_cmplxr4(&m_lp.l2x2mn[0],
					    static_cast<int64_t>(4*m_tsc.nleaves),
					    {0.0f,0.0f});
	 gms::common::init_unroll8x_cmplxr4(&m_lp.scat2x2m[0],
					    static_cast<int64_t>(4*m_tsc.nleaves),
					    {0.0f,0.0f});
 #pragma omp parallel for  schedule(static)   shared(m_lp.l4x4phm,m_lp.l2x2mp,m_lp.l2x2mn,m_lp.stokes4x4m,m_lp.scat2x2m)  private(i,off64,off16,off4)                                                                        
	 for(i=0; i != m_tsc.nleaves; ++i) {
             Leaf_phase_matrices(m_lp.l4x4phm,
	                         m_lp.l2x2mp,
				 m_lp.l2x2mn,
				 m_lp.stokes4x4m,
				 m_lp.scat2x2m,
				 ldiam[i],
				 lthick[i],
				 lmg[i],
				 lrho[i],
				 ldens[i],
				 theta,
				 ctheta,
				 stheta,
				 epsr[i],
				 rad_freq,
				 rad_wv,
				 rad_k0,
				 lorient[i]);
		     m_lp.l4x4phm    += off64;
		     m_lp.l2x2mp     += off4;
		     m_lp.l2x2mn     += off4;
		     m_lp.stokes4x4m += off16;
		     m_lp.scat2x2m   += off4
	 }
}
