
#if defined __GNUC__ && !defined __INTEL_COMPILER
#include <omp.h>
#endif
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

#include "GMS_indices.h"
#include "GMS_common.h"

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
     m_tsc.leaves_incang       = NULL;
     m_tsc.leaves_xparam       = NULL;
     m_tsc.leaves_yparam       = NULL;
     m_tsc.branches_thick      = NULL;
     m_tsc.branches_dens       = NULL;
     m_tsc.branches_incang     = NULL;
     m_tsc.branches_xparam     = NULL;
     m_tsc.branches_yparam     = NULL;
     m_tsc.branches_zparam     = NULL;
     m_tsh.tree_xangle         = 0.0f;
     m_tsh.tree_yangle         = 0.0f;
     m_tsh.sin_xangle          = 0.0f;
     m_tsh.cos_xangle          = 0.0f;
     m_tsh.sin_yangle          = 0.0f;
     m_tsh.cos_yangle          = 0.0f;
     m_tsh.tree_rcs            = 0.0f;
     m_tsh.crown_rcs           = 0.0f;
     m_tsh.trunk_rcs           = 0.0f;
     m_tsh.leaves_rcs          = NULL;
     m_tsh.leaves_reflect      = NULL;
     m_tsh.branches_rcs        = NULL;
     m_tsh.branches_reflect    = NULL;
     m_tsh.leaves_xang         = NULL;
     m_tsh.leaves_sin_xang     = NULL;
     m_tsh.leaves_cos_xang     = NULL;
     m_tsh.leaves_yang         = NULL;
     m_tsh.leaves_sin_yang     = NULL;
     m_tsh.leaves_cos_yang     = NULL;
     m_tsh.branches_xang       = NULL;
     m_tsh.branches_sin_xang   = NULL;
     m_tsh.branches_cos_xang   = NULL;
     m_tsh.branches_yang       = NULL;
     m_tsh.branches_sin_yang   = NULL;
     m_tsh.branches_cos_yang   = NULL;

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
      m_tsc.leaves_incang       = gms_avx512vec16_emalloca(static_cast<size_t>(2*m_tsc.nleaves),64);                        \
      m_tsc.leaves_xparam       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nleaves*m_tsc.leaves_param_npts),64);  \
      m_tsc.leaves_yparam       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nleaves*m_tsc.leaves_param_npts),64);  \
      m_tsc.branches_thick      = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches),64);                        \
      m_tsc.branches_dens       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches),64);                        \
      m_tsc.branches_incang     = gms_avx512vec16_emalloca(static_cast<size_t>(2*m_tsc.nbranches),64);                      \
      m_tsc.branches_xparam     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches*m_tsc.branches_param_npts),64);   \
      m_tsc.branches_yparam     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches*m_tsc.branches_param_npts),64);   \
      m_tsc.branches_zparam     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nbranches*m_tsc.branches_param_npts),64); 
#endif

#if !defined(GMS_TREE_SCATTERER_AVX512_HOT_ALLOC_CTOR)
    #define GMS_TREE_SCATTERER_AVX512_HOT_ALLOC_CTOR                                                                      \
      m_tsh.leaves_rcs          = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_reflect      = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \ 
      m_tsh.branches_rcs        = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_reflect    = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.leaves_xang         = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_sin_xang     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_cos_xang     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_yang         = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_sin_yang     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_cos_yang     = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.branches_xang       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \ 
      m_tsh.branches_sin_xang   = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_cos_xang   = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_yang       = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_sin_yang   = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_cos_yang   = gms_avx512vec16_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);
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
      m_tsh.tree_xangle         = 0.0f;
      m_tsh.tree_yangle         = 0.0f;
      m_tsh.sin_xangle          = 0.0f;
      m_tsh.cos_xangle          = 0.0f;
      m_tsh.sin_yangle          = 0.0f;
      m_tsh.cos_yangle          = 0.0f;
      m_tsh.tree_rcs            = 0.0f;
      m_tsh.crown_rcs           = 0.0f;
      m_tsh.trunk_rcs           = 0.0f;
      GMS_TREE_SCATTERER_AVX512_HOT_ALLOC_CTOR
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
     _mm_free(m_tsc.leaves_incang);
     m_tsc.leaves_incang   = NULL;
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
     _mm_free(m_tsh.leaves_xang);
     m_tsh.leaves_xang      = NULL;
     _mm_free(m_tsh.leaves_sin_xang);
     m_tsh.leaves_sin_xang  = NULL;
     _mm_free(m_tsh.leaves_cos_xang);
     m_tsh.leaves_cos_xang  = NULL;
     _mm_free(m_tsh.leaves_yang);
     m_tsh.leaves_yang      = NULL;
     _mm_free(m_tsh.leaves_sin_yang);
     m_tsh.leaves_sin_yang  = NULL;
     _mm_free(m_tsh.leaves_cos_yang);
     m_tsh.leaves_cos_yang   = NULL;
     _mm_free(m_tsh.branches_xang);
     m_tsh.branches_xang     = NULL;
     _mm_free(m_tsh.branches_sin_xang);
     m_tsh.branches_sin_xang = NULL;
     _mm_free(m_tsh.branches_cos_xang);
     m_tsh.branches_cos_xang = NULL;
     _mm_free(m_tsh.branches_yang);
     m_tsh.branches_yang     = NULL;
     _mm_free(m_tsh.branches_sin_yang);
     m_tsh.branches_sin_yang = NULL;
     _mm_free(m_tsh.branches_cos_yang);
     m_tsh.branches_cos_yang = NULL;
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
ComputeTrunkParamEq_zmm16r4(const int32_t zpoints) {

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
     gms::common::avx512vec16_init_unroll8x(&m_tsc.leaves_incang[0],
                                            static_cast<int64_t>(2*m_tsc.nleaves),
					    ZERO);
     uniform3 = svrng_new_uniform_distribution_float(0.3f,0.7f);
     for(int32_t i = 0; i != 1; ++i) {

        __assume_aligned(m_tsc.leaves_incang,64);
#pragma vector always
#pragma vectorlength(16)
#pragma code_align(32)
#endif
        for(int32_t j = 0; j != m_tsc.nleaves; ++j) {
            vrand3 = svrng_generate16_float(engine,uniform3);
	    m_tsc.leaves_incang[Ix2D(i,m_tsc.nleaves,j)] = *(AVX512Vec16*)&vrand3;
	}
     }
     uniform4 = svrng_new_uniform_distribution_float(0.75f,1.5f);
     for(int32_t i = 1; i != 2; ++i) {

        __assume_aligned(m_tsc.leaves_incang,64);
#pragma vector always
#pragma vectorlength(16)
#pragma code_align(32)
#endif
        for(int32_t j = 0; i != m_tsc.nleaves; ++j) {
            vrand4 = svrng_generate16_float(engine,uniform4);
	    m_tsc.leaves_incang[Ix2D(i,m_tsc.nleaves,j)] = *(AVX512Vec16*)&vrand4;
	}
     }
     gms::common::avx512vec16_init_unroll8x(&m_tsc.branches_incang[0],
                                            static_cast<int64_t>(2*m_tsc.nbranches),
				            ZERO);
     for(int32_t i = 0; i != 1; ++i) {

        __assume_aligned(m_tsc.branches_incang,64);
#pragma vector always
#pragma vectorlength(16)
#pragma code_align(32)
#endif
        for(int32_t j = 0; j != m_tsc.nbranches; ++j) {
            vrand5 = svrng_generate16_float(engine,uniform4);
	    m_tsc.branches_incang[Ix2D(i,m_tsc.nbranches,j)] = *(AVX512Vec16*)&vrand5;
	}
     }
     uniform5 = svrng_new_uniform_distribution_float(0.75f,1.0f);
     for(int32_t i = 1; i != 2; ++i) {

       __assume_aligned(m_tsc.branches_incang,64);
#pragma vector always
#pragma vectorlength(16)
#endif
        for(int32_t j = 0; j != m_tsc.nbranches; ++j) {
            vrand6 = svrng_generate16_float(engine,uniform5);
	    m_tsc.branches_incang[Ix2D(i,m_tsc.nbranches,j)] = *(AVX512Vec16*)&vrand6;
	}
     }
    
#elif defined __GNUC__ && !defined __INTEL_COMPILER
      float * __restrict __ATTR_ALIGN__(64) plthick  = NULL;
      float * __restrict __ATTR_ALIGN__(64) pldense  = NULL;
      float * __restrict __ATTR_ALIGN__(64) plincang = NULL;
      float * __restrict __ATTR_ALIGN__(64) pbincang = NULL;
      const int32_t leaves_len = 16*m_tsc.nleaves;
      const int32_t branch_len = 16*m_tsc.nbranches;
     
      plthick  = gms::common::gms_efmalloca(static_cast<size_t>(leaves_len),64);
      pldense  = gms::common::gms_efmalloca(static_cast<size_t>(leaves_len),64);
      plincang = gms::common::gms_efmalloca(static_cast<size_t>(2*leaves_len),64);
      pbincang = gms::common::gms_efmalloca(static_cast<size_t>(2*branch_len),64);
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
      avx256_init_unroll4x_ps(&plincang[0],
                              static_cast<int64_t>(2*16*m_tsc.nleaves),
			      0.0f);
     
      auto srand3 = std::bind(std::uniform_real_distribution<float>(0.3f,0.7f),
                              std::mt19937(seed));
      for(int32_t i = 0; i != 1; ++i) {
          for(int32_t j = 0; j != leaves_len; ++j) {
              float rf = srand3();
	      plincang[Ix2D(i,leaves_len,j)] = rf;
	  }
      }
      auto srand4 = std::bind(std::uniform_real_distribution<float>(0.75f,1.5f),
                              std::mt19937(seed));
      for(int32_t i = 1; i != 2; ++i) {
          for(int32_t j = 0; j != leaves_len; ++j) {
              float rf = srand4();
	      plincang[Ix2D(i,leaves_len,j)] = rf;
	  }
      }
      avx256_init_unroll4x_ps(&pbincang[0],
                              static_cast<int64_t>(2*16*m_tsc.nbranches),
			      0.0f);
      for(int32_t i = 0; i != 1; ++i) {
          for(int32_t j = 0; j != branch_len; ++j) {
              float rf = srand4();
	      pbincang[Ix2D(i,branch_len,j)] = rf;
	  }
      }
      auto srand5 = std::bind(std::uniform_real_distribution<float>(0.75f,1.0f),
                             std::mt19937(seed));
      for(int32_t i = 1; i != 2; ++i) {
          for(int32_t j = 0; j != branch_len; ++j) {
              float rf = srand5();
	      pbincang[Ix2D(i,branch_len,j)] = rf;
	  }
      }
      gms::common::avx512vec16_copy_from_r4(&m_tsc.leaves_thick[0],
                                           &plthick[0],
					   leaves_len);
      gms::common::avx512vec16_copy_from_r4(&m_tsc.leaves_dens[0],
                                           &pldense[0],
					   leaves_len);
      gms::common::avx512vec16_copy_from_r4(&m_tsc.mleaves_incang[0],
                                           &plincang[0],
					   2*leaves_len);
      gms::common::avx512vec16_copy_from_r4(m_tsc.branches_incang[0],
                                           &pbincang[0],
					   2*branch_len);
      _mm_free(plthick);
      plthick = NULL;
      _mm_free(pldense);
      pldense = NULL;
      _mm_free(plincang);
      plincang = NULL;
      _mm_free(pbincang);
      pbincang = NULL;
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

void
gms::math::TreeScattererAVX512::
ComputeLeavesParamEq_zmm16r4( const AVX512Vec16 va,
                              const AVX512Vec16 vb){






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
     
}


