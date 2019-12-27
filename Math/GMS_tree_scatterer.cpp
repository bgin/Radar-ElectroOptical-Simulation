
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
#include "GMS_tree_scatterer.h"
//
#include "GMS_malloc.h"

#include "GMS_indices.h"
#include "GMS_common.h"


gms::math::
TreeScatterer::TreeScatterer() {

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

#if !defined(GMS_TREE_SCATTERER_COLD_ALLOC_CTOR)
    #define GMS_TREE_SCATTERER_COLD_ALLOC_CTOR                                                                          \
      m_tsc.leaves_moist        = gms_eimalloca4(static_cast<size_t>(m_tsc.nleaves),64);                                \
      m_tsc.branches_moist      = gms_eimalloca4(static_cast<size_t>(m_tsc.nbranches),64);                              \          
      m_tsc.trunk_xparam        = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.trunk_param_npts),64);                 \
      m_tsc.trunk_yparam        = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.trunk_param_npts),64);                 \
      m_tsc.trunk_zparam        = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.trunk_param_npts),64);                 \
      m_tsc.leaves_thick        = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nleaves),64);                          \
      m_tsc.leaves_dens         = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nleaves),64);                          \
      m_tsc.leaves_incang       = gms_avxvec8_emalloca(static_cast<size_t>(2*m_tsc.nleaves),64);                        \
      m_tsc.leaves_xparam       = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.leaves_param_npts*m_tsc.nleaves),64);  \
      m_tsc.leaves_yparam       = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.leaves_param_npts*m_tsc.nleaves),64);  \
      m_tsc.branches_thick      = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nbranches),64);                        \
      m_tsc.branches_dens       = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nbranches),64);                        \
      m_tsc.branches_incang     = gms_avxvec8_emalloca(static_cast<size_t>(2*m_tsc.nbranches),64);                      \
      m_tsc.branches_xparam     = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.branches_param_npts*m_tsc.nbranches),64);   \
      m_tsc.branches_yparam     = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.branches_param_npts*m_tsc.nbranches),64);   \
      m_tsc.branches_zparam     = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.branches_param_npts*m_tsc.nbranches),64);         
#endif

#if !defined(GMS_TREE_SCATTERER_HOT_ALLOC_CTOR)
    #define GMS_TREE_SCATTERER_HOT_ALLOC_CTOR                                                                         \
      m_tsh.leaves_rcs          = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_reflect      = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \ 
      m_tsh.branches_rcs        = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_reflect    = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.leaves_xang         = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_sin_xang     = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_cos_xang     = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_yang         = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_sin_yang     = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.leaves_cos_yang     = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nleaves),64);           \
      m_tsh.branches_xang       = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \ 
      m_tsh.branches_sin_xang   = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_cos_xang   = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_yang       = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_sin_yang   = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);         \
      m_tsh.branches_cos_yang   = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);
#endif

gms::math::
TreeScatterer::TreeScatterer(const int32_t nleaves,
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
      GMS_TREE_SCATTERER_COLD_ALLOC_CTOR
      m_tsh.tree_xangle         = 0.0f;
      m_tsh.tree_yangle         = 0.0f;
      m_tsh.sin_xangle          = 0.0f;
      m_tsh.cos_xangle          = 0.0f;
      m_tsh.sin_yangle          = 0.0f;
      m_tsh.cos_yangle          = 0.0f;
      m_tsh.tree_rcs            = 0.0f;
      m_tsh.crown_rcs           = 0.0f;
      m_tsh.trunk_rcs           = 0.0f;
      GMS_TREE_SCATTERER_HOT_ALLOC_CTOR
      
}

gms::math::
TreeScatterer::~TreeScatterer() {

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
gms::math::TreeScatterer::
SetMoistness_scalar() {

    
    const uint32_t cutoff_hi = 1U<<16U;
    int32_t result = 0;
    // Memory first touch here!!
    for(int32_t i = 0; i != m_tsc.nleaves; ++i){
        m_tsc.leaves_moist[i] = 0;
    }
    for(int32_t i = 0; i != m_tsc.nbranches; ++i) {
        m_tsc.branches_moist[i] = 0;
    }

    for(int32_t i = 0; i != m_tsc.nleaves; ++i) {
        uint32_t random = 0U;
        result = _rdrand32_step(&random);
        if(!result) {
	   continue;
	}
	else {
	   if(random>cutoff_hi) {
	     m_tsc.leaves_moist[i] = 1;
	   }
	   else{
	     m_tsc.leaves_moist[i] = 0;
	   }
	}
	   
    }
    result = 0;
    for(int32_t i = 0; i != m_tsc.nbranches; ++i) {
        uint32_t random = 0U;
	result = _rdrand32_step(&random);
	if(!result) {
	   continue;
	}
	 else {
             if(random>cutoff_hi) {
                m_tsc.branches_moist[i] = 1;
	     }
	     else {
                m_tsc.branches_moist[i] = 0;
	     }
	 }
    }
}

void
gms::math::TreeScatterer::
ComputeTrunkParamEq_ymm8r4(const int32_t zpoints) {

     static const AVXVec8 twopi = AVXVec8{6.283185307179586f};
     static const AVXVec8 vzero = AVXVec8{};
     
     struct _T0_ {
        AVXVec8 vtheta0;
	AVXVec8 vtheta1;
	AVXVec8 vtheta2;
	AVXVec8 vtheta3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
        AVXVec8 vthinc0;
	AVXVec8 vthinc1;
	AVXVec8 vthinc2;
	AVXVec8 vthinc3;
     } __ATTR_ALIGN__(64) t1;

     struct _T2_ {
        AVXVec8 vhinc0;
	AVXVec8 vhinc1;
	AVXVec8 vhinc2;
	AVXVec8 vhinc3;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
        AVXVec8 vhinit0;
	AVXVec8 vhinit1;
	AVXVec8 vhinit2;
	AVXVec8 vhinit3;
     } __ATTR_ALIGN__(64) t3;

     AVXVec8 tmp1, tmp2;
     AVXVec8 vrad, vNPTS;
     // Locals first-touch
     t0.vtheta0 = vzero;
     t0.vtheta1 = vzero;
     t0.vtheta2 = vzero;
     t0.vtheta3 = vzero;
     t1.vthinc0 = vzero;
     t1.vthinc1 = vzero;
     t1.vthinc2 = vzero;
     t1.vthinc3 = vzero;
     t2.vhinc0  = vzero;
     t2.vhinc1  = vzero;
     t2.vhinc2  = vzero;
     t2.vhinc3  = vzero;
     t3.vhinit0 = vzero;
     t3.vhinit1 = vzero;
     t3.vhinit2 = vzero;
     t3.vhinit3 = vzero;
     vrad       = vzero;
     vNPTS      = vzero;
     tmp1       = vzero;
     tmp2       = vzero;
     vNPTS = AVXVec8{static_cast<float>(m_tsc.trunk_param_npts)};
     tmp1  = twopi/vNPTS;
     t1.vthinc0 = tmp1;
     t1.vthinc0 = t1.vthinc0*VINC0;
     t1.vthinc1 = tmp1;
     t1.vthinc1 = t1.vthinc1*VINC1;
     t1.vthinc2 = tmp1;
     t1.vthinc2 = t1.vthinc2*VINC2;
     t1.vthinc3 = tmp1;
     t1.vthinc3 = t1.vthinc3*VINC3;
     vrad = AVXVec8{m_tsc.trunk_radius};
     zpoints = zpoints+m_tsc.trunk_param_npts;
     tmp2 = AVXVec8{height/static_cast<float>(zpoints)};
     t2.vhinc0 = tmp2;
     t2.vhinc0 = t2.vhinc0*VINC0;
     t2.vhinc1 = tmp2;
     t2.vhinc1 = t2.vhinc1*VINC1;
     t2.vhinc2 = tmp2;
     t2.vhinc2 = t2.vhinc2*VINC2;
     t2.vhinc3 = tmp2;
     t2.vhinc3 = t2.vhinc3*VINC3;
     // First memory touch.
     gms::common::avxvec8_init_unroll8x(&m_tsc.trunk_xparam[0],
                                        static_cast<int64_t>(m_tsc.trunk_param_npts),
					vzero);
     gms::common::avxvec8_init_unroll8x(&m_tsc.trunk_yparam[0],
                                        static_cast<int64_t>(m_tsc.trunk_param_npts),
					vzero);
     gms::common::avxvec8_init_unroll8x(&m_tsc.trunk_zparam[0],
                                        static_cast<int64_t>(m_tsc.trunk_param_npts),
					vzero);
#if defined __GNUC__ && !defined __INTEL_COMPILER
     
#pragma omp simd aligned(m_tsc.trunk_xparam,m_tsc.trunk_yparam,m_tsc.trunk_zparam:64)
#elif defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(m_tsc.trunk_xparam,64);
     __assume_aligned(m_tsc.trunk_yparam,64);
     __assume_aligned(m_tsc.trunk_zparam,64);
#pragma vector always
#pragma vector vectorlength(8)
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
gms::math::TreeScatterer::
SetThicknessDensAng_ymm8r4(const AVXVec8 * __restrict bradii) {
    
   
#if defined __ICC || defined __INTEL_COMPILER   
    svrng_float8_t vrand1,vrand2,vrand3,vrand4,vrand5,vrand6;
    svrng_engine_t engine;
    svrng_distribution_t uniform1,uniform2,uniform3,uniform4,
                         uniform5;
    uint32_t seed;
    int32_t result;
    
    // Memory first-touch
    gms::common::avxvec8_init_unroll8x(&m_tsc.leaves_thick[0],
                                       static_cast<int64_t>(m_tsc.nleaves),
				       AVXVec8{});
    gms::common::avxvec8_init_unroll8x(&m_tsc.leaves_dens[0],
                                       static_cast<int64_t>(m_tsc.nleaves),
				       AVXVec8{});
    result = _rdrand_32_step(&seed)
    if(!result) seed = 1458963254U;
    engine = svrng_new_mt19937_engine(seed);
    uniform1 = svrng_new_uniform_distribution_float(0.1f,0.7f);
    uniform2 = svrng_new_uniform_distribution_float(0.1f,0.6f);
#if defined __ICC || defined __INTEL_COMPILER
    __assume_aligned(m_tsc.leaves_thick,64);
    __assume_aligned(m_tsc.leaves_dens,64);
#pragma vector always
#pragma vectorlength(8)
#endif
     for(int32_t i = 0; i != m_tsc.nleaves; ++i) {
         vrand1 = svrng_generate8_float(engine,uniform1);
         m_tsc.leaves_thick[i] = *(AVXVec8*)&vrand1;
	 vrand2 = svrng_generate8_float(engine,uniform2);
	 m_tsc.leaves_dens[i]  = *(AVXVec8*)&vrand2;
     }
     gms::common::avxvec8_init_unroll8x(&m_tsc.leaves_incang[0],
                                        static_cast<int64_t>(2*m_tsc.nleaves),
					AVXVec8{});
     uniform3 = svrng_new_uniform_distribution_float(0.3f,0.7f);
     for(int32_t i = 0; i != 1; ++i) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(m_tsc.leaves_incang,64);
#pragma vector always
#pragma vectorlength(8)
#endif
        for(int32_t j = 0; j != m_tsc.nleaves; ++j) {
            vrand3 = svrng_generate8_float(engine,uniform3);
	    m_tsc.leaves_incang[Ix2D(i,m_tsc.nleaves,j)] = *(AVXVec8*)&vrand3;
	}
     }
     uniform4 = svrng_new_uniform_distribution_float(0.75f,1.5f);
     for(int32_t i = 1; i != 2; ++i) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(m_tsc.leaves_incang,64);
#pragma vector always
#pragma vectorlength(8)
#endif
        for(int32_t j = 0; i != m_tsc.nleaves; ++j) {
            vrand4 = svrng_generate8_float(engine,uniform4);
	    m_tsc.leaves_incang[Ix2D(i,m_tsc.nleaves,j)] = *(AVXVec8*)&vrand4;
	}
     }
     gms::common::avxvec8_init_unroll8x(&m_tsc.branches_incang[0],
                                        static_cast<int64_t>(2*m_tsc.nbranches),
					AVXVec8{});
     for(int32_t i = 0; i != 1; ++i) {
#if defined __ICC || defined __INTEL_COMPILER
        __assume_aligned(m_tsc.branches_incang,64);
#pragma vector always
#pragma vectorlength(8)
#endif
        for(int32_t j = 0; j != m_tsc.nbranches; ++j) {
            vrand5 = svrng_generate8_float(engine,uniform4);
	    m_tsc.branches_incang[Ix2D(i,m_tsc.nbranches,j)] = *(AVXVec8*)&vrand5;
	}
     }
     uniform5 = svrng_new_uniform_distribution_float(0.75f,1.0f);
     for(int32_t i = 1; i != 2; ++i) {
#if defined __ICC || defined __INTEL_COMPILER
       __assume_aligned(m_tsc.branches_incang,64);
#pragma vector always
#pragma vectorlength(8)
#endif
        for(int32_t j = 0; j != m_tsc.nbranches; ++j) {
            vrand6 = svrng_generate8_float(engine,uniform5);
	    m_tsc.branches_incang[Ix2D(i,m_tsc.nbranches,j)] = *(AVXVec8*)&vrand6;
	}
     }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
      float * __restrict __ATTR_ALIGN__(64) plthick  = NULL;
      float * __restrict __ATTR_ALIGN__(64) pldense  = NULL;
      float * __restrict __ATTR_ALIGN__(64) plincang = NULL;
      float * __restrict __ATTR_ALIGN__(64) pbincang = NULL;
      const int32_t leaves_len = 8*m_tsc.nleaves;
      const int32_t branch_len = 8*m_tsc.nbranches;
      plthick  = _mm_malloc(static_cast<size_t>(leaves_len)*sizeof(float),64);
      if(plthick == NULL && leaves_len != 0) {
         std::exit(EXIT_FAILURE);
      }
      pldense = _mm_malloc(static_cast<size_t>(leaves_len)*sizeof(float),64);
      if(pldense == NULL && leaves_len != 0) {
         std::exit(EXIT_FAILURE);
      }
      plincang = _mm_malloc(static_cast<size_t>(2*8*m_tsc.nleaves)*sizeof(float),64);
      if(plincang == NULL) {
         std::exit(EXIT_FAILURE);
      }
      pbincang = _mm_malloc(static_cast<size_t>(2*8*m_tsc.nbranches)*sizeof(float),64);
      if(pbincang == NULL) {
         std::exit(EXIT_FAILURE);
      }
      //
      //
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
                              static_cast<int64_t>(2*8*m_tsc.nleaves),
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
                              static_cast<int64_t>(2*8*m_tsc.nbranches),
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
   
      for(int32_t i = 0; i != m_tsc.nleaves; ++i) {
#pragma omp simd aligned(m_tsc.leaves_thick,plthick,m_tsc.leaves_dens,pldense,64)
          for(int32_t j = 0; j != 8; ++j) {
              m_tsc.leaves_thick[i].m256_f32[j] = plthick[i*8+j];
	      m_tsc.leaves_dens[i].m256_f32[j]  = pldense[i*8+j];
	  }
      }
      for(int32_t i = 0; i != 2*m_tsc.nleaves; ++i) {
#pragma omp simd aligned(m_tsc.leaves_incang,plincang,64)
          for(int32_t j = 0; j != 8; ++j) {
              m_tsc.leaves_incang[i].m256_f32[j] = plincang[i*8+j];
	  }
      }
      for(int32_t i = 0; i != 2*m_tsc.nbranches; ++j) {
#pragma omp simd aligned(m_tsc.branches_incang,pbincang,64)
          for(int32_t j = 0; j != 8; ++j) {
              m_tsc.branches_incang[i].m256_f32[j] = pbincang[i*8+j];
	  }
      }
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
       avxvec8_init_unroll8x(&m_tsc.branches_thick[0],
                             static_cast<int64_t>(m_tsc.nbranches),
			     AVXVec8{});
       avxvec8_init_unroll8x(&m_tsc.branches_dens[0],
                             static_cast<int64_t>(m_tsc.nbranches),
			     AVXVec8{});
       avxvec8_copy_unroll8x(&m_tsc.branches_thick[0],
                             &bradii[0],
			     static_cast<int64_t>(m_tsc.nbranches));

}

void
gms::math::TreeScatterer::
ComputeLeavesParamEq_ymm8r4(const AVXVec8 va,
                            const AVXVec8 vb) {

     struct _T0_ {
       AVXVec8 vthinc0;
       AVXVec8 vthinc1;
       AVXVec8 vthinc2;
       AVXVec8 vthinc3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
       AVXVec8 vtheta0;
       AVXVec8 vtheta1;
       AVXVec8 vtheta2;
       AVXVec8 vtheta3;
     } __ATTR_ALIGN__(64) t1

     struct _T2_ {
       AVXVec8 vsqrt;
       AVXVec8 vsqrtarg;
       AVXVec8 vC;
       AVXVec8 tmp;
     } __ATTR_ALIGN__(64) t2;
#if defined __ICC || defined __INTEL_COMPILER
     svrng_float8_t vrand1, vrand2;
     svrng_engine_t engine;
     svrng_distribution_t uniform1, uniform2;
     uint32_t seed;
     int32_t result;
     // Locals first memory-touch
     t0.vthinc0 = ZERO;
     t0.vthinc1 = ZERO;
}
