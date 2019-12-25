

#include "GMS_tree_scatterer.h"
//
#include "GMS_malloc.h"
#include "GMS_avxvecf32.h"
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
      m_tsh.branches_sin_yang   = gms_avxvec8_emalloca(static_cast<size_t>(m_tsc.nsteps*m_tsc.nbranches),64);
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
