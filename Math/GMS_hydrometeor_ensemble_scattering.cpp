
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

#include "GMS_hydrometeor_ensemble_scattering.h"
#include "GMS_malloc.h"
#include "GMS_indices.h"
#include "GMS_common.h"

gms::math::HydroMeteorScatterers::
HydroMeteorScatterers() {

     m_hsc.m_np    = -1;
     m_hsc.m_ID    = -1;
     m_hsc.m_nt    = -1;
     m_hsc.m_nxpts = -1;
     m_hsc.m_nypts = -1;
     m_hsc.m_nzpts = -1;
     m_hsc.m_tpv   = 0.0f;
     m_hsc.m_tpsa  = 0.0f;
     m_hsc.m_tpm   = 0.0f;
     m_hsc.m_pcs   = NULL;
     m_hsc.m_radii = NULL;
     m_hsc.m_pes   = NULL;
     m_hsc.m_ppx   = NULL;
     m_hsc.m_ppy   = NULL;
     m_hsc.m_ppz   = NULL;
     m_hsc.m_type  = " ";
     m_hsc.m_shape = " ";
     m_hsh.m_dang[1856] = {};
     m_hsh.m_imat[1856] = {};
     m_hsh.m_pol[1856]  = {};
     m_hsh.m_i11[1856]  = {};
     m_hsh.m_i21[1856]  = {};
     m_hsh.m_i12[1856]  = {};
     m_hsh.m_i22[1856]  = {};
     m_hsh.m_mue[4*4*1856] = {};
     m_hsh.m_cext  = 0.0;
     m_hsh.m_cabs  = 0.0;
     m_hsh.m_csca  = 0.0;
     m_hsh.m_assym = 0.0;
     m_hsh.m_cextv = 0.0;
     m_hsh.m_cabsv = 0.0;
     m_hsh.m_cscav = 0.0;
     m_hsh.m_cbakv = 0.0;
     m_hsh.m_cprv  = 0.0;
     m_hsh.m_cexts = 0.0;
     m_hsh.m_cabss = 0.0;
     m_hsh.m_cscas = 0.0;
     m_hsh.m_cbaks = 0.0;
     m_hsh.m_cprs  = 0.0;
     m_hsh.m_cexti = NULL;
     m_hsh.m_cabsi = NULL;
     m_hsh.m_cscai = NULL;
     m_hsh.m_assymi = NULL;
     m_hsh.m_cpri  = NULL;
     m_hsh.m_rdist = NULL;
     m_hsh.m_theta = NULL;
     m_hsh.m_phi   = NULL;
     m_hsh.m_vfall = NULL;
}

gms::math::HydroMeteorScatterers::
HydroMeteorScatterers(const int32_t np,
                      const int32_t shpts,
                      const int32_t ID,
		      const int32_t nt,
		      const int32_t nxpts,
		      const int32_t nypts,
		      const int32_t nzpts,
                      const char * type,
		      const char * shape) {
     using namespace gms::common;
     m_hsc.m_np    = np
     m_hsc.nshpts  = shpts;
     m_hsc.m_ID    = ID
     m_hsc.m_nt    = nt
     m_hsc.m_nxpts = nxpts;
     m_hsc.m_nypts = nypts
     m_hsc.m_nzpts = nzpts;   
     m_hsc.m_tpv   = 0.0f;
     m_hsc.m_tpsa  = 0.0f;
     m_hsc.m_tpm   = 0.0f;
     m_hsc.m_pcs   = gms_avxvec8_emalloca(static_cast<size_t>(m_hsc.m_np*m_nshpts),64);
     m_hsc.m_radii = gms_efmalloca(static_cast<size_t>(m_hsc.m_np),64);
     m_hsc.m_pes   = gms_efmalloca(static_cast<size_t>(3*m_hsc.m_np),64);
     m_hsc.m_ppx   = gms_avxvec8_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.nxpts),64);
     m_hsc.m_ppy   = gms_avxvec8_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.nypts),64);
     m_hsc.m_ppz   = gms_avxvec8_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.nzpts),64);
     m_type        = type;
     m_shape       = shape;
     m_hsh.m_dang[1856]    = {};
     m_hsh.m_imat[1856]    = {};
     m_hsh.m_pol[1856]     = {};
     m_hsh.m_i11[1856]     = {};
     m_hsh.m_i21[1856]     = {};
     m_hsh.m_i12[1856]     = {};
     m_hsh.m_i22[1856]     = {};
     m_hsh.m_mue[4*4*1856] = {};
     m_hsh.m_cext   = 0.0;
     m_hsh.m_cabs   = 0.0;
     m_hsh.m_csca   = 0.0;
     m_hsh.m_assym  = 0.0;
     m_hsh.m_cextv  = 0.0;
     m_hsh.m_cabsv  = 0.0;
     m_hsh.m_cscav  = 0.0;
     m_hsh.m_cbakv  = 0.0;
     m_hsh.m_cprv   = 0.0;
     m_hsh.m_cexts  = 0.0;
     m_hsh.m_cabss  = 0.0;
     m_hsh.m_cscas  = 0.0;
     m_hsh.m_cbaks  = 0.0;
     m_hsh.m_cprs   = 0.0;
     m_hsh.m_cexti  = gms_edmalloca(static_cast<size_t>(m_hsc.m_np),64);
     m_hsh.m_cabsi  = gms_edmalloca(static_cast<size_t>(m_hsc.m_np),64);
     m_hsh.m_cscai  = gms_edmalloca(static_cast<size_t>(m_hsc.m_np),64);
     m_hsh.m_assymi = gms_edmalloca(static_cast<size_t>(m_hsc.m_np),64);
     m_hsh.m_cpri   = gms_edmalloca(static_cast<size_t>(m_hsc.m_np),64);
     m_hsh.m_rdist  = gms_efmalloca(static_cast<size_t>(m_hsc.m_nt),64);
     m_hsh.m_theta  = gms_efmalloca(static_cast<size_t>(m_hsc.m_nt),64);
     m_hsh.m_phi    = gms_efmalloca(static_cast<size_t>(m_hsc.m_nt),64);
     m_hsh.m_vfall  = gms_efmalloca(static_cast<size_t>(m_hsc.m_nt),64);
}

gms::math::HydroMeteorScatterers::
~HydroMeteorScatterers() {

    _mm_free(m_hsc.m_pcs);
    m_hsc.m_pcs = NULL;
    _mm_free(m_hsc.m_radii);
    m_hsc.m_radii = NULL;
    _mm_free(m_hsc.m_pes);
    m_hsc.m_pes = NULL;
    _mm_free(m_hsc.m_ppx);
    m_hsc.m_ppx = NULL;
    _mm_free(m_hsc.m_ppy);
    m_hsc.m_ppy = NULL;
    _mm_free(m_hsc.m_ppz);
    m_hsc.m_ppz = NULL;
    m_hsc.m_type = NULL;
    m_hsc.m_shape = NULL;
    _mm_free(m_hsh.m_cexti);
    m_hsh.m_cexti = NULL;
    _mm_free(m_hsh.m_cabsi);
    m_hsh.m_cabsi = NULL;
    _mm_free(m_hsh.m_cscai);
    m_hsh.m_cscai = NULL;
    _mm_free(m_hsh.m_assymi);
    m_hsh.m_assymi = NULL;
    _mm_free(m_hsh.m_cpri);
    m_hsh.m_cpri = NULL;
    _mm_free(m_hsh.m_rdist);
    m_hsh.m_rdist = NULL;
    _mm_free(m_hsh.m_theta);
    m_hsh.m_theta = NULL;
    _mm_free(m_hsh.m_phi);
    m_hsh.m_phi = NULL;
    _mm_free(m_hsh.m_vfall);
    m_hsh.m_vfall = NULL;
}

void
gms::math::HydroMeteorScatterers::
ComputeShape_ymm8r4(float * __restrict cn,
                    float * __restrict cdef ) {

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
         AVXVec8 term0;
	 AVXVec8 term1;
	 AVXVec8 term2;
	 AVXVec8 term3;
      } __ATTR_ALIGN__(64) t2;
      AVXVec8 __ATTR_ALIGN__(32) vNPTS;
      AVXVec8 __ATTR_ALIGN__(32) vC;
      AVXVec8 __ATTR_ALIGN__(32) tmp;
      svrng_float8_t vrandx, vrandy;
      svrng_engine_t enginex,enginey;
      svrng_distribution_t uniformx, uniformy;
      uint32_t seedx,seedy;
      int32_t resultx,resulty;
      // locals first-touch
      t0.vtheta0 = vzero;
      t0.vtheta1 = vzero;
      t0.vtheta2 = vzero;
      t0.vtheta3 = vzero;
      t1.vthinc0 = vzero;
      t1.vthinc1 = vzero;
      t1.vthinc2 = vzero;
      t1.vthinc3 = vzero;
      t2.term0   = vzero;
      t2.term1   = vzero;
      t2.term2   = vzero;
      t2.term3   = vzero;
      vNPTS      = vzero;
      vC         = vzero;
      tmp        = vzero;
      
}


