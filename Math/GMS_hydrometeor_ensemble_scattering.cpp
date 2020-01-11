
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
#include <math.h>
#include "GMS_hydrometeor_ensemble_scattering.h"
#include "GMS_malloc.h"
#include "GMS_indices.h"
#include "GMS_common.h"

gms::math::HydroMeteorScatterers::
HydroMeteorScatterers() {

     m_hsc.m_np     = -1;
     m_hsc.m_nshpts = -1;
     m_hsc.m_ID     = -1;
     m_hsc.m_nt     = -1;
     m_hsc.m_nxpts  = -1;
     m_hsc.m_nypts  = -1;
     m_hsc.m_nzpts  = -1;
     m_hsc.m_tpv    = 0.0f;
     m_hsc.m_tpsa   = 0.0f;
     m_hsc.m_tpm    = 0.0f;
     m_hsc.m_pcs    = NULL;
     m_hsc.m_radii  = NULL;
     m_hsc.m_pes    = NULL;
     m_hsc.m_ppx    = NULL;
     m_hsc.m_ppy    = NULL;
     m_hsc.m_ppz    = NULL;
     m_hsc.m_type   = " ";
     m_hsc.m_shape  = " ";
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
     m_hsc.m_nshpts  = shpts;
     m_hsc.m_ID    = ID
     m_hsc.m_nt    = nt
     m_hsc.m_nxpts = nxpts;
     m_hsc.m_nypts = nypts
     m_hsc.m_nzpts = nzpts;   
     m_hsc.m_tpv   = 0.0f;
     m_hsc.m_tpsa  = 0.0f;
     m_hsc.m_tpm   = 0.0f;
     m_hsc.m_pcs   = gms_avxvec8_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.m_nshpts),64);
     m_hsc.m_radii = gms_avxvec8_emalloca(static_cast<size_t>(m_hsc.m_np),64);
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

bool
gms::math::HydroMeteorScatterers::
ComputeShape_ymm8r4(
#if defined __ICC || defined __INTEL_COMPILER
AVXVec8 * __restrict cn,
AVXVec8 * __restrict cdef
#elif defined __GNUC__ && !defined __INTEL_COMPILER
float * __restrict cn,
float * __restrict cdef
#endif
                     ) {
#if defined __ICC || defined __INTEL_COMPILER
      static const AVXVec8 Scale10 = AVXVec8{10.0f};
      static const AVXVec8 Scale3  = AVXVec8{3.0f};
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
      AVXVec8 __ATTR_ALIGN__(32) def;
      AVXVec8 __ATTR_ALIGN__(32) n;
      AVXVec8 __ATTR_ALIGN__(32) rad;
      svrng_float8_t vrandn, vrands,vrandd;
      svrng_engine_t enginen,engines,engined;
      svrng_distribution_t uniformn, uniforms,uniformd;
      uint32_t seedn,seeds,seedd;
      int32_t resultn,results,resultd;
      int32_t errn,errs,errd;
      // locals first-touch
      t0.vtheta0 = ZERO;
      t0.vtheta1 = ZERO;
      t0.vtheta2 = ZERO;
      t0.vtheta3 = ZERO;
      t1.vthinc0 = ZERO;
      t1.vthinc1 = ZERO;
      t1.vthinc2 = ZERO;
      t1.vthinc3 = ZERO;
      t2.term0   = ZERO;
      t2.term1   = ZERO;
      t2.term2   = ZERO;
      t2.term3   = ZERO;
      vNPTS      = ZERO;
      vC         = ZERO;
      tmp        = ZERO;
      def        = ZERO;
      n          = ZERO;
      rad        = ZERO;
      //===============
      vNPTS = AVXVec8{static_cast<float>(m_hsc.m_nshpts)};
      // Arrays first-touch
      gms::common::avxvec8_init_unroll8x(&m_hsc.m_pcs[0],
                                         static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nshpts),
					 ZERO);
      gms::common::avxvec8_init_unroll8x_ps(&m_hsc.m_radii[0],
                                           static_cast<int64_t>(m_hsc.m_np),
					   ZERO);
      for(int32_t i = 0; i != m_hsc.m_np; ++i) {
          __assume_aligned(cn,64);
	  __assume_aligned(m_hsc.m_radii,64);
	  __assume_aligned(cdef,64);
          // loop over a coupled particles
          seedn   = 0U;
	  resultn = -9999;
	  errn    = -9999;
	  resultn = _rdrand32_step(&seedn);
	  if(!resultn) seedn = 321564789U;
	  enginen = svrng_new_mt19937_engine(seedn);
	  errn = svrng_get_status();
	  if(errn != SVRNG_STATUS_OK) {
             return (false);
	  }
	  uniformn = svrng_new_uniform_distribution_float(0.1f,1.0f);
	  vrandn   = svrng_generate8_float(enginen,uniformn);
	  cn[i] = *(AVXVec8*)&vrandn;
	  cn[i] *= Scale10;
	  errn = -9999;
	  errn = svrng_get_status();
	  if(errn != SVRNG_STATUS_OK) {
             svrng_delete_engine(enginen);
	     return (false);
	  }
	  seeds   = 0U;
	  results = -9999;
	  errs    = -9999;
	  results = _rdrand32_step(&seeds);
	  if(!results) seeds = 189654123U;
	  engines = svrng_new_mt19937_engine(seeds);
	  errs = svrng_get_status();
	  if(errs != SVRNG_STATUS_OK) {
             return (false);
	  }
	  uniforms = svrng_new_uniform_distribution_float(0.1f,1.0f);
	  vrads    = svrng_generate8_float(engines,uniforms);
	  m_hsc.m_radii[i] = *(AVXVec8*)&vrads;
	  m_hsc.m_radii[i] *= Scale3;
	  errs = -9999;
	  errs = svrng_get_status();
	  if(errs != SVRNG_STATUS_OK) {
             svrng_delete_engine(engines);
	     return (false);
	  }
	  seedd   = 0U;
	  resultd = -9999;
	  errd    = -9999;
	  resultd = _rdrand32_step(&seedd);
	  if(!resultd) seedd = 254987632U;
	  engined = svrng_new_mt19937_engine(seedd);
	  errd = svrng_get_status();
	  if(errd != SVRNG_STATUS_OK) {
             return (false);
	  }
	  uniformd = svrng_new_uniform_distribution_float(0.1f,1.0f);
	  vradd    = svrng_generate8_float(engined,uniformd);
	  cdef[i]  = *(AVXVec8*)&vradd;
	  errd = -9999;
	  errd = svrng_get_status();
	  if(errd != SVRNG_STATUS_OK) {
             svrng_delete_engine(engined);
	     return (false);
	  }
	  vC  = TWO_PI*m_hsc.m_radii[i];
	  tmp = vC/vNTPS;
	  t1.vthinc0 = tmp;
	  t1.vthinc0 *= VINC0;
	  t1.vthinc1 = tmp;
	  t1.vthinc1 *= VINC1;
	  t1.vthinc2 = tmp;
	  t1.vthinc2 *= VINC2;
	  t1.vthinc3 = tmp;
	  t1.vthinc3 *= VINC3;
	  t0.vtheta0 = ZERO;
	  t2.term0   = ZERO;
	  t0.vtheta1 = ZERO;
	  t2.term1   = ZERO;
	  t0.vtheta2 = ZERO;
	  t2.term2   = ZERO;
	  t0.vtheta3 = ZERO;
	  t2.term3   = ZERO;
	  def = cdef[i];
	  n   = cn[i];
	  rad = m_hsc.m_radii[i];
#pragma vector always
#pragma vector vectorlength(8)
          for(int32_t j = 0; j != m_hsc.m_nshpts-3; j += 4) {
              __assume_aligned(m_hsc.m_pcs,64);
              t0.vtheta0 += t1.vthinc0;
	      t2.term0 = ONE+def*cos(n+t0.vtheta0);
	      m_hsc.m_pcs[Ix2D(i,m_hsc.m_nshpts,j+0)] = rad*t2.term0;
	      t0.vtheta1 += t1.vthinc1;
	      t2.term1 = ONE+def*cos(n+t0.vtheta1);
	      m_hsc.m_pcs[Ix2D(i,m_hsc.m_nshpts,j+1)] = rad+t2.term1;
	      t0.vtheta2 += t1.vthinc2;
	      t2.term2 = ONE+def*cos(n+t0.vtheta2);
	      m_hsc.m_pcs[Ix2D(i,m_hsc.m_nshpts,j+2)] = rad+t2.term2;
	      t0.vtheta3 += t1.vthinc3;
	      t2.term3 = ONE+def*cos(n+t0.vtheta3);
	      m_hsc.m_pcs[Ix2D(i,m_hsc.m_nshpts,j+3)] = rad+t2.term3;
	  }
      }
      svrng_delete_engine(enginen);
      svrng_delete_engine(engines);
      svrng_delete_engine(engined);
      return (true);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
       // Scalar version for GCC compiler
       //
       float theta;
       float thinc;
       float inc;
       float term;
       float NPTS;
       float C;
       float tmp;
       float n;
       float def;
       float rad;
       const int64_t totpts = static_cast<int64_t>(8*m_hsc.m_np*8*m_hsc.m_nshpts);
       const int64_t npts   = static_cast<int64_t>(8*m_hsc.m_np);
       const int64_t nshpts = static_cast<int64_t>(8*m_hsc.m_nshpts);
       float * __restrict __ATTR_ALIGN__(64) radii = NULL;
       float * __restrict __ATTR_ALIGN__(64) pcs   = NULL;
       std::clock_t seedn,seeds,seedd;
       // Locals first-touch
       theta = 0.0f;
       thinc = 0.0f;
       inc   = 0.0f;
       term  = 0.0f;
       NPTS  = 0.0f;
       C     = 0.0f;
       tmp   = 0.0f;
       n     = 0.0f;
       def   = 0.0f;
       rad   = 0.0f;
       radii = gms::common::gms_efmalloca(static_cast<size_t>(npts),64);
       pcs   = gms::common::gms_efmalloca(static_cast<size_t>(totpts),64);
       // Arrays first-touch
       gms::common::avx256_init_unroll8x_ps(&radii[0],
                                            npts,
					    0.0f);
       NPTS = static_cast<float>(nshpts);
       inc = 1.0f;
       gms::common::avx256_init_unroll8x_ps(&pcs[0],
                                            totpts,
					    0.0f);
       cn   = (float*)__builtin_assume_aligned(cn,64);
       cdef = (float*)__builtin_assume_aligned(cdef,64);
       for(int32_t i = 0; i != 8*m_hsc.m_np; ++i) {
           seedn = std::clock();
	   auto srandn =  std::bind(std::uniform_real_distribution<float>(0.1f,1.0f),
	                              std::mt19937(seedn);
	   cn[i] = 10.0f*srandn();			      
	   seeds = std::clock();
	   auto srads  =  std::bind(std::uniform_real_distribution<float>(0.1f,1.0f),
	                              std::mt19937(seeds);
	   radii[i] = 3.0f*srads();
	   seedd = std::clock();
	   auto sradd  =  std::bind(std::uniform_real_distribution<float>(0.1f,1.0f),
	                              std::mt19937(seedd);
	   cdef[i] = seedd();
	   C = 6.283185307179586f*radii[i];
	   tmp = C/NPTS;
	   thinc = tmp;
	   thinc *= inc;
	   theta = 0.0f;
	   term  = 0.0f;
	   n     = cn[i];
	   rad   = radii[i];
	   def   = cdef[i];
#pragma omp simd aligned(radii,pcs:64)
           for(int32_t j = 0; j != 8*m_hsc.m_nshpts; ++j) {
               theta += thinc;
	       term  = 1.0f+def*cos(n+theta);
	       pcs[Ix2D(i,8*m_hsc.m_nshpts,j)] = rad*term;
	   }

        }	  
       // Memory first-touch
        gms::common::avxvec8_init_unroll8x(&m_hsc.m_pcs[0],
                                         static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nshpts),
					 ZERO);
        gms::common::avxvec8_init_unroll8x_ps(&m_hsc.m_radii[0],
                                           static_cast<int64_t>(m_hsc.m_np),
					   ZERO);
	 // Copy results back to member arrays.
	for(int32_t i = 0; i != m_hsc.m_np; ++i) {
#pragma omp simd aligned(m_hsc.m_radii,radii:64)
            for(int32_t j = 0; j != 8; ++j) {
                m_hsc.m_radii[i].m256_f32[j] = radii[i*8+j];
	    }
	}
	//
	for(int32_t i = 0; i != m_hsc.m_np*m_hsc.m_nshpts; ++i) {
#pragma omp simd aligned(m_hsc.m_pcs,pcs:64)
            for(int32_t j = 0; j != 8; ++j) {
                m_hsc.m_pcs[i].m256_f32[j] = pcs[i*8+j];
	    }
	}
	_mm_free(radii);
	radii = NULL;
	_mm_free(pcs);
	pcs = NULL;
	return (true);
#endif
}

bool
gms::math::HydroMeteorScatterers::
ComputeEnsembleShape( const float inz,
                      const float incz,
		      const EnsembleShapes type,
		      const float r,
		      const float inphi,
		      const float inth,
		      const float incphi,
		      const float incth,
		      const float sphrad,
		      const float chebn,
		      const float cdeform) {

     float term1,phi,theta,x,y,z,u;
     term1 = 0.0f;
     phi   = 0.0f;
     theta = 0.0f;
     x     = 0.0f;
     y     = 0.0f;
     z     = 0.0f;
     u     = 0.0f;
     switch(type) {

          case EnsembleShapes::Cylindrical : {
               z = inz;
               theta = inth;
#if defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(m_hsc.m_pes,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              m_hsc.m_pes = (float*)__builtin_assume_aligned(m_hsc.m_pes,64);
#endif
	       for(int32_t i = 0; i != 3; ++i) {
                   for(int32_t j = 0; j != m_hsc.m_np; ++j) {
                       theta += incth;
		       z += incz;
		       if(i == 0) {
                          m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = r*cos(theta);
		       }
		       else if(i == 1) {
                          m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = r*sin(theta);
		       }
		       else if(i == 2) {
                          m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = z;
		       }
		   }
	       }
	       break;
	   }
	  case EnsembleShapes::Spheroidal : {
               theta = inth;
	       phi   = inphi;
#if defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(m_hsc.m_pes,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              m_hsc.m_pes = (float*)__builtin_assume_aligned(m_hsc.m_pes,64);
#endif	       
	       for(int32_t i = 0; i != 3; ++i) {
                   for(int32_t j = 0; j != m_hsc.m_np; ++j) {
                       theta += incth;
		       phi += incphi;
		       u = r*cos(phi);
		       if(i == 0) {
                          m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = std::sqrt(r*r-u*u)*cos(theta);
		       }
		       else if(i == 1) {
                          m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = std::sqrt(r*r-u*u)*sin(theta);
		       }
		       else if(i == 2) {
                          m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = u;
		       } 
		   }
	       }
	       break;
	  }
        case EnsembleShapes::Chebyshev : {
             theta = inth;
	     phi = inphi;
#if defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(m_hsc.m_pes,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              m_hsc.m_pes = (float*)__builtin_assume_aligned(m_hsc.m_pes,64);
#endif
              for(int32_t i = 0; i != 3; ++i) {
                  for(int32_t j = 0; j != m_hsc.m_np; ++j) {
                      theta += incth;
		      phi += incphi;
		      term1 = sphrad*(1.0f+cdeform+cos(chebn*theta));
		      if(i == 0) {
                         x = term1*sin(theta)*cos(phi);
			 m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = x;
		      }
		      else if(i == 1) {
                         y = term1*sin(theta)*sin(phi);
			 m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = y;
		      }
		      else if(i == 2) {
                         z = term1*cos(theta);
			  m_hsc.m_pes[Ix2D(i,m_hsc.m_np,j)] = z;
		      }
		  }
	      }
	      break;
	  }
	default : {
                    return (false);
	   }
      }
     return (true);
}

bool
gms::math::HydroMeteorScatterers::
ComputeXparam_ymm8r4(const AVXVec8 * __restrict cn,
                     const AVXVec8 * __restrict cdef) {

     struct _T0_ {
        AVXVec8 vtheta0;
	AVXVec8 vtheta1;
	AVXVec8 vtheta2;
	AVXVec8 vtheta3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
        AVXVec8 vphi0;
	AVXVec8 vphi1;
	AVXVec8 vphi2;
	AVXVec8 vphi3;
     } __ATTR_ALIGN__(64) t1;

     struct _T2_ {
        AVXVec8 vthinc0;
	AVXVec8 vthinc1;
	AVXVec8 vthinc2;
	AVXVec8 vthinc3;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
        AVXVec8 vphinc0;
	AVXVec8 vphinc1;
	AVXVec8 vphinc2;
	AVXVec8 vphinc3;
     } __ATTR_ALIGN__(64) t3;

     struct _T4_ {
        AVXVec8 term0;
	AVXVec8 term1;
	AVXVec8 term2;
	AVXVec8 term3;
     } __ATTR_ALIGN__(64) t4;

     AVXVec8 __ATTR_ALIGN__(32) cn_rand;
     AVXVec8 __ATTR_ALIGN__(32) sphr_rand;
     AVXVec8 __ATTR_ALIGN__(32) cdef_rand;
     AVXVec8 __ATTR_ALIGN__(32) vNPTS;
     AVXVec8 __ATTR_ALIGN__(32) vC;
     AVXVec8 __ATTR_ALIGN__(32) tmp1;
     AVXVec8 __ATTR_ALIGN__(32) tmp2;
     // Locals first-touch
     t0.vtheta0 = ZERO;
     t0.vtheta1 = ZERO;
     t0.vtheta2 = ZERO;
     t0.vtheta3 = ZERO;
     t1.vphi0   = ZERO;
     t1.vphi1   = ZERO;
     t1.vphi2   = ZERO;
     t1.vphi3   = ZERO;
     t2.vthinc0 = ZERO;
     t2.vthinc1 = ZERO;
     t2.vthinc2 = ZERO;
     t2.vthinc3 = ZERO;
     t3.vphinc0 = ZERO;
     t3.vphinc1 = ZERO;
     t3.vphinc2 = ZERO;
     t3.vphinc3 = ZERO;
     t4.term0   = ZERO;
     t4.term1   = ZERO;
     t4.term2   = ZERO;
     t4.term3   = ZERO
     cn_rand    = ZERO;
     sphr_rand  = ZERO;
     cdef_rand  = ZERO;
     vNPTS      = ZERO;
     vC         = ZERO;
     tmp1       = ZERO;
     tmp2       = ZERO;
     vNPTS = AVXVec8{static_cast<float>(m_hsc.m_nxpts)};
     // Array first touch
     gms::common::avxvec8_init_unroll8x(&m_hsc.m_ppx[0],
                                        static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nxpts),
					ZERO);
#if defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(m_hsc.m_ppx,64);
     __assume_aligned(m_hsc.m_radii,64);
     __assume_aligned(cn,64);
     __assume_aligned(cdef,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     m_hsc.m_ppx   = (AVXVec8*)__builtin_assume_aligned(m_hsc.m_ppx,64);
     m_hsc.m_radii = (AVXVec8*)__builtin_assume_aligned(m_hsc.m_radii,64);
     cn            = (AVXVec8*)__builtin_assume_aligned(cn,64);
     cdef          = (AVXVec8*)__builtin_assume_aligned(cdef,64);
#endif
     for(int32_t i = 0; i != m_hsc.m_np; ++i) {
         cn_rand = cn[i];
         sphr_rand = m_hsc.m_radii[i];
	 cdef_rand = cdef[i];
	 vC = TWO_PI*sphr_rand;
	 tmp1 = vC/vNPTS;
	 tmp2 = tmp1;
	 t2.vthinc0 = tmp1;
	 t2.vthinc0 += VINC0;
	 t3.vphinc0 = tmp2;
	 t3.vphinc0 += VINC0;
	 t2.vthinc1 = tmp1;
	 t2.vthinc1 += VINC1;
	 t3.vphinc1 = tmp2;
	 t3.vphinc1 += VINC1;
	 t2.vthinc2 = tmp1;
	 t2.vthinc2 += VINC2;
	 t3.vphinc2 = tmp2;
	 t3.vphinc2 += VINC2;
	 t2.vthinc3 = tmp1;
	 t2.vthinc3 += VINC3;
	 t3.vphinc3 = tmp2;
	 t3.vphinc3 += VINC3;
	 t0.vtheta0 = ZERO;
	 t1.vphi0   = ZERO;
	 t4.term0   = ZERO;
	 t0.vtheta1 = ZERO;
	 t1.vphi1   = ZERO;
	 t4.term1   = ZERO;
	 t0.vtheta2 = ZERO;
	 t1.vphi2   = ZERO;
	 t4.term2   = ZERO;
	 t0.vtheta3 = ZERO;
	 t1.vphi3   = ZERO;
	 t4.term3   = ZERO;
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector always
#pragma vector vectorlength(8)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned( m_hsc.m_ppx,m_hsc.m_radii,cn,cdef:64)
#endif
        for(int32_t j = 0; j != m_hsc.m_nxpts-3; j += 4) {
            t0.vtheta0 += t2.vthinc0;
	    t1.vphi0   += t3.vphinc0;
	    t4.term0   = sphr_rad*(ONE+cdef_rand*cos(cn_rand*t0.vtheta0));
	    t4.term0   = t4.term0*sin(t0.vtheta0)*cos(t1.vphi0);
	    m_hsc.m_ppx[Ix2D(i,m_hsc.m_nxpts,j+0)] = t4.term0;
	    t0.vtheta1 += t2.vthinc1;
	    t1.vphi1   += t3.vphinc1;
	    t4.term1   = sphr_rad*(ONE+cdef_rand*cos(cn_rand*t0.vtheta1));
	    t4.term1   = t4.term1*sin(t0.vtheta1)*cos(t1.vphi1);
	    m_hsc.m_ppx[Ix2D(i,m_hsc.m_nxpts,j+1)] = t4.term1;
	    t0.vtheta2 += t2.vthinc2;
	    t1.vphi2   += t3.vphinc2;
	    t4.term2   = sphr_rad*(ONE+cdef_rand*cos(cn_rand*t0.vtheta2));
	    t4.term2   = t4.term2*sin(t0.vtheta2)*cos(t1.vphi2);
	    m_hsc.m_ppx[Ix2D(i,m_hsc.m_nxpts,j+2)] = t4.term2;
	    t0.vtheta3 += t2.vthinc3;
	    t1.vphi3   += t3.vphinc3;
	    t4.term3   = sphr_rand*(ONE+cdef_rand*cos(cn_rand*t0.vtheta3));
	    t4.term3   = t4.term3*sin(t0.vtheta3)*cos(t1.vphi3);
	    m_hsc.m_ppx[Ix2D(i,m_hsc.m_nxpts,j+3)] = t4.term3;
	}
     }
     
}


				    

				    


