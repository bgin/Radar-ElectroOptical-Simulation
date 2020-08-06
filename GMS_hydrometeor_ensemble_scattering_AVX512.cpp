
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
#include "GMS_hydrometeor_ensemble_scattering_AVX512.h"
#include "GMS_malloc.h"
#include "GMS_indices.h"
#include "GMS_common.h"


gms::math::HMScatterersAVX512::
HMScatterersAVX512() {

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

gms::math::HMScatterersAVX512::
HMScatterersAVX512(   const int32_t np,
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
     m_hsc.m_pcs   = gms_avx512vec16_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.m_nshpts),64);
     m_hsc.m_radii = gms_avx512vec16_emalloca(static_cast<size_t>(m_hsc.m_np),64);
     m_hsc.m_pes   = gms_efmalloca(static_cast<size_t>(3*m_hsc.m_np),64);
     m_hsc.m_ppx   = gms_avx512vec16_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.nxpts),64);
     m_hsc.m_ppy   = gms_avx512vec16_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.nypts),64);
     m_hsc.m_ppz   = gms_avx512vec16_emalloca(static_cast<size_t>(m_hsc.m_np*m_hsc.nzpts),64);
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

gms::math::HMScatterersAVX512::
~HMScatterersAVX512() {

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
gms::math::HMScatterersAVX512::
ComputeShape_zmm16r4(
#if defined __ICC || defined __INTEL_COMPILER
   AVX512Vec16 * __restrict cn,
   AVX512Vec16 * __restrict cdef
#elif defined __GNUC__ && !defined __INTEL_COMPILER
   float * __restrict cn,
   float * __restrict cdef
#endif
                     ) {

#if defined __ICC || defined __INTEL_COMPILER
      static const AVX512Vec16 Scale10 = AVX512Vec16{10.0f};
      static const AVX512Vec16 Scale3  = AVX512Vec16{3.0f};
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
         AVX512Vec16 term0;
	 AVX512Vec16 term1;
	 AVX512Vec16 term2;
	 AVX512Vec16 term3;
      } __ATTR_ALIGN__(64) t2;
      AVX512Vec16 __ATTR_ALIGN__(64) vNPTS;
      AVX512Vec16 __ATTR_ALIGN__(64) vC;
      AVX512Vec16 __ATTR_ALIGN__(64) tmp;
      AVX512Vec16 __ATTR_ALIGN__(64) def;
      AVX512Vec16 __ATTR_ALIGN__(64) n;
      AVX512Vec16 __ATTR_ALIGN__(64) rad;
      svrng_float16_t vrandn, vrands,vrandd;
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
      vNPTS = AVX512Vec16{static_cast<float>(m_hsc.m_nshpts)};
      // Arrays first-touch
      gms::common::avx512vec16_init_unroll8x(&m_hsc.m_pcs[0],
                                         static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nshpts),
					 ZERO);
      gms::common::avx512vec16_init_unroll8x_ps(&m_hsc.m_radii[0],
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
	  vrandn   = svrng_generate16_float(enginen,uniformn);
	  cn[i] = *(AVX512Vec16*)&vrandn;
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
	  vrads    = svrng_generate16_float(engines,uniforms);
	  m_hsc.m_radii[i] = *(AVX512Vec16*)&vrads;
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
	  vradd    = svrng_generate16_float(engined,uniformd);
	  cdef[i]  = *(AVX512Vec16*)&vradd;
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
#pragma code_align(32)
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
       const int64_t totpts = static_cast<int64_t>(16*m_hsc.m_np*m_hsc.m_nshpts);
       const int64_t npts   = static_cast<int64_t>(16*m_hsc.m_np);
       const int64_t nshpts = static_cast<int64_t>(16*m_hsc.m_nshpts);
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
       gms::common::avx512_init_unroll8x_ps(&radii[0],
                                            npts,
					    0.0f);
       NPTS = static_cast<float>(nshpts);
       inc = 1.0f;
       gms::common::avx512_init_unroll8x_ps(&pcs[0],
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
	   C     = 6.283185307179586f*radii[i];
	   tmp   = C/NPTS;
	   thinc = tmp;
	   thinc *= inc;
	   theta = 0.0f;
	   term  = 0.0f;
	   n     = cn[i];
	   rad   = radii[i];
	   def   = cdef[i];
#pragma omp simd aligned(radii,pcs:64)
           for(int32_t j = 0; j != 16*m_hsc.m_nshpts; ++j) {
               theta += thinc;
	       term  = 1.0f+def*cos(n+theta);
	       pcs[Ix2D(i,16*m_hsc.m_nshpts,j)] = rad*term;
	   }

        }	  
       // Memory first-touch
       gms::common::avx512vec16_init_unroll8x(&m_hsc.m_pcs[0],
                                         static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nshpts),
					 ZERO);
       gms::common::avx512vec16_init_unroll8x(&m_hsc.m_radii[0],
                                           static_cast<int64_t>(m_hsc.m_np),
					   ZERO);
	 // Copy results back to member arrays.
       gms::common::avx512vec16_copy_from_r4(&m_hsc.m_radii[0],
	                                  &radii[0],
					  npts);
       gms::common::avx512vec16_copy_from_r4(&m_hsc.m_pcs[0],
	                                  &pcs[0],
					  totpts);
  
	_mm_free(radii);
	radii = NULL;
	_mm_free(pcs);
	pcs = NULL;
	return (true);
#endif

}


bool
gms::math::HMScatterersAVX512::
ComputeEnsembleShape( const float inz,
                      const float incz,
		      const EnsembleShapesAVX512 type,
		      const float r,
		      const float inphi,
		      const float inth,
		      const float incphi,
		      const float incth,
		      const float sphrad,
		      const float chebn,
		      const float cdeform ) {

     float term1,phi,theta,x,y,z,u;
     term1 = 0.0f;
     phi   = 0.0f;
     theta = 0.0f;
     x     = 0.0f;
     y     = 0.0f;
     z     = 0.0f;
     u     = 0.0f;
       switch(type) {

          case EnsembleShapesAVX512::Cylindrical : {
               z = inz;
               theta = inth;
#if defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(m_hsc.m_pes,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              m_hsc.m_pes = (float*)__builtin_assume_aligned(m_hsc.m_pes,64);
#endif
	       for(int32_t i = 0; i != 3; ++i) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(m_hsc.m_pes:64)
#elif defined __ICC || defined __INTEL_COMPILER
#pragma vector always
#pragma vector vectorlength(16)
#pragma code_align(32)
#endif
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
	  case EnsembleShapesAVX512::Spheroidal : {
               theta = inth;
	       phi   = inphi;
#if defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(m_hsc.m_pes,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              m_hsc.m_pes = (float*)__builtin_assume_aligned(m_hsc.m_pes,64);
#endif	       
	       for(int32_t i = 0; i != 3; ++i) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(m_hsc.m_pes:64)
#elif defined __ICC || defined __INTEL_COMPILER
#pragma vector always
#pragma vector vectorlength(16)
#pragma code_align(32)
#endif
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
        case EnsembleShapesAVX512::Chebyshev : {
             theta = inth;
	     phi = inphi;
#if defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(m_hsc.m_pes,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              m_hsc.m_pes = (float*)__builtin_assume_aligned(m_hsc.m_pes,64);
#endif
              for(int32_t i = 0; i != 3; ++i) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(m_hsc.m_pes:64)
#elif defined __ICC || defined __INTEL_COMPILER
#pragma vector always
#pragma vector vectorlength(16)
#pragma code_align(32)
#endif	      
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
gms::math::HMScatterersAVX512::
ComputeXparam_zmm16r4(const AVX512Vec16 * __restrict cn,
                      const AVX512Vec16 * __restrict cdef,
		      const char * __restrict pmc_event1,
		      const char * __restrict pmc_event2,
		      const char * __restrict pmc_event3,
		      const char * __restrict pmc_event4) {
     if(__builtin_expect(NULL == cn,0) ||
        __builtin_expect(NULL == cdef,0)) {
           return (false);
     }
 
       struct _T0_ {
        AVX512Vec16 vtheta0;
	AVX512Vec16 vtheta1;
	AVX512Vec16 vtheta2;
	AVX512Vec16 vtheta3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
        AVX512Vec16 vphi0;
	AVX512Vec16 vphi1;
	AVX512Vec16 vphi2;
	AVX512Vec16 vphi3;
     } __ATTR_ALIGN__(64) t1;

     struct _T2_ {
        AVX512Vec16 vthinc0;
	AVX512Vec16 vthinc1;
	AVX512Vec16 vthinc2;
	AVX512Vec16 vthinc3;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
        AVX512Vec16 vphinc0;
	AVX512Vec16 vphinc1;
	AVX512Vec16 vphinc2;
	AVX512Vec16 vphinc3;
     } __ATTR_ALIGN__(64) t3;

     struct _T4_ {
        AVX512Vec16 term0;
	AVX512Vec16 term1;
	AVX512Vec16 term2;
	AVX512Vec16 term3;
     } __ATTR_ALIGN__(64) t4;

     AVX512Vec16 __ATTR_ALIGN__(64) cn_rand;
     AVX512Vec16 __ATTR_ALIGN__(64) sphr_rand;
     AVX512Vec16 __ATTR_ALIGN__(64) cdef_rand;
     AVX512Vec16 __ATTR_ALIGN__(64) vNPTS;
     AVX512Vec16 __ATTR_ALIGN__(64) vC;
     AVX512Vec16 __ATTR_ALIGN__(64) tmp1;
     AVX512Vec16 __ATTR_ALIGN__(64) tmp2;
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
     gms::common::avx512vec16_init_unroll8x(&m_hsc.m_ppx[0],
                                        static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nxpts),
					ZERO);
#if defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(m_hsc.m_ppx,64);
     __assume_aligned(m_hsc.m_radii,64);
     __assume_aligned(cn,64);
     __assume_aligned(cdef,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     m_hsc.m_ppx   = (AVX512Vec16*)__builtin_assume_aligned(m_hsc.m_ppx,64);
     m_hsc.m_radii = (AVX512Vec16*)__builtin_assume_aligned(m_hsc.m_radii,64);
     cn            = (AVX512Vec16*)__builtin_assume_aligned(cn,64);
     cdef          = (AVX512Vec16*)__builtin_assume_aligned(cdef,64);
#endif
#if (SAMPLE_HW_PMC) == 1
           if(pmc_event1 != NULL &&
	      pmc_event2 != NULL &&
	      pmc_event3 != NULL &&
	      pmc_event4 != NULL )       {
            
	      // For now -- only single batch of 4 events is supported
	      const PFC_CNT ZERO_CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CNT CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CFG CFG[7] = {2,2,2,0,0,0,0};
	      CFG[3] = pfcParseCfg(pmc_event1);
	      CFG[4] = pfcParseCfg(pmc_event2);
	      CFG[5] = pfcParseCfg(pmc_event3);
	      CFG[6] = pfcParseCfg(pmc_event4);
	      // Reconfigure PMC and clear their count
	      pfcWrCfgs(0,7,CFG);
	      pfcWrCnts(0,7,ZERO_CNT);
	      memset(CNT,0,sizeof(CNT));
	      // Hot section
	      PFCSTART(CNT);
        }
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
#pragma vector vectorlength(16)
#pragma code_align(32)
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
#if (SAMPLE_HW_PMC) == 1
            PFCEND(CNT);
	    pfcRemoveBias(CNT,1);
	    // Print the results
	    printf("%-10s:\n", __PRETTY_FUNCTION__);
	    printf("Instructions Issued                  : %20lld\n", (signed long long)CNT[0]);
	    printf("Unhalted core cycles                 : %20lld\n", (signed long long)CNT[1]);
	    printf("Unhalted reference cycles            : %20lld\n", (signed long long)CNT[2]);
	    printf("%-37s: %20lld\n", pmc_event1                    , (signed long long)CNT[3]);
	    printf("%-37s: %20lld\n", pmc_event2                    , (signed long long)CNT[4]);
	    printf("%-37s: %20lld\n", pmc_event3                    , (signed long long)CNT[5]);
	    printf("%-37s: %20lld\n", pmc_event4                    , (signed long long)CNT[6]);
#endif
    return (true);
}

bool
gms::math::HMScatterersAVX512::
ComputeYparam_zmm16r4( const AVX512Vec16 * __restrict cn,
		       const AVX512Vec16 * __restrict cdef,
		       const char * __restrict pmc_event1,
		       const char * __restrict pmc_event2,
		       const char * __restrict pmc_event3,
		       const char * __restrict pmc_event4) {

     if(__builtin_expect(NULL == cn,0) ||
        __builtin_expect(NULL == cdef,0)) {
          return (false);
     }
     struct _T0_ {
        AVX512Vec16 vtheta0;
	AVX512Vec16 vtheta1;
	AVX512Vec16 vtheta2;
	AVX512Vec16 vtheta3;
     } __ATTR_ALIGN__(64) t0;

     struct _T1_ {
        AVX512Vec16 vphi0;
	AVX512Vec16 vphi1;
	AVX512Vec16 vphi2;
	AVX512Vec16 vphi3;
     } __ATTR_ALIGN__(64) t1;

     struct _T2_ {
        AVX512Vec16 vthinc0;
	AVX512Vec16 vthinc1;
	AVX512Vec16 vthinc2;
	AVX512Vec16 vthinc3;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
        AVX512Vec16 vphinc0;
	AVX512Vec16 vphinc1;
	AVX512Vec16 vphinc2;
	AVX512Vec16 vphinc3;
     } __ATTR_ALIGN__(64) t3;

     struct _T4_ {
        AVX512Vec16 term0;
	AVX512Vec16 term1;
	AVX512Vec16 term2;
	AVX512Vec16 term3;
     } __ATTR_ALIGN__(64) t4;

     struct _T5_ {
        AVX512Vec16 cn_rand;
	AVX512Vec16 sphr_rand;
	AVX512Vec16 cdef_rand;
     } __ATTR_ALIGN__(64) t5;
     AVX512Vec16 __ATTR_ALIGN__(64) vNPTS;
     AVX512Vec16 __ATTR_ALIGN__(64) vC;
     AVX512Vec16 __ATTR_ALIGN__(64) tmp1;
     AVX512Vec16 __ATTR_ALIGN__(64) tmp2;
     // Automatics first-touch
     t0.vtheta0    = ZERO;
     t0.vtheta1    = ZERO;
     t0.vtheta2    = ZERO;
     t0.vtheta3    = ZERO;
     t1.vphi0      = ZERO;
     t1.vphi1      = ZERO;
     t1.vphi2      = ZERO;
     t1.vphi3      = ZERO;
     t2.vthinc0    = ZERO;
     t2.vthinc1    = ZERO;
     t2.vthinc2    = ZERO;
     t2.vthinc3    = ZERO;
     t3.vphinc0    = ZERO;
     t3.vphinc1    = ZERO;
     t3.vphinc2    = ZERO;
     t3.vphinc3    = ZERO;
     t4.term0      = ZERO;
     t4.term1      = ZERO;
     t4.term2      = ZERO;
     t4.term3      = ZERO
     t5.cn_rand    = ZERO;
     t5.sphr_rand  = ZERO;
     t5.cdef_rand  = ZERO;
     vNPTS         = ZERO;
     vC            = ZERO;
     tmp1          = ZERO;
     tmp2          = ZERO;
     //----------------------//
     vNPTS = AVX512Vec16{static_cast<float>(m_hsc.m_nypts)};
     // Array first-touch
     gms::common::avx512vec16_init_unroll8x(&m_hsc.m_ppy[0],
                                        static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nypts),
					ZERO);
#if defined __ICC || __INTEL_COMPILER
    __assume_aligned(m_hsc.m_ppy,64);
    __assume_aligned(m_hsc.m_radii,64);
    __assume_aligned(cn,64);
    __assume_aligned(cdef,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    m_hsc.m_ppy   = (AVX512Vec16*)__builtin_assume_aligned(m_hsc.m_ppy,64);
    m_hsc.m_radii = (AVX512Vec16*)__builtin_assume_aligned(m_hsc.m_radii,64);
    cn            = (AVX512Vec16*)__builtin_assume_aligned(cn,64);
    cdef          = (AVX512Vec16*)__builtin_assume_aligned(cdef,64);
#endif
#if (SAMPLE_HW_PMC) == 1
           if(pmc_event1 != NULL &&
	      pmc_event2 != NULL &&
	      pmc_event3 != NULL &&
	      pmc_event4 != NULL )       {
            
	      // For now -- only single batch of 4 events is supported
	      const PFC_CNT ZERO_CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CNT CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CFG CFG[7] = {2,2,2,0,0,0,0};
	      CFG[3] = pfcParseCfg(pmc_event1);
	      CFG[4] = pfcParseCfg(pmc_event2);
	      CFG[5] = pfcParseCfg(pmc_event3);
	      CFG[6] = pfcParseCfg(pmc_event4);
	      // Reconfigure PMC and clear their count
	      pfcWrCfgs(0,7,CFG);
	      pfcWrCnts(0,7,ZERO_CNT);
	      memset(CNT,0,sizeof(CNT));
	      // Hot section
	      PFCSTART(CNT);
        }
#endif
      for(int32_t i = 0; i != m_hsc.m_np; ++i) {

         t5.cn_rand   = cn[i];
         t5.sphr_rand = m_hsc.m_radii[i];
	 vC           = TWO_PI*t5.sphr_rand;
	 t5.cdef_rand = cdef[i];
	 tmp1         = vC/vNPTS;
	 tmp2         = tmp1;
	 t2.vthinc0   = tmp1;
	 t2.vthinc0   += VINC0;
	 t3.vphinc0   = tmp2;
	 t3.vphinc0   += VINC0;
	 t2.vthinc1   = tmp1;
	 t2.vthinc1   += VINC1;
	 t3.vphinc1   = tmp2;
	 t3.vphinc1   += VINC1;
	 t2.vthinc2   = tmp1;
	 t2.vthinc2   += VINC2;
	 t3.vphinc2   = tmp2;
	 t3.vphinc2   += VINC2;
	 t2.vthinc3   = tmp1;
	 t2.vthinc3   += VINC3;
	 t3.vphinc3   = tmp2;
	 t3.vphinc3   += VINC3;
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
#pragma vector vectorlength(16)
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(m_hsc.m_ppy,m_hsc.m_radii,cn,cdef:64)
#endif
        for(int32_t j = 0; j != m_hsc.m_nypts-3; j += 4) {

	    t0.vtheta0  += t2.vthinc0;
	    t1.vphi0    += t3.vphinc0;
	    t4.term0    =  t5.sphr_rand*(ONE+t5.cdef_rand*cos(t5.cn_rand*t0.vtheta0));
	    t4.term0    =  t4.term0*sin(t0.vtheta0)*sin(t1.vphi0);
	    m_hsc.m_ppy[Ix2D(i,m_hsc.m_nypts,j+0)] = t4.term0;
	    t0.vtheta1  += t2.vthinc1;
	    t1.vphi1    += t3.vphinc1;
	    t4.term1    =  t5.sphr_rand*(ONE+t5.cdef_rand*cos(t5.cn_rand*t0.vtheta1));
	    t4.term1    =  t4.term1*sin(t0.vtheta1)*sin(t1.vphi1);
	    m_hsc.m_ppy[Ix2D(i,m_hsc.m_nypts,j+1)] = t4.term1;
	    t0.vtheta2  += t2.vthinc2;
	    t1.vphi2    += t3.vphinc2;
	    t4.term2    =  t5.sphr_rand*(ONE+t5.cdef_rand*cos(t5.cn_rand*t0.vtheta2));
	    t4.term2    =  t4.term2*sin(t0.vtheta2)*sin(t1.vphi2);
	    m_hsc.m_ppy[Ix2D(i,m_hsc.m_nypts,j+2)] = t4.term2;
	    t0.vtheta3  += t2.vthinc3;
	    t1.vphi3    += t3.vphinc3;
	    t4.term3    =  t5.sphr_rand*(ONE+t5.cdef_rand*cos(t5.cn_rand*t0.vtheta3));
	    t4.term3    =  t4.term3*sin(t0.vtheta3)*sin(t1.vphi3);
	    m_hsc.m_ppy[Ix2D(i,m_hsc.m_nypts,j+3)] = t4.term3;
	}

    }
#if (SAMPLE_HW_PMC) == 1
            PFCEND(CNT);
	    pfcRemoveBias(CNT,1);
	    // Print the results
	    printf("%-10s:\n", __PRETTY_FUNCTION__);
	    printf("Instructions Issued                  : %20lld\n", (signed long long)CNT[0]);
	    printf("Unhalted core cycles                 : %20lld\n", (signed long long)CNT[1]);
	    printf("Unhalted reference cycles            : %20lld\n", (signed long long)CNT[2]);
	    printf("%-37s: %20lld\n", pmc_event1                    , (signed long long)CNT[3]);
	    printf("%-37s: %20lld\n", pmc_event2                    , (signed long long)CNT[4]);
	    printf("%-37s: %20lld\n", pmc_event3                    , (signed long long)CNT[5]);
	    printf("%-37s: %20lld\n", pmc_event4                    , (signed long long)CNT[6]);
#endif    
    return (true);
				    
}

bool
gms::math::HMScatterersAVX512::
ComputeZparam_zmm16r4(const AVX512Vec16 * __restrict cn,
                      const AVX512Vec16 * __restrict cdef,
		      const char * __restrict pmc_event1,
		      const char * __restrict pmc_event2,
		      const char * __restrict pmc_event3,
		      const char * __restrict pmc_event4) {

    if(__builtin_expect(NULL == cn,0) ||
        __builtin_expect(NULL == cdef,0)) {
          return (false);
     }
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
        AVX512Vec16 term0;
	AVX512Vec16 term1;
	AVX512Vec16 term2;
	AVX512Vec16 term3;
     } __ATTR_ALIGN__(64) t2;

     struct _T3_ {
        AVX512Vec16 cn_rand;
	AVX512Vec16 sphr_rand;
	AVX512Vec16 cdef_rand;
     } __ATTR_ALIGN__(64) t3;

     AVX512Vec16 __ATTR_ALIGN__(64) vC;
     AVX512Vec16 __ATTR_ALIGN__(64) vNPTS;
     AVX512Vec16 __ATTR_ALIGN__(64) tmp;
     //============================//
     // Autoimatics first-touch
     t0.vtheta0   = ZERO;
     t0.vtheta1   = ZERO;
     t0.vtheta2   = ZERO;
     t0.vtheta3   = ZERO;
     t1.vthinc0   = ZERO;
     t1.vthinc1   = ZERO;
     t1.vthinc2   = ZERO;
     t1.vthinc3   = ZERO;
     t2.term0     = ZERO;
     t2.term1     = ZERO;
     t2.term2     = ZERO;
     t2.term3     = ZERO;
     t3.cn_rand   = ZERO;
     t3.sphr_rand = ZERO;
     t3.cdef_rand = ZERO;
     vC           = ZERO;
     vNPTS        = ZERO;
     tmp          = ZERO;
     vNPTS        = AVX512Vec16{static_cast<float>(m_hsc.nzpts)};
     // Array first-touch
     gms::common::avx512vec16_init_unroll8x(&m_hsc.m_ppz[0],
                                        static_cast<int64_t>(m_hsc.m_np*m_hsc.m_nzpts),
					ZERO);
#if defined __ICC || defined __INTEL_COMPILER
    __assume_aligned(m_hsc.m_ppz,64);
    __assume_aligned(m_hsc.m_radii,64);
    __assume_aligned(cn,64);
    __assume_aligned(cdef,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    m_hsc.m_ppz   = (AVX512Vec16*)__builtin_assume_aligned(m_hsc.m_ppz,64);
    m_hsc.m_radii = (AVX512Vec16*)__builtin_assume_aligned(m_hsc.m_radii,64);
    cn            = (AVX512Vec16*)__builtin_assume_aligned(cn,64);
    cdef          = (AVX512Vec16*)__builtin_assume_aligned(cdef,64);
#endif
#if (SAMPLE_HW_PMC) == 1
           if(pmc_event1 != NULL &&
	      pmc_event2 != NULL &&
	      pmc_event3 != NULL &&
	      pmc_event4 != NULL )       {
            
	      // For now -- only single batch of 4 events is supported
	      const PFC_CNT ZERO_CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CNT CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CFG CFG[7] = {2,2,2,0,0,0,0};
	      CFG[3] = pfcParseCfg(pmc_event1);
	      CFG[4] = pfcParseCfg(pmc_event2);
	      CFG[5] = pfcParseCfg(pmc_event3);
	      CFG[6] = pfcParseCfg(pmc_event4);
	      // Reconfigure PMC and clear their count
	      pfcWrCfgs(0,7,CFG);
	      pfcWrCnts(0,7,ZERO_CNT);
	      memset(CNT,0,sizeof(CNT));
	      // Hot section
	      PFCSTART(CNT);
        }
#endif
       for(int32_t i = 0; i != m_hsc.m_np; ++i) {

        t3.cn_rand   = cn[i];
	t3.sphr_rand = m_hsc.m_radii[i];
	vC           = TWO_PI*t3.sphr_rand;
	t3.cdef_rand = cdef[i];
	tmp          = vC/vNPTS;
	t1.vthinc0   = tmp;
	t1.vthinc0   += VINC0;
	t1.vthinc1   = tmp;
	t1.vthinc1   += VINC1;
	t1.vthinc2   = tmp;
	t1.vthinc2   += VINC2;
	t1.vthinc3   = tmp;
	t1.vthinc3   += VINC3;
	t0.vtheta0   = ZERO;
	t2.term0     = ZERO;
	t0.vtheta1   = ZERO;
	t2.term1     = ZERO;
	t0.vtheta2   = ZERO;
	t2.term2     = ZERO;
	t0.vtheta3   = ZERO;
	t2.term3     = ZERO;
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector always
#pragma vector vectorlength(16)
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(m_hsc.m_ppz,m_hsc.m_radii,cn,cdef:64)
#endif
          for(int32_t j = 0; j != m_hsc.m_nzpts-3; j += 4) {

            t0.vtheta0  += t1.vthinc0;
	    t2.term0    =  t3.sphr_rand*(ONE+t3.cdef_rand*cos(t3.cn_rand*t0.vtheta0));
	    t2.term0    =  t2.term0*cos(t0.vtheta0);
	    m_hsc.m_ppz[Ix2D(i,m_hsc.m_nzpts,j+0)] = t2.term0;
	    t0.vtheta1  += t1.vthinc1;
	    t2.term1    =  t3.sphr_rand*(ONE+t3.cdef_rand*cos(t3.cn_rand*t0.vtheta1));
	    t2.term1    =  t2.term1*cos(t0.vtheta1);
	    m_hsc.m_ppz[Ix2D(i,m_hsc.m_nzpts,j+1)] = t2.term1;
	    t0.vtheta2  += t1.vthinc2;
	    t2.term2    =  t3.sphr_rand*(ONE+t3.cdef_rand*cos(t3.cn_rand*t0.vtheta2));
	    t2.term2    =  t2.term2*cos(t0.vtheta2);
	    m_hsc.m_ppz[Ix2D(i,m_hsc.m_nzpts,j+2)] = t2.term2;
	    t0.vtheta3  += t1.vthinc3;
	    t2.term3    =  t3.sphr_rand*(ONE+t3.cdef_rand*cos(t3.cn_rand*vtheta3));
	    t2.term3    =  t3.term3*cos(t0.vtheta3);
	    m_hsc.m_ppz[Ix2D(i,m_hsc.m_nzpts,j+3)] = t2.term3;
	}
     }
#if (SAMPLE_HW_PMC) == 1
            PFCEND(CNT);
	    pfcRemoveBias(CNT,1);
	    // Print the results
	    printf("%-10s:\n", __PRETTY_FUNCTION__);
	    printf("Instructions Issued                  : %20lld\n", (signed long long)CNT[0]);
	    printf("Unhalted core cycles                 : %20lld\n", (signed long long)CNT[1]);
	    printf("Unhalted reference cycles            : %20lld\n", (signed long long)CNT[2]);
	    printf("%-37s: %20lld\n", pmc_event1                    , (signed long long)CNT[3]);
	    printf("%-37s: %20lld\n", pmc_event2                    , (signed long long)CNT[4]);
	    printf("%-37s: %20lld\n", pmc_event3                    , (signed long long)CNT[5]);
	    printf("%-37s: %20lld\n", pmc_event4                    , (signed long long)CNT[6]);
#endif   
     return (true);
}

bool
gms::math::HMScatterersAVX512::
ComputeEnsembleVolume(const float * __restrict sphrad,
                      const float * __restrict chebn,
		      const float * __restrict cdeform,
		      const int32_t totlen) {

     if(__builtin_expect(totlen != (m_hsc.m_np*16),1)) {
        return (false);
     }
     float term1,term1a,term2,term3,term4;
     term1  = 0.0f;
     term1a = 0.0f;
     term2  = 0.0f;
     term3  = 0.0f;
     term4  = 0.0f;
#if defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(sphrad,64);
     __assume_aligned(chebn,64);
     __assume_aligned(cdeform,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     sphrad  = (const float*)__builtin_assume_aligned(sphrad,64);
     chebn   = (const float*)__builtin_assume_aligned(chebn,64);
     cdeform = (const float*)__builtin_assume_aligned(cdeform,64);
#endif
     for(int32_t i = 0; i != totlen; ++i) {
         const float sphrad_p3 = sphrad[i]*sphrad[i]*sphrad[i];
	 term1 = 8.377580409564404f*sphrad_p3;
         const float chebn_p2   = chebn[i]*chebn[i];
	 const float cdeform_p2 = cdeform[i]*cdeform[i];
	 term1a = 1.0f+1.5f*cdeform_p2*(4.0f*chebn_p2-2.0f*0.25f*chebn_p2-1.0f);
	 if((static_cast<int32_t>(chebn[i]) & 1)  == 0) {
	     const float cdeform_p3 = cdeform[i]*cdeform[i]*cdeform[i];
             term2 = 3.0f*cdeform[i]*(1.0f+cdeform_p2*0.25f)/
	             (chebn_p2-1.0f);
	     term3 = 0.25f*cdeform_p3/(9.0f*chebn_p2-1.0f);
	     term4 = term1*(term1a-term2-term3);
	     m_hsc.m_tpv += term4;
	 }
	 else {
              term2 = term1*term1a;
	      m_hsc.m_tpv += term2;
	 }
     }
     return (true);
}

bool
gms::math::HMScatterersAVX512::
ComputeEnsembleSurface(const float * __restrict sphrad,
                       const float * __restrict chebn,
		       const float * __restrict cdeform,
		       const int32_t totlen) {

     if(__builtin_expect(totlen != (m_hsc.m_np*16),1)) {
          return (false);
     }
     float term1,term2,term3,term4,term5,term5a,tmp;
     term1  = 0.0f;
     term2  = 0.0f;
     term3  = 0.0f;
     term4  = 0.0f;
     term5  = 0.0f;
     term5a = 0.0f;
     tmp    = 0.0f;
#if defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(sphrad,64);
     __assume_aligned(chebn,64);
     __assume_aligned(cdeform,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     sphrad  = (const float*)__builtin_assume_aligned(sphrad,64);
     chebn   = (const float*)__builtin_assume_aligned(chebn,64);
     cdeform = (const float*)__builtin_assume_aligned(cdeform,64);
#endif
     for(int32_t i = 0; i != totlen; ++i) {
         const float sphrad_p2 = sphrad[i]*sphrad[i];
         term1 = 25.132741228718346f*sphrad_p2;
	 if((static_cast<int32_t>(chebn[i]) & 1) == 0) {
             const float chebn_p2   = chebn[i]*chebn[i];
	     term2 = 1.0f-2.0f*cdeform[i]/(chebn_p2-1.0f);
	     const float cdeform_p2 = cdeform[i]*cdeform[i];
	     //
	     const float chebn_p4   = chebn_p2*chebn_p2;
	     term3 = cdeform_p2*(chebn_p4+2.0f*chebn_p2-1.0f)/
	             (4.0f*chebn_p2-1.0f);
	     const float cdeform_p4 = cdeform_p2*cdeform_p2;
	     term4 = 3.0f*cdeform_p4*chebn_p4*chebn_p4/
	             (64.0f*chebn_p4-12.0f*chebn_p2-1.0f);
	     term5 = -6.0f*cdeform_p5*cdeform[i]*chebn_p4*chebn_p4;
	     term5a = 1.0f/(chebn_p2-1.0f*9.0f*chebn_p2-1.0f*25.0f*
	              chebn_p2-1.0f);
	     tmp = term1*(term2+term3-term4-term5*term5a);
	     m_hsc.m_tps += tmp;
	 }
	  else {
	      const float chebn_p2   = chebn[i]*chebn[i];
              const float cdeform_p2 = cdeform[i]*cdeform[i];
	      const float chebn_p4   = chebn_p2*chebn_p2;
	      term2 = 1.0f+cdeform_p2*chebn_p4+2.0f*chebn_p2-1.0f/
	              (4.0f*chebn_p2-1.0f);
	      term3 = 3.0f*cdeform_p2*cdeform_p2*chebn_p4*0.015625f;
	      term4 = 1.0f+20.0f*chebn_p2-1.0f/
	              (16.0f*chebn_p2-1.0f*4.0f*chebn_p2-1.0f);
	      tmp = term1 * (term2-term3*term4);
	      m_hsc.m_tps += tmp;
	  }
     }
     return (true);
}

void
gms::math::HMScatterersAVX512::
ComputeEnsembleVfall(const float * __restrict aRe,
                     const float * __restrict bRe,
		     const float vb,
		     const float * __restrict kvisc,
		     const int32_t nx,
		     const int32_t ny,
		     const int32_t nz,
		     const float A,
		     const float rho_b,
		     const float * __restrict rho_f,
		     const float mD,
		     const float Re,
		     const float * __restrict aRet,
		     const float * __restrict bRet) {

     float term1,term2,term2a,term3,inv,t1,t2;
     term1 = 0.0f;
     term2 = 0.0f;
     term2a = 0.0f;
     term3 = 0.0f;
     t1    = 0.0f;
     t2    = 0.0f;
     term2 = (2.0f*vb*9.81f)/A;
#if defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(aRe,64);
     __assume_aligned(bRe,64);
     __assume_aligned(kvisc,64);
     __assume_aligned(rho_f,64);
     __assume_aligned(aRet,64);
     __assume_aligned(bRet,64);
     __assume_aligned(m_hsc.m_pfv,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     aRe         = (const float*)__builtin_assume_aligned(aRe,64);
     bRe         = (const float*)__builtin_assume_aligned(bRe,64);
     kvisc       = (const float*)__builtin_assume_aligned(kvisc,64);
     rho_f       = (const float*)__buitlin_assume_aligned(rho_f,64);
     aRet        = (const float*)__builtin_assume_aligned(aRet,64);
     bRet        = (const float*)__builtin_assume_aligned(bRet,64);
     m_hsc.m_pfv = (float*)__builtin_assume_aligned(m_hsc.m_pfv,64)
#endif
       if((std::abs(Re) - 999.0f) <= std::numeric_limits<float>::epsilon()) {

         for(int32_t i = 0; i != m_hsc.m_nt; ++i) {
             t1 = aRet[i];
	     t2 = bRet[i];
	     tmp = 1.0f-2.0f*t2;
	     for(int32_t ix = 0; ix != nx; ++ix) {
                 for(int32_t iy = 0; iy != ny; ++iy) {
#if defined __INTEL_COMPILER
#pragma vector always
#pragma unroll(2)
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(kvisc,rho_f,aRet,bRet,m_hsc.m_pfv:64)
#endif
                     for(int32_t iz = 0; iz != nz; ++iz) {
		         
                         term1 = t1*std::pow(kvisc[Ix3D(ix,ny,iy,nz,iz)],tmp);
			 if(rho_b > rho_f[Ix3D(ix,ny,iy,nz,iz)]) {
                             term2a = std::abs(rho_b/rho_f[Ix3D(ix,ny,iy,nz,iz)]);
			 }
			 else {
                             term2a = std::abs(rho_b/rho_f[Ix3D(ix,ny,iy,nz,iz)]-1.0f);
			 }
			 
		     }
		 }
	     }
	     term3 = mD*mD*t2-1.0f;
	     m_hsc.m_pfv[i] = term1*std::pow((term2*term2a),t2)-1.0f;
	 }
     }
     else {
           for(int32_t i = 0; i != m_hsc.m_nt; ++i) {
             t1 = aRe[i];
	     t2 = bRe[i];
	     tmp = 1.0f-2.0f*t2;
	     for(int32_t ix = 0; ix != nx; ++ix) {
                 for(int32_t iy = 0; iy != ny; ++iy) {
#if defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#pragma unroll(2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(kvisc,rho_f,aRe,bRe,m_hsc.m_pfv:64)
#endif		 
                     for(int32_t iz = 0; iz != nz; ++iz) {
		         
                         term1 = t1*std::pow(kvisc[Ix3D(ix,ny,iy,nz,iz)],tmp);
			 if(rho_b > rho_f[Ix3D(ix,ny,iy,nz,iz)]) {
                             term2a = std::abs(rho_b/rho_f[Ix3D(ix,ny,iy,nz,iz)]);
			 }
			 else {
                             term2a = std::abs(rho_b/rho_f[Ix3D(ix,ny,iy,nz,iz)]-1.0f);
			 }
			 
		     }
		 }
	     }
	     term3 = mD*mD*t2-1.0f;
	     m_hsc.m_pfv[i] = term1*std::pow((term2*term2a),t2)-1.0f;
	 }
     }
}

#include "GMS_tmatrix_ensemble_iface.h"

void
gms::math::HMScatterersAVX512::
ComputeHydroMeteorScattering(
			     const char * __restrict pmc_event1,
			     const char * __restrict pmc_event2,
			     const char * __restrict pmc_event3,
			     const char * __restrict pmc_event4,
			     int32_t idMie,
			     double  small,
			     int32_t MXINT,
			     int32_t NADD,
			     int32_t idscmt,
			     double  sang,
			     double  w,
			     int32_t irat,
			     int32_t nL,
			     int32_t * __restrict idshp,
			     double  * __restrict shp,
			     double  * __restrict r0) {
#if defined __ICC || defined __INTEL_COMPILER
#if (SAMPLE_HW_PMC) == 1
           if(pmc_event1 != NULL &&
	      pmc_event2 != NULL &&
	      pmc_event3 != NULL &&
	      pmc_event4 != NULL )       {
            
	      // For now -- only single batch of 4 events is supported
	      const PFC_CNT ZERO_CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CNT CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CFG CFG[7] = {2,2,2,0,0,0,0};
	      CFG[3] = pfcParseCfg(pmc_event1);
	      CFG[4] = pfcParseCfg(pmc_event2);
	      CFG[5] = pfcParseCfg(pmc_event3);
	      CFG[6] = pfcParseCfg(pmc_event4);
	      // Reconfigure PMC and clear their count
	      pfcWrCfgs(0,7,CFG);
	      pfcWrCnts(0,7,ZERO_CNT);
	      memset(CNT,0,sizeof(CNT));
	      // Hot section
	      PFCSTART(CNT);
        }
#endif
            mod_tmatrix_mps_mp_tmatrix_mps_driver_(
						  &idMie,
						  &small,
						  &MXINT,
						  &NADD,
						  &idscmt,
						  &sang,
						  &w,
						  &irat,
						  &nL,
						  &idshp[0],
						  &shp[0],
						  &r0[0],
						  &m_hsh.m_cext,
						  &m_hsh.m_cabs,
						  &m_hsh.m_csca,
						  &m_hsh.m_assym,
						  &m_hsh.m_cextv,
						  &m_hsh.m_cabsv,
						  &m_hsh.m_cscav,
						  &m_hsh.m_cbakv,
						  &m_hsh.m_cprv,
						  &m_hsh.m_cexts,
						  &m_hsh.m_cabss,
						  &m_hsh.m_cscas,
						  &m_hsh.m_cbaks,
						  &m_hsh.m_cprs,
						  &m_hsh.m_dang[0],
						  &m_hsh.m_inat[0],
						  &m_hsh.m_pol[0],
						  &m_hsh.m_i11[0],
						  &m_hsh.m_i21[0],
						  &m_hsh.m_i12[0],
						  &m_hsh.m_i22[0],
						  &m_hsh.m_cexti[0],
						  &m_hsh.m_cabsi[0],
						  &m_hsh.m_cscai[0],
						  &m_hsh.m_assymi[0],
						  &m_hsh.m_cpri[0],
						  &m_hsh.m_mue[0]);
#if (SAMPLE_HW_PMC) == 1
            PFCEND(CNT);
	    pfcRemoveBias(CNT,1);
	    // Print the results
	    printf("%-10s:\n", __PRETTY_FUNCTION__);
	    printf("Instructions Issued                  : %20lld\n", (signed long long)CNT[0]);
	    printf("Unhalted core cycles                 : %20lld\n", (signed long long)CNT[1]);
	    printf("Unhalted reference cycles            : %20lld\n", (signed long long)CNT[2]);
	    printf("%-37s: %20lld\n", pmc_event1                    , (signed long long)CNT[3]);
	    printf("%-37s: %20lld\n", pmc_event2                    , (signed long long)CNT[4]);
	    printf("%-37s: %20lld\n", pmc_event3                    , (signed long long)CNT[5]);
	    printf("%-37s: %20lld\n", pmc_event4                    , (signed long long)CNT[6]);
#endif

#elif defined __GNUC__ || defined __GFORTRAN__ && (!defined __INTEL_COMPILER)
#if (SAMPLE_HW_PMC) == 1
           if(pmc_event1 != NULL &&
	      pmc_event2 != NULL &&
	      pmc_event3 != NULL &&
	      pmc_event4 != NULL )       {
            
	      // For now -- only single batch of 4 events is supported
	      const PFC_CNT ZERO_CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CNT CNT[7] = {0,0,0,0,0,0,0};
	      PFC_CFG CFG[7] = {2,2,2,0,0,0,0};
	      CFG[3] = pfcParseCfg(pmc_event1);
	      CFG[4] = pfcParseCfg(pmc_event2);
	      CFG[5] = pfcParseCfg(pmc_event3);
	      CFG[6] = pfcParseCfg(pmc_event4);
	      // Reconfigure PMC and clear their count
	      pfcWrCfgs(0,7,CFG);
	      pfcWrCnts(0,7,ZERO_CNT);
	      memset(CNT,0,sizeof(CNT));
	      // Hot section
	      PFCSTART(CNT);
        }
#endif
            __mod_tmatrix_mps_MOD_tmatrix_mps_driver(
						  &idMie,
						  &small,
						  &MXINT,
						  &NADD,
						  &idscmt,
						  &sang,
						  &w,
						  &irat,
						  &nL,
						  &idshp[0],
						  &shp[0],
						  &r0[0],
						  &m_hsh.m_cext,
						  &m_hsh.m_cabs,
						  &m_hsh.m_csca,
						  &m_hsh.m_assym,
						  &m_hsh.m_cextv,
						  &m_hsh.m_cabsv,
						  &m_hsh.m_cscav,
						  &m_hsh.m_cbakv,
						  &m_hsh.m_cprv,
						  &m_hsh.m_cexts,
						  &m_hsh.m_cabss,
						  &m_hsh.m_cscas,
						  &m_hsh.m_cbaks,
						  &m_hsh.m_cprs,
						  &m_hsh.m_dang[0],
						  &m_hsh.m_inat[0],
						  &m_hsh.m_pol[0],
						  &m_hsh.m_i11[0],
						  &m_hsh.m_i21[0],
						  &m_hsh.m_i12[0],
						  &m_hsh.m_i22[0],
						  &m_hsh.m_cexti[0],
						  &m_hsh.m_cabsi[0],
						  &m_hsh.m_cscai[0],
						  &m_hsh.m_assymi[0],
						  &m_hsh.m_cpri[0],
						  &m_hsh.m_mue[0]);
#if (SAMPLE_HW_PMC) == 1
            PFCEND(CNT);
	    pfcRemoveBias(CNT,1);
	    // Print the results
	    printf("%-10s:\n", __PRETTY_FUNCTION__);
	    printf("Instructions Issued                  : %20lld\n", (signed long long)CNT[0]);
	    printf("Unhalted core cycles                 : %20lld\n", (signed long long)CNT[1]);
	    printf("Unhalted reference cycles            : %20lld\n", (signed long long)CNT[2]);
	    printf("%-37s: %20lld\n", pmc_event1                    , (signed long long)CNT[3]);
	    printf("%-37s: %20lld\n", pmc_event2                    , (signed long long)CNT[4]);
	    printf("%-37s: %20lld\n", pmc_event3                    , (signed long long)CNT[5]);
	    printf("%-37s: %20lld\n", pmc_event4                    , (signed long long)CNT[6]);
#endif 
#endif

}



				    
				    
