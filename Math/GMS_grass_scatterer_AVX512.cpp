
#if defined __GNUC__ && !defined __INTEL_COMPILER
#include <omp.h>
#endif


#if defined __GNUC__ || defined __INTEL_COMPILER
#include <random>
#endif
#include <ctime>
#include <cstdlib>
#endif
#include "GMS_grass_scatterer_AVX512.h"
#include "GMS_malloc.h"
#include "GMS_indices.h"
#include "GMS_common.h"
#include "GMS_error_macros.h"


gms::math::
GrassScattererAVX512::GrassScattererAVX512() {

     m_gsc.nplants    = -1;
     m_gsc.nsteps     = -1;
     m_gsc.ordinal    = -1;
     m_gsc.param_npts = -1;
     m_gsc.tot_area   = 0.0f;
     m_gsc.lat        = 367.0f;
     m_gsc.lon        = 367.0f;
     m_gsc.elev       = -1.0f;
     m_gsc.Tah        = 0.0f;
     m_gsc.Tav        = 0.0f;
     m_gsc.epsilon    = {0.0f,0.0f};
     m_gsc.mem_lock   = false;
     m_gsc.A          = NULL;
     m_gsc.moistness  = NULL;
     //
     m_gsc.xparam     = NULL;
     m_gsc.yparam     = NULL;
     m_gsc.zparam     = NULL;
     m_gsh.Polv[96]   = {};
     m_gsh.Polh[96]   = {};
     m_gsh.xang       = NULL;
     m_gsh.sin_xang   = NULL;
     m_gsh.cos_xang   = NULL;
     m_gsh.yang       = NULL;
     m_gsh.sin_yang   = NULL;
     m_gsh.cos_yang   = NULL;

}

gms::math::
GrassScattererAVX512::
GrassScattererAVX512(const int32_t nsteps,
                     const int32_t ordinal,
		     const int32_t param_npts,
		     const float lat,
		     const float lon,
		     const float elev,
		     const std::complex<float> cepsilon,
		     const bool lockm) {
     using namespace gms::common;
     //generate random number of grass cylinders
     // per unit area.
     // Min value is 50 max value is 500
     const uint32_t lo   = 50U;
     const uint32_t hi   = 500U;
     int32_t result      = 0;
     uint32_t rand       = 0U;
     result = _rdrand32_step(&rand);
     if(!result) {
     	  ABORT_ON_ERROR("GrassScattererAVX::GrassScattererAVX -- !!! _rdrand32_step failure !!! ", result)
     }
     if(rand < lo) rand = 50U;
     if(rand > hi) rand = 500U;
     m_gsc.nplants = static_cast<int32_t>(rand);
     m_gsc.nsteps = nsteps;
     m_gsc.ordinal = ordinal;
     m_gsc.param_npts = param_npts;
     m_gsc.tot_area  = 0.0f;
     m_gsc.lat       = lat;
     m_gsc.lon       = lon;
     m_gsc.elev      = elev;
     m_gsc.Tah       = 0.0f;
     m_gsc.Tav       = 0.0f;
     m_gsc.epsilon   = cepsilon;
     m_gsc.mem_lock  = lockm;
#if (USE_MMAP_4KiB) == 1
     m_gsc.A         = gms_efmmap_4KiB(static_cast<size_t>(m_gsc.nplants),PROT_READ | PROT_WRITE,
				       MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsc.moistness = gms_immap_4KiB(static_cast<size_t>(m_gsc.nplants),PROT_READ | PROT_WRITE,
				      MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
#elif (USE_MMAP_2MiB) == 1
     m_gsc.A         = gms_efmmap_2MiB(static_cast<size_t>(m_gsc.nplants),PROT_READ | PROT_WRITE,
				       MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsc.moistness = gms_immap_2MiB(static_cast<size_t>(m_gsc.nplants),PROT_READ | PROT_WRITE,
				      MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB,-1,0,m_gsc.mem_lock);
#elif (USE_MMAP_1GiB) == 1
     m_gsc.A         = gms_efmmap_1GiB(static_cast<size_t>(m_gsc.nplants),PROT_READ | PROT_WRITE,
				       MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB ,-1,0,m_gsc.mem_lock);
     m_gsc.moistness = gms_immap_1GiB(static_cast<size_t>(m_gsc.nplants),PROT_READ | PROT_WRITE,
				      MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
#else
     m_gsc.A         = gms_efmalloca(static_cast<size_t>(m_gsc.nplants),64,m_gsc.mem_lock);
     m_gsc.moistness = gms_eimalloca4(static_cast<size_t>(m_gsc.nplants),64,m_gsc.mem_lock);
#endif
     m_gsc.epsilon   = cepsilon;
#if   (USE_MMAP_4KiB) == 1
     m_gsc.xparam    = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsc.yparam    = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsc.zparam    = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
#elif (USE_MMAP_2MiB) == 1
     m_gsc.xparam    = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsc.yparam    = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsc.zparam    = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
#elif (USE_MMAP_1GiB) == 1
     m_gsc.xparam    = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
     m_gsc.yparam    = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
     m_gsc.zparam    = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.param_npts),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
#else
     m_gsc.xparam    = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.param_npts),64,m_gsc.mem_lock);
     m_gsc.yparam    = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.param_npts),64,m_gsc.mem_lock);
     m_gsc.zparam    = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.param_npts),64,m_gsc.mem_lock);
#endif
     m_gsh.Polv[96]  = {0.0f,0.0f};
     m_gsh.Polh[96]  = {0.0f,0.0f};
#if (USE_MMAP_4KiB) == 1
     m_gsh.xang      = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsh.sin_xang  = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsh.cos_xang  = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsh.yang      = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsh.sin_yang  = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
     m_gsh.cos_yang  = gms_avx512vec16_mmap_4KiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE,-1,0,m_gsc.mem_lock);
#elif (USE_MMAP_2MiB) == 1
     m_gsh.xang      = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsh.sin_xang  = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsh.cos_xang  = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsh.yang      = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsh.sin_yang  = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
     m_gsh.cos_yang  = gms_avx512vec16_mmap_2MiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB ,-1,0,m_gsc.mem_lock);
#elif (USE_MMAP_1GiB) == 1
     m_gsh.xang      = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
     m_gsh.sin_xang  = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
     m_gsh.cos_xang  = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
     m_gsh.yang      = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
     m_gsh.sin_yang  = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
     m_gsh.cos_yang  = gms_avx512vec16_mmap_1GiB(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						 PROT_READ | PROT_WRITE,
				                 MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_1GB,-1,0,m_gsc.mem_lock);
#else
     m_gsh.xang      = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						64,m_gsc.mem_lock);
     m_gsh.sin_xang  = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						64,m_gsc.mem_lock);
     m_gsh.cos_xang  = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						64,m_gsc.mem_lock);
     m_gsh.yang      = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						64,m_gsc.mem_lock);
     m_gsh.sin_yang  = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						64,m_gsc.mem_lock);
     m_gsh.cos_yang  = gms_avx512vec16_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),
						64,m_gsc.mem_lock);
#endif
}

gms::math::GrassScattererAVX512
::~GrassScattererAVX512() {

  if(m_gsc.mem_lock) {
    int32_t dummy;
    dummy = munlock(m_gsc.A,static_cast<size_t>(m_gsc.nplants));
    dummy = munlock(m_gsc.moistness,static_cast<size_t>(m_gsc.nplants));
    dummy = munlock(m_gsc.xparam,static_cast<size_t>(m_gsc.param_npts));
    dummy = munlock(m_gsc.yparam,static_cast<size_t>(m_gsc.param_npts));
    dummy = munlock(m_gsc.zparam,static_cast<size_t>(m_gsc.param_npts));
    dummy = munlock(m_gsh.xang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
    dummy = munlock(m_gsh.sin_xang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
    dummy = munlock(m_gsh.cos_xang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
    dummy = munlock(m_gsh.yang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
    dummy = munlock(m_gsh.sin_yang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
    dummy = munlock(m_gsh.cos_yang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
  }
#if (USE_MMAP_4KiB) == 1 || (USE_MMAP_2MiB) == 1 || (USE_MMAP_1GiB) == 1
     int32_t dummy = 0;
     dummy = munmap(m_gsc.A,static_cast<size_t>(m_gsc.nplants));
     m_gsc.A = NULL;
     dummy = munmap(m_gsc.moistness,static_cast<size_t>(m_gsc.nplants));
     m_gsc.moistness = NULL;
     dummy = munmap(m_gsc.xparam,static_cast<size_t>(m_gsc.param_npts));
     m_gsc.xparam = NULL;
     dummy = munmap(m_gsc.yparam,static_cast<size_t>(m_gsc.param_npts));
     m_gsc.yparam = NULL;
     dummy = munmap(m_gsc.zparam,static_cast<size_t>(m_gsc.param_npts));
     m_gsc.zparam = NULL;
     dummy = munmap(m_gsh.xang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
     m_gsh.xang = NULL;
     dummy = munmap(m_gsh.sin_xang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
     m_gsh.sin_xang = NULL;
     dummy = munmap(m_gsh.cos_xang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
     m_gsh.cos_xang = NULL;
     dummy = munmap(m_gsh.yang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
     m_gsh.yang = NULL;
     dummy = munmap(m_gsh.sin_yang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
     m_gsh.sin_yang = NULL;
     dummy = munmap(m_gsh.cos_yang,static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants));
     m_gsh.cos_yang = NULL;
#else
 
  _mm_free(m_gsc.A);
  m_gsc.A = NULL;
  _mm_free(m_gsc.moistness);
  m_gsc.moistness = NULL;
  _mm_free(m_gsc.xparam);
  m_gsc.xparam = NULL;
  _mm_free(m_gsc.yparam);
  m_gsc.yparam = NULL;
  _mm_free(m_gsc.zparam);
  m_gsc.zparam = NULL;
  _mm_free(m_gsh.xang);
  m_gsh.xang = NULL;
  _mm_free(m_gsh.sin_xang);
  m_gsh.sin_xang = NULL;
  _mm_free(m_gsh.cos_xang);
  m_gsh.cos_xang = NULL;
  _mm_free(m_gsh.yang);
  m_gsh.yang = NULL;
  _mm_free(m_gsh.sin_yang);
  m_gsh.sin_yang = NULL;
  _mm_free(m_gsh.cos_yang);
  m_gsh.cos_yang = NULL;
#endif
}

void
gms::math::GrassScattererAVX512
::SetGrassMoistnessMask() {

     std::random_device rd;
  // Memory first touch
  for(int32_t i = 0; i != m_gsc.nplants-7; i += 8) {
      m_gsc.moistness[i+0] = 0U;
      m_gsc.moistness[i+1] = 0U;
      m_gsc.moistness[i+2] = 0U;
      m_gsc.moistness[i+3] = 0U;
      m_gsc.moistness[i+4] = 0U;
      m_gsc.moistness[i+5] = 0U;
      m_gsc.moistness[i+6] = 0U;
      m_gsc.moistness[i+7] = 0U;
  }
  std::mt19937 rgen(rd());
  std::uniform_int_distribution<> distro(0,1);
  for(int32_t i = 0; i != m_gsc.nplants-3; i += 4) {
      m_gsc.moistness[i+0] = distro(rgen);
      m_gsc.moistness[i+1] = distro(rgen);
      m_gsc.moistness[i+2] = distro(rgen);
      m_gsc.moistness[i+3] = distro(rgen);
  }
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
gms::math::GrassScattererAVX512::
ComputeGrassParamEq_zmm16r4() {


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
       //
       AVX512Vec16 tmp1;
       AVX512Vec16 tmp2;
     } __ATTR_ALIGN__(64) t4;

     struct _T5_ {
       AVX512Vec16 tvrad;
       AVX512Vec16 tvz;
     } __ATTR_ALIGN__(64) t5;
     AVX512Vec16 rScale = AVX512Vec16{10.0f}; // unit of mm
     AVX512Vec16 hScale = AVX512Vec16{100.0f}; // unit of cm
     const int64_t totpts = static_cast<int64_t>(m_gsc.nplants*m_gsc.param_npts);
     const static float n2PI = 6.283185307179586f;
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
     
     t4.tmp1    = ZERO;
     t4.tmp2    = ZERO;
     t4.vNPTS   = AVX512Vec16{static_cast<float>(m_tsc.param_npts)};
     
     t4.tmp1    = PI/t4.vNPTS;
     t5.tvrad   = ZERO;
     t5.tvz     = ZERO;
     // Memory first touch
     gms::common::avx512vec16_init_unroll8x(&m_gsc.xparam[0],
                                        totpts,
					ZERO);
     t1.vthinc0 = t4.tmp1;
     t1.vthinc0 += VINC0;
     t1.vthinc1 = t4.tmp1;
     t1.vthinc1 += VINC1;
     gms::common::avx512vec16_init_unroll8x(&m_gsc.yparam[0],
                                        totpts,
					ZERO);
     t1.vthinc2 = t4.tmp1;
     t1.vthinc2 += VINC2;
     t1.vthinc3 = t4.tmp1;
     t1.vthinc3 += VINC3;
     gms::common::avx512vec16_init_unroll8x(&m_gsc.zparam[0],
                                        totpts,
					ZERO);
     // Loop over grass cylinders
#if defined __ICC || defined __INTEL_COMPILER
    
     __assume_aligned(m_gsc.A,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    
     m_gsc.A = (float*)__builtin_assume_aligned(m_gsc.A,64);
#endif
#if (SAMPLE_HW_PMC) == 1
     uint64_t cycles_lvl2 = 0xFFFFFFFFFFFFFFFFULL;
    #if (CHECK_L2_LICENSE) == 1
         const int32_t niter = 128;
	 cycles_lvl2 = rdpmc(0); // assume first PMC is 0
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
     uint64_t prog_counters_start[5] = {};
     uint64_t prog_counters_end[5]   = {};
     uint64_t tsc_start,tsc_stop;
     uint64_t act_cyc_start,act_cyc_stop;
     uint64_t ref_cyc_start,ref_cyc_stop;
     uint64_t volatile dummy1,dummy2,dummy3,dummy4,dummy5,dummy6;
     double utilization,nom_ghz,avg_ghz;
     // Forcing the IF and decoding of PMC functions  before the measurement loop
     dummy1 = rdpmc_actual_cycles();
     dummy2 = rdpmc_reference_cycles();
     dummy3 = rdtscp();
     dummy4 = rdpmc(0);
     // Actual measurements
     // CPUID to force serialization.
     dummy5 = get_core_counter_width();
     if(cycles_lvl2) {
        for(int32_t i = 1; i != 5; ++i) {
            prog_counters_start[i] = rdpmc(i);
        }
     }
     else {
         for(int32_t i = 0; i != 4; ++i) {
            prog_counters_start[i] = rdpmc(i);
        }
     }
     act_cyc_start = rdpmc_actual_cycles();
     ref_cyc_start = rdpmc_reference_cycles();
     tsc_start = rdtscp();
#endif
      for(int32_t i = 0; i != m_gsc.nplants; ++i) {
         seedr = std::clock();
	 auto rand_r = std::bind(std::uniform_real_distribution<float>(0.5f,0.9f),
	                         std::mt19937(seedr));
	 const float rtemp = rand_r();
	 //
	 t5.tvrad = rScale * AVX512Vec16{rtemp};
	 seedz = std::clock();
	 auto rand_z = std::bind(std::uniform_real_distribution<float>(0.3f,0.5f),
	                         std::mt19937(seedz));
	 const float ztemp = rand_z();
	 const float c0 = n2PI*rtemp;
	 //
	 // Compute cylinders surface.
	 m_gsc.A[i] = c0*ztemp+c0*rtemp;
	 t5.tvz   = hScale * AVX512Vec16{ztemp};
	 t4.tmp2  = t4.vNPTS;
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
#if defined __ICC || defined __INTEL_COMPILER
	 __assume_aligned(m_gsc.xparam,64);
	 __assume_aligned(m_gsc.yparam,64);
	 __assume_aligned(m_gsc.zparam,64);
#pragma vector always
#pragma vector vectorlength(16)
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        m_gsc.xparam = (AVX512Vec16*)__builtin_assume_aligned(m_gsc.xparam,64);
	m_gsc.yparam = (AVX512Vec16*)__builtin_assume_aligned(m_gsc.yparam,64);
	m_gsc.zparam = (AVX512Vec16*)__builtin_assume_aligned(m_gsc.zparam,64);
#pragma omp simd aligned(m_gsc.xparam,m_gsc.yparam,m_gsc.zparam:64)
#endif
        for(int32_t j = 0; j != m_gsc.param_npts-3; j += 4) {
             // for every point -- do...
	    t0.vtheta0 = t0.vtheta0+t1.vthinc0;
	    m_gsc.xparam[Ix2D(i,m_gsc.param_npts,j+0)] = t5.tvrad*cos(t0.vtheta0);
	    m_gsc.yparam[Ix2D(i,m_gsc.param_npts,j+0)] = t5.tvrad*sin(t0.vtheta0);
	    t3.vhinit0 = t3.vhinit0+t2.vhinc0;
	    m_gsc.zparam[Ix2D(i,m_gsc.param_npts,j+0)] = t3.vhinit0;
	    t0.vtheta1 = t0.vtheta1+t1.vthinc1;
	    m_gsc.xparam[Ix2D(i,m_gsc.param_npts,j+1)] = t5.tvrad*cos(t0.vtheta1);
	    m_gsc.yparam[Ix2D(i,m_gsc.param_npts,j+1)] = t5.tvrad*sin(t0.vtheta1);
	    t3.vhinit1 = t3.vhinit1+t2.vhinc1;
	    m_gsc.zparam[Ix2D(i,m_gsc.param_npts,j+1)] = t3.vhinit1;
	    t0.vtheta2 = t0.vtheta2+t1.vthinc2;
	    m_gsc.xparam[Ix2D(i,m_gsc.param_npts,j+2)] = t5.tvrad*cos(t0.vtheta2);
	    m_gsc.yparam[Ix2D(i,m_gsc.param_npts,j+2)] = t5.tvrad*sin(t0.vtheta2);
	    t3.vhinit2 = t3.vhinit2+t2.vhinc2;
	    m_gsc.zparam[Ix2D(i,m_gsc.param_npts,j+2)] = t3.vhinit2;
	    t0.vtheta3 = t0.vtheta3+t1.vthinc3;
	    m_gsc.xparam[Ix2D(i,m_gsc.param_npts,j+3)] = t5.tvrad*cos(t0.vtheta3);
	    m_gsc.yparam[Ix2D(i,m_gsc.param_npts,j+3)] = t5.tvrad*sin(t0.vtheta3);
	    t3.vhinit3 = t3.vhinit3+t2.vhinc3;
	    m_gsc.zparam[Ix2D(i,m_gsc.param_npts,j+3)] = t3.vhinit3;
	}
     }
#if (SAMPLE_HW_PMC) == 1
            dummy6 = get_core_counter_width()
	    if(cycles_lvl2) {
	       for(int32_t i = 1; i != 5; ++i) {
                   prog_counters_stop[i] = rdpmc(i);
               }
	    }
	    else {
               for(int32_t i = 0; i != 4; ++i) {
                   prog_counters_stop[i] = rdpmc(i);
               }
	    }
            act_cyc_stop = rdpmc_actual_cycles();
            ref_cyc_stop = rdpmc_reference_cycles();
            tsc_stop = rdtscp();
	    nom_ghz = get_TSC_frequency()/1.0e9;
	    utilization = (double)(ref_cyc_end-ref_cyc_start)/(double)(tsc_end-tsc_start);
	    avg_ghz = (double)(act_cyc_end-act_cyc_start)/(double)(tsc_end-tsc_start)*nom_ghz;
	    syslog(LOG_INFO,"%-10s:\n", __PRETTY_FUNCTION__);
	    syslog(LOG_INFO, "*************** Hardware Counters -- Dump Begin **************");
	    syslog(LOG_INFO,"Core utilization                      : %f\n",utilization );
	    syslog(LOG_INFO,"Core average frequency                : %f\n",avg_ghz);
	    syslog(LOG_INFO,"Reference cycles                      : %20lld\n",ref_cyc_end-ref_cyc_start);
	    syslog(LOG_INFO,"Actual cycles                         : %20lld\n",act_cyc_stop-act_cyc_start);
	    if(!cycles_lvl2) {
	       syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event1                    , prog_counters_stop[0]-prog_counters_start[0]);
	       syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event2                    , prog_counters_stop[1]-prog_counters_start[1]);
	       syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event3                    , prog_counters_stop[2]-prog_counters_start[2]);
	       syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event4                    , prog_counters_stop[3]-prog_counters_start[3]);
	       syslog(LOG_INFO, "*************** Hardware Counters -- Dump End   **************");
	    }
	     else {
                 syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event1                    , prog_counters_stop[1]-prog_counters_start[1]);
	         syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event2                    , prog_counters_stop[2]-prog_counters_start[2]);
	         syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event3                    , prog_counters_stop[3]-prog_counters_start[3]);
	         syslog(LOG_INFO,"%-37s: %20lld\n", pmc_event4                    , prog_counters_stop[4]-prog_counters_start[4]);
	         syslog(LOG_INFO, "*************** Hardware Counters -- Dump End   **************");
	     }
#endif
}

#include <math.h>
void
gms::math::GrassScattererAVX512::
ComputeGrassHVPolarization(const float gamma,
                           const float ah,
			   const float av) {
        const static float n2PI = 6.283185307179586f;
     const static float n1rad = 0.0174444444444444444f;
     std::complex<float> dielectric = {0.0f,0.0f};
     std::complex<float> ch0        = {0.0f,0.0f};
     //
     std::complex<float> cv0        = {0.0f,0.0f};
     //
     const float n28PI = 14.0f*TWOPI;
     const float K = n2PI/gamma;
     const float sarea = 0.0f;
     float term1 = 0.0f;
     float term2 = 0.0f;
     float term3 = 0.0;
     float term4 = 0.0f;
     float term5 = 0.0f;
     float t     = 0.0f;
     float re    = 0.0f;
     float im    = 0.0f;
     float c0    = 0.0f;
     float c1    = 0.0f;
#if defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(m_gsc.A,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     m_gsc.A = (float*)__builtin_assume_aligned(m_gsc.A,64);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)

#endif
     for(int32_t i = 0; i != m_gsc.nplants-3; i += 4) {
     
         float t0 += m_gsc.A[i+0];
	 float t1 += m_gsc.A[i+1];
	 float t2 += m_gsc.A[i+2];
	 float t3 += m_gsc.A[i+3]
	// const float tot = t0+t1;
	// m_gsc.tot_area += tot;
     }
     m_gsc.tot_area = t0+t1+t2+t3;
     sarea = m_gsc.tot_area*m_gsc.tot_area;
     term1 = static_cast<float>(m_gsc.nplants)*asqrt*K*K;
     c0 = m_gsc.epsilon.real();
     t = 1.0f/(1.0f+c0);
     c1 = m_gsc.epsilon.imag();
     re = c0-1.0f*c0-1.0f;
     im = c1*c1;
     dielectric = {c0,c1};
     term2 = 3.0f+16.0f*t+96.0f*t*t;
     term3 = 3.0f*((ah/K)*(ah/K));
     term4 = 12+8.0f*t-64.0f*t*t;
     term5 = 3.0f*((av/K)*(av/K));
     for(int32_t i = 0; i != 90; ++i) {
         // Full grazing angle sweep.
	 const float theta = n1rad*static_cast<float>(i);
	 const float t0    = sin(theta);
	 const float a0    = 4.0f*(1.0f+2.0f*t0*t0);
	 ch0  = dielectric*term2;
	 const float t1    = term1/(n28PI*t0);
	 const float t2    = term3+a0;
	 m_gsh.Polh[i]     = t1*ch0/t2;
	 const float t3    = cos(theta);
	 const float t4    = term2+t3*t3*term4;
	 cv1 = dielectric*t4;
	 const float t5    = term5+a0;
	 m_gsh.Polv[i]     = t1*cv1/t5;
     }
     
}


