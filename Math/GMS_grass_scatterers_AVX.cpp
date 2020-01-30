#if defined __GNUC__ && !defined __INTEL_COMPILER
#include <omp.h>
#endif


#if defined __GNUC__ || defined __INTEL_COMPILER
#include <random>
#endif
#include <ctime>
#include <cstdlib>
#endif
#include "GMS_grass_scatterers_AVX.h"
#include "GMS_malloc.h"
#include "GMS_indices.h"
#include "GMS_common.h"
#include "GMS_error_macros.h"

gms::math::
GrassScattererAVX::GrassScattererAVX() {

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

gms::math::GrassScattererAVX::
GrassScattererAVX(const int32_t nsteps,
                  const int32_t ordinal,
		  const int32_t param_npts,
		  const float lat,
		  const float lon,
		  const float elev,
		  const std::complex<float> cepsilon) {

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
  if(rand < lo) rand = 50;
  if(rand > hi) rand = 500;
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
  m_gsc.A         = gms_efmalloca(static_cast<size_t>(m_gsc.nplants),64);
  m_gsc.moistness = gms_eimalloca4(static_cast<size_t>(m_gsc.nplants),64);
  m_gsc.xparam    = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.param_npts),64);
  m_gsc.yparam    = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.param_npts),64);
  m_gsc.zparam    = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.param_npts),64);
  m_gsh.Polv[96]  = {0.0f,0.0f};
  m_gsh.Polh[96]  = {0.0f,0.0f};
  m_gsh.xang      = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),64);
  m_gsh.sin_xang  = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),64);
  m_gsh.cos_xang  = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),64);
  m_gsh.yang      = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),64);
  m_gsh.sin_yang  = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),64);
  m_gsh.cos_yang  = gms_avxvec8_emalloca(static_cast<size_t>(m_gsc.nsteps*m_gsc.nplants),64);
}

gms::math::GrassScattererAVX
::~GrassScattererAVX() {

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
}

void
gms::math::GrassScattererAVX::
SetGrassMoistnessMask() {

  std::random_device rd;
  // Memory first touch
  for(int32_t i = 0; i != m_gsc.nplants; ++i) {
      m_gsc.moistness[i] = 0U;
  }
  std::mt19937 rgen(rd());
  std::uniform_int_distribution<> distro(0,1);
  for(int32_t i = 0; i != m_gsc.nplants; ++i) {
      m_gsc.moistness[i] = distro(rgen);
  }
}

void
gms::math::GrassScattererAVX::
ComputeGrassParamEq_ymm8r4() {
                           
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

     struct _T4_ {
       AVXVec8 vNPTS;
       //
       AVXVec8 tmp1;
       AVXVec8 tmp2;
     } __ATTR_ALIGN__(64) t4;

     struct _T5_ {
       AVXVec8 tvrad;
       AVXVec8 tvz;
     } __ATTR_ALIGN__(64) t5;
     AVXVec8 rScale = AVXVec8{10.0f}; // unit of mm
     AVXVec8 hScale = AVXVec8{100.0f}; // unit of cm
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
     t4.vNPTS   = AVXVec8{static_cast<float>(m_tsc.param_npts)};
     
     t4.tmp1    = PI/t4.vNPTS;
     t5.tvrad   = ZERO;
     t5.tvz     = ZERO;
     // Memory first touch
     gms::common::avxvec8_init_unroll8x(&m_gsc.xparam[0],
                                        totpts,
					ZERO);
     t1.vthinc0 = t4.tmp1;
     t1.vthinc0 += VINC0;
     t1.vthinc1 = t4.tmp1;
     t1.vthinc1 += VINC1;
     gms::common::avxvec8_init_unroll8x(&m_gsc.yparam[0],
                                        totpts,
					ZERO);
     t1.vthinc2 = t4.tmp1;
     t1.vthinc2 += VINC2;
     t1.vthinc3 = t4.tmp1;
     t1.vthinc3 += VINC3;
     gms::common::avxvec8_init_unroll8x(&m_gsc.zparam[0],
                                        totpts,
					ZERO);
     // Loop over grass cylinders
#if defined __ICC || defined __INTEL_COMPILER
    
     __assume_aligned(m_gsc.A,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    
     m_gsc.A = (float*)__builtin_assume_aligned(m_gsc.A,64);
#endif
     for(int32_t i = 0; i != m_gsc.nplants; ++i) {
         seedr = std::clock();
	 auto rand_r = std::bind(std::uniform_real_distribution<float>(0.5f,0.9f),
	                         std::mt19937(seedr));
	 const float rtemp = rand_r();
	 //
	 t5.tvrad = rScale * AVXVec8{rtemp};
	 seedz = std::clock();
	 auto rand_z = std::bind(std::uniform_real_distribution<float>(0.3f,0.5f),
	                         std::mt19937(seedz));
	 const float ztemp = rand_z();
	 const float c0 = n2PI*rtemp;
	 //
	 // Compute cylinders surface.
	 m_gsc.A[i] = c0*ztemp+c0*rtemp;
	 t5.tvz   = hScale * AVXVec8{ztemp};
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
#pragma vector vectorlength(8)
#pragma code_align(32)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        m_gsc.xparam = (AVXVec8*)__builtin_assume_aligned(m_gsc.xparam,64);
	m_gsc.yparam = (AVXVec8*)__builtin_assume_aligned(m_gsc.yparam,64);
	m_gsc.zparam = (AVXVec8*)__builtin_assume_aligned(m_gsc.zparam,64);
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
}

#include <math.h>
void
gms::math::GrassScattererAVX::
ComputeGrassHVPolarization(
			   const float gamma,
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



