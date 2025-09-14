
#ifndef __GMS_SLEEFSIMDSP_H__
#define __GMS_SLEEFSIMDSP_H__


//          Copyright Naoki Shibata 2010 - 2017.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Always use -ffp-contract=off option to compile SLEEF.

#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "misc.h"
#include "GMS_config.h"

void Sleef_x86CpuID(int32_t out[4], uint32_t eax, uint32_t ecx);

//#if (defined(_MSC_VER))
//#pragma fp_contract (off)
//#endif

//#ifdef ENABLE_SSE2
//#define CONFIG 2
//#include "helpersse2.h"
//#ifdef DORENAME
//#include "renamesse2.h"
//#endif
//#endif

//#ifdef ENABLE_AVX
//#define CONFIG 1
//#include "helperavx.h"
//#//ifdef DORENAME
//#include "renameavx.h"
//#endif
//#endif

/*#ifdef ENABLE_FMA4
#define CONFIG 4
#include "helperavx.h"
#ifdef DORENAME
#include "renamefma4.h"
#endif
#endif*/

#ifdef ENABLE_AVX2
#define CONFIG 1
#include "GMS_helperavx2.h"
//#ifdef DORENAME
//#include "renameavx2.h"
//#endif
#endif

#ifdef ENABLE_AVX512F
#define CONFIG 1
#include "GMS_helperavx512f.h"
//#ifdef DORENAME
//#include "renameavx512f.h"
//#endif
#endif

/*#ifdef ENABLE_VECEXT
#define CONFIG 1
#include "helpervecext.h"
#ifdef DORENAME
#include "renamevecext.h"
#endif
#endif*/

/*#ifdef ENABLE_PUREC
#define CONFIG 1
#include "helperpurec.h"
#ifdef DORENAME
#include "renamepurec.h"
#endif
#endif*/

/*#ifdef ENABLE_NEON32
#define CONFIG 1
#include "helperneon32.h"
#ifdef DORENAME
#include "renameneon32.h"
#endif
#endif*/

//

#include "GMS_df.h"

//

#define PI4_Af 0.78515625f
#define PI4_Bf 0.00024187564849853515625f
#define PI4_Cf 3.7747668102383613586e-08f
#define PI4_Df 1.2816720341285448015e-12f

#define PI_Af 3.140625f
#define PI_Bf 0.0009670257568359375f
#define PI_Cf 6.2771141529083251953e-07f
#define PI_Df 1.2154201256553420762e-10f

#define PI_XDf 1.2141754268668591976e-10f
#define PI_XEf 1.2446743939339977025e-13f

#define TRIGRANGEMAXf 1e+7 // 39000
#define SQRT_FLT_MAX 18446743523953729536.0

#define L2Uf 0.693145751953125f
#define L2Lf 1.428606765330187045e-06f
#define R_LN2f 1.442695040888963407359924681001892137426645954152985934135449406931f

//

static INLINE vopmask visnegzero_vo_vf(vfloat d) {
  return veq_vo_vi2_vi2(vreinterpret_vi2_vf(d), vreinterpret_vi2_vf(vcast_vf_f(-0.0)));
}

static INLINE vmask vsignbit_vm_vf(vfloat f) {
  return vand_vm_vm_vm(vreinterpret_vm_vf(f), vreinterpret_vm_vf(vcast_vf_f(-0.0f)));
}

static INLINE vfloat vmulsign_vf_vf_vf(vfloat x, vfloat y) {
  return vreinterpret_vf_vm(vxor_vm_vm_vm(vreinterpret_vm_vf(x), vsignbit_vm_vf(y)));
}

static INLINE vfloat vsign_vf_vf(vfloat f) {
  return vreinterpret_vf_vm(vor_vm_vm_vm(vreinterpret_vm_vf(vcast_vf_f(1.0f)), vand_vm_vm_vm(vreinterpret_vm_vf(vcast_vf_f(-0.0f)), vreinterpret_vm_vf(f))));
}

static INLINE vopmask vsignbit_vo_vf(vfloat d) {
  return veq_vo_vi2_vi2(vand_vi2_vi2_vi2(vreinterpret_vi2_vf(d), vcast_vi2_i(0x80000000)), vcast_vi2_i(0x80000000));
}

static INLINE vint2 vsel_vi2_vf_vf_vi2_vi2(vfloat f0, vfloat f1, vint2 x, vint2 y) {
  return vsel_vi2_vo_vi2_vi2(vlt_vo_vf_vf(f0, f1), x, y);
}

static INLINE vint2 vsel_vi2_vf_vi2(vfloat d, vint2 x) {
  return vand_vi2_vo_vi2(vsignbit_vo_vf(d), x);
}

#ifndef ENABLE_AVX512F
static INLINE vint2 vilogbk_vi2_vf(vfloat d) {
  vopmask o = vlt_vo_vf_vf(d, vcast_vf_f(5.421010862427522E-20f));
  d = vsel_vf_vo_vf_vf(o, vmul_vf_vf_vf(vcast_vf_f(1.8446744073709552E19f), d), d);
  vint2 q = vand_vi2_vi2_vi2(vsrl_vi2_vi2_i(vcast_vi2_vm(vreinterpret_vm_vf(d)), 23), vcast_vi2_i(0xff));
  q = vsub_vi2_vi2_vi2(q, vsel_vi2_vo_vi2_vi2(o, vcast_vi2_i(64 + 0x7f), vcast_vi2_i(0x7f)));
  return q;
}
#endif

//

/*EXPORT*/
  
  __ATTR_ALWAYS_INLINE__
  static
  inline
  vint2 xilogbf(vfloat d) {
  vint2 e = vilogbk_vi2_vf(vabs_vf_vf(d));
  e = vsel_vi2_vo_vi2_vi2(veq_vo_vf_vf(d, vcast_vf_f(0.0f)), vcast_vi2_i(FP_ILOGB0), e);
  e = vsel_vi2_vo_vi2_vi2(visnan_vo_vf(d), vcast_vi2_i(FP_ILOGBNAN), e);
  e = vsel_vi2_vo_vi2_vi2(visinf_vo_vf(d), vcast_vi2_i(INT_MAX), e);
  return e;
}

static inline vfloat vpow2i_vf_vi2(vint2 q) {
  return vreinterpret_vf_vm(vcast_vm_vi2(vsll_vi2_vi2_i(vadd_vi2_vi2_vi2(q, vcast_vi2_i(0x7f)), 23)));
}

  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat vldexp_vf_vf_vi2(vfloat x, vint2 q); 

/*EXPORT*/ 
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat 
  xldexpf(vfloat x, vint2 q);

/*EXPORT*/ 
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
 
  vfloat xsinf(vfloat d); 

 /*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
 
  vfloat xcosf(vfloat d); 
  

  /*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat2 xsincosf(vfloat d); 
  
  
/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xtanf(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xsinf_u1(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xcosf_u1(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat2 xsincosf_u1(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat2 xsincospif_u05(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat2 xsincospif_u35(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xtanf_u1(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xatanf(vfloat d); 
  
  
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat atan2kf(vfloat y, vfloat x); 
  
  
  
  __ATTR_ALWAYS_INLINE__
  static
  inline
  vfloat visinf2_vf_vf_vf(vfloat d, vfloat m) {
  return vreinterpret_vf_vm(vand_vm_vo32_vm(visinf_vo_vf(d), vor_vm_vm_vm(vsignbit_vm_vf(d), vreinterpret_vm_vf(m))));
}

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xatan2f(vfloat y, vfloat x); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xasinf(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xacosf(vfloat d); 
  

//

  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat2 atan2kf_u1(vfloat2 y, vfloat2 x);
  
/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xatan2f_u1(vfloat y, vfloat x); 
  
  
/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xasinf_u1(vfloat d); 
  

 /*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xacosf_u1(vfloat d); 
  

/*EXPORT*/
  
  __ATTR_ALWAYS_INLINE__
  static
  inline
  vfloat xatanf_u1(vfloat d) {
  vfloat2 d2 = atan2kf_u1(vcast_vf2_vf_vf(vabs_vf_vf(d), vcast_vf_f(0)), vcast_vf2_f_f(1, 0));
  vfloat r = vadd_vf_vf_vf(d2.x, d2.y);
  r = vsel_vf_vo_vf_vf(visinf_vo_vf(d), vcast_vf_f(1.570796326794896557998982), r);
  return vmulsign_vf_vf_vf(r, d);
}

//

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xlogf(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xexpf(vfloat d); 

/*#ifdef ENABLE_NEON32
EXPORT vfloat xsqrtf(vfloat d) {
  vfloat e = vreinterpret_vf_vi2(vadd_vi2_vi2_vi2(vcast_vi2_i(0x20000000), vand_vi2_vi2_vi2(vcast_vi2_i(0x7f000000), vsrl_vi2_vi2_i(vreinterpret_vi2_vf(d), 1))));
  vfloat m = vreinterpret_vf_vi2(vadd_vi2_vi2_vi2(vcast_vi2_i(0x3f000000), vand_vi2_vi2_vi2(vcast_vi2_i(0x01ffffff), vreinterpret_vi2_vf(d))));
  float32x4_t x = vrsqrteq_f32(m);
  x = vmulq_f32(x, vrsqrtsq_f32(m, vmulq_f32(x, x)));
  float32x4_t u = vmulq_f32(x, m);
  u = vmlaq_f32(u, vmlsq_f32(m, u, u), vmulq_f32(x, vdupq_n_f32(0.5)));
  e = vreinterpret_vf_vm(vandnot_vm_vo32_vm(veq_vo_vf_vf(d, vcast_vf_f(0)), vreinterpret_vm_vf(e)));
  u = vmul_vf_vf_vf(e, u);

  u = vsel_vf_vo_vf_vf(visinf_vo_vf(d), vcast_vf_f(INFINITYf), u);
  u = vreinterpret_vf_vm(vor_vm_vo32_vm(vor_vo_vo_vo(visnan_vo_vf(d), vlt_vo_vf_vf(d, vcast_vf_f(0))), vreinterpret_vm_vf(u)));
  u = vmulsign_vf_vf_vf(u, d);

  return u;
}
#elif defined(ENABLE_CLANGVEC)
EXPORT vfloat xsqrtf(vfloat d) {
  vfloat q = vsqrt_vf_vf(d);
  q = vsel_vf_vo_vf_vf(visnegzero_vo_vf(d), vcast_vf_f(-0.0), q);
  return vsel_vf_vo_vf_vf(vispinf_vo_vf(d), INFINITYf, q);
}
#else*/
/*EXPORT*/
  
  __ATTR_ALWAYS_INLINE__
  static
  inline
vfloat xsqrtf(vfloat d) { return vsqrt_vf_vf(d); }
//#endif

/*EXPORT*/

  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xcbrtf(vfloat d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xcbrtf_u1(vfloat d); 
  

  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
   vfloat2 logkf(vfloat d); 
   

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xlogf_u1(vfloat d); 

  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat expkf(vfloat2 d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
 
  vfloat xpowf(vfloat x, vfloat y); 
  

  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat2 expk2f(vfloat2 d); 
  
/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xsinhf(vfloat x); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xcoshf(vfloat x); 
  

 /*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xtanhf(vfloat x); 
  
  
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat2 logk2f(vfloat2 d); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xasinhf(vfloat x); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xacoshf(vfloat x); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xatanhf(vfloat x); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xexp2f(vfloat a); 
  
/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xexp10f(vfloat a); 
  

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xexpm1f(vfloat a); 

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
  
  vfloat xlog10f(vfloat a); 

/*EXPORT*/
  __ATTR_VECTORCALL__
  __ATTR_ALIGN__(32)
 
  vfloat xlog1pf(vfloat a); 























#endif /*__GMS_SLEEFSIMDSP_H__*/
