
#ifndef __GMS_SIMD_MACROS_H__
#define __GMS_SIMD_MACROS_H__

namespace file_info {
#if defined _WIN64  
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif
	const unsigned int gGMS_SIMD_MACROS_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_SIMD_MACROS_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_SIMD_MACROS_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_SIMD_MACROS_FULLVER = 
	  1000U*gGMS_SIMD_MACROS_MAJOR+100U*gGMS_SIMD_MACROS_MINOR+10U*gGMS_SIMD_MACROS_MICRO;

	const char * const pgGMS_SIMD_MACROS_CREATE_DATE = "01-10-2018 10:41 + 00200 (MON 01 OCT 2018 GMT + 2)";

	const char * const pgGMS_SIMD_MACROS_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_SIMD_MACROS_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_SIMD_MACROS_SYNOPSIS = "#defines of AVX and AVX512 intrinsics.";
}



#include <immintrin.h>


namespace gms {
	namespace math {

#define Vec4_SETZERO  _mm256_setzero_pd()
#define Vec4_SET1(x)  _mm256_set1_pd((x))
#define Vec4_SET(x0,x1,x2,x3) _mm256_set_pd((x0),(x1),(x2),(x3))
#define Vec4_SETI4(x)    _mm_set1_epi32((x))
#define Vec4_SQRT(x)  _mm256_sqrt_pd((x))
#define Vec4_POW(x,y)   _mm256_pow_pd((x),(y))
#define Vec4_MUL(x,y)   _mm256_mul_pd((x),(y))
#define Vec4_CMP(x,y,p) _mm256_cmp_pd((x),(y),(p))
#define Vec4_SUB(x,y)   _mm256_sub_pd((x),(y))
#define Vec4_ADD(x,y)   _mm256_add_pd((x),(y))
#define Vec4_TESTZ(x,y) _mm256_testz_pd((x),(y))
#define Vec4_CVTI4(x)   _mm256_cvtepi32_pd((x))
#define Vec4_DIV(x,y)   _mm256_div_pd((x),(y))
#define Vec4_AND(x,y)   _mm256_and_pd((x),(y))    
#define Vec4_SIN(x)     _mm256_sin_pd((x))
#define Vec4_COS(x)     _mm256_cos_pd((x))
#define Vec4_EXP(x)     _mm256_exp_pd((x))
#define Vec4_LOAD(x)    _mm256_load_pd((x))
#define Vec4_LOADU(x)   _mm256_loadu_pd((x))
#define Vec4_LOG(x)     _mm256_log_pd((x))
#define Vec4_FMAD(x,y,z) _mm256_fmadd_pd((x),(y),(z))
#define Vec4_FMSUB(x,y,z) _mm256_fmusb((x),(y),(z))
#define Vec4_ATAN(x)      _mm256_atan_pd((x))
#define Vec4_TAN(x)       _mm256_tan_pd((x))

#define Vec8_SETZERO  _mm256_set1_pd(0.0)
#define Vec8_SET1(x)  _mm512_set1_pd((x))
#define Vec8_SET(x0,x1,x2,x3,x4,x5,x6,x7) _mm512_set_pd((x0),(x1),(x2),(x3),(x4),(x5),(x6),(x7))
#define Vec8_SQRT(x)  _mm512_sqrt_pd((x))
#define Vec8_POW(x,y) _mm512_pow_pd((x),(y))
#define Vec8_MUL(x,y) _mm512_mul_pd((x),(y))
#define Vec8_CMP(x,y,imm) _mm512_cmp_pd_mask((x),(y),(imm))
#define Vec8_SUB(x,y)    _mm512_sub_pd((x),(y))
#define Vec8_ADD(x,y)    _mm512_sub_pd((x),(y))
#define Vec8_DIV(x,y)    _mm512_div_pd((x),(y))
#define Vec8_AND(x,y)    _mm512_and_pd((x),(y))
#define Vec8_SIN(x)      _mm512_sin_pd((x))
#define Vec8_SINH(x)     _mm512_sinh_pd((x))
#define Vec8_COS(x)      _mm512_cos_pd((x))
#define Vec8_COSH(x)     _mm512_cosh_pd((x))
#define Vec8_EXP(x)      _mm512_exp_pd((x))
#define Vec8_SETI4(x)    _mm256_set1_epi32((x))
#define Vec8_CVTI4(x)    _mm512_cvtepi32_pd((x))
#define Vec8_ABS(x)      _mm512_abs_pd((x))
#define Vec8_BLEND(imm,x,y) _mm512_mask_blend_pd((imm),(x),(y));
#define Vec8_LOAD(x)      _mm512_load_pd((x))
#define Vec8_LOADU(x)     _mm512_loadu_pd((x))
#define Vec8_LOG(x)       _mm512_log_pd((x))
#define Vec8_FMAD(x,y,z)  _mm512_fmadd_pd((x),(y),(z))
#define Vec8_FMSUB(x,y,z) _mm512_fmsub_pd((x),(y),(z))
#define Vec8_TAN(x)       _mm512_tan_pd((x))
#define Vec8_ATAN(x)      _mm512_atan((x))

	}
}


#endif /*__GMS_SIMD_MACROS_H__*/
