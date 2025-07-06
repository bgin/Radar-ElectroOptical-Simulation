



/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/




#include <immintrin.h>
#include "GMS_csinhv_svml_zmm16r4.h"
//#include "GMS_cephes.h"
#include <cmath>


__ATTR_ALWAYS_INLINE__
static inline
float ceph_sinf( const float xx ) {
constexpr float FOPI = 1.27323954473516;
/* Note, these constants are for a 32-bit significand: */
/*
static float DP1 =  0.7853851318359375;
static float DP2 =  1.30315311253070831298828125e-5;
static float DP3 =  3.03855025325309630e-11;
static float lossth = 65536.;
*/

/* These are for a 24-bit significand: */
constexpr float PIO4F =  0.7853981633974483096f;
constexpr float DP1 = 0.78515625;
constexpr float DP2 = 2.4187564849853515625e-4;
constexpr float DP3 = 3.77489497744594108e-8;
constexpr float lossth = 8192.;
constexpr float T24M1 = 16777215.;

float sincof[] = {
-1.9515295891E-4,
 8.3321608736E-3,
-1.6666654611E-1
};
float coscof[] = {
 2.443315711809948E-005,
-1.388731625493765E-003,
 4.166664568298827E-002
};
float *p;
float x, y, z;
register unsigned long j;
register int sign;
sign = 1;
x = xx;
if( xx < 0 )
	{
	sign = -1;
	x = -xx;
	}
if( x > T24M1 )
	{
	//mtherr( "sinf", TLOSS );
	return(0.0);
	}
j = FOPI * x; /* integer part of x/(PI/4) */
y = j;
/* map zeros to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}
j &= 7; /* octant modulo 360 degrees */
/* reflect in x axis */
if( j > 3)
	{
	sign = -sign;
	j -= 4;
	}

if( x > lossth )
	{
	//mtherr( "sinf", PLOSS );
	x = x - y * PIO4F;
	}
else
	{
/* Extended precision modular arithmetic */
	x = ((x - y * DP1) - y * DP2) - y * DP3;
	}
/*einits();*/
z = x * x;
if( (j==1) || (j==2) )
	{
/* measured relative error in +/- pi/4 is 7.8e-8 */
/*
	y = ((  2.443315711809948E-005 * z
	  - 1.388731625493765E-003) * z
	  + 4.166664568298827E-002) * z * z;
*/
	p = coscof;
	y = *p++;
	y = y * z + *p++;
	y = y * z + *p++;
	y *= z * z;
	y -= 0.5 * z;
	y += 1.0;
	}
else
	{
/* Theoretical relative error = 3.8e-9 in [-pi/4, +pi/4] */
/*
	y = ((-1.9515295891E-4 * z
	     + 8.3321608736E-3) * z
	     - 1.6666654611E-1) * z * x;
	y += x;
*/
	p = sincof;
	y = *p++;
	y = y * z + *p++;
	y = y * z + *p++;
	y *= z * x;
	y += x;
	}
/*einitd();*/
if(sign < 0)
	y = -y;
return( y);
}

__ATTR_ALWAYS_INLINE__
static inline
float ceph_cosf( const float xx ) {
/* Note, these constants are for a 32-bit significand: */
/*
static float DP1 =  0.7853851318359375;
static float DP2 =  1.30315311253070831298828125e-5;
static float DP3 =  3.03855025325309630e-11;
static float lossth = 65536.;
*/

/* These are for a 24-bit significand: */
constexpr float PIO4F =  0.7853981633974483096f;
constexpr float FOPI = 1.27323954473516f;
constexpr float DP1 = 0.78515625f;
constexpr float DP2 = 2.4187564849853515625e-4f;
constexpr float DP3 = 3.77489497744594108e-8f;
constexpr float lossth = 8192.f;
constexpr float T24M1 = 16777215.f;
float x, y, z;
int j, sign;

/* make argument positive */
sign = 1;
x = xx;
if( x < 0 )
	x = -x;

if( x > T24M1 )
	{
	//mtherr( "cosf", TLOSS );
	return(0.0);
	}

j = FOPI * x; /* integer part of x/PIO4 */
y = j;
/* integer and fractional part modulo one octant */
if( j & 1 )	/* map zeros to origin */
	{
	j += 1;
	y += 1.0;
	}
j &= 7;
if( j > 3)
	{
	j -=4;
	sign = -sign;
	}

if( j > 1 )
	sign = -sign;

if( x > lossth )
	{
	//mtherr( "cosf", PLOSS );
	x = x - y * PIO4F;
	}
else
/* Extended precision modular arithmetic */
	x = ((x - y * DP1) - y * DP2) - y * DP3;

z = x * x;

if( (j==1) || (j==2) )
	{
	y = (((-1.9515295891E-4 * z
	     + 8.3321608736E-3) * z
	     - 1.6666654611E-1) * z * x)
	     + x;
	}
else
	{
	y = ((  2.443315711809948E-005 * z
	  - 1.388731625493765E-003) * z
	  + 4.166664568298827E-002) * z * z;
	y -= 0.5 * z;
	y += 1.0;
	}
if(sign < 0)
	y = -y;
return( y );
}





void gms::math::csinhv_svml_zmm16r4_u10x_u(const float * __restrict xre,
                                           const float * __restrict xim,
                                           float * __restrict csre,
                                           float * __restrict csim,
                                           int32_t n) 
{

                       if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                           __m512 zmm24,zmm25,zmm26,zmm27;
                           __m512 zmm28,zmm29,zmm30,zmm31;
                          int32_t i;
#if (CSINHV_SVML_ZMM16R4_PEEL_LOOP) == 1
                            while(((uintptr_t)&csre & 63  &&
                                   (uintptr_t)&csim & 63) && n) 
                            {
                                    const float z0 = *xre;
                                    const float z1 = *xim;
                                    const float z2 = std::sinh(z0)*ceph_cosf(z1);
                                    *csre          = z2;
                                    const float z3 = std::cosh(z0)*ceph_sinf(z1);
                                    *csim          = z3;
                                    xre++;
                                    xim++;
                                    csre++;
                                    csim++;
                                    n--;
                            }
#endif 
 
                          for(i = 0; (i+159) < n; i += 160) 
                          {
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
#endif 
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+16],_MM_HINT_T0);
#endif 
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
#endif                               
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
#endif 
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
#endif                             
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
#endif 
                              zmm20 = _mm512_loadu_ps(&xre[i+80]);
                              zmm21 = _mm512_loadu_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_storeu_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_storeu_ps(&csim[i+80],zmm23);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
#endif                             
                              zmm24 = _mm512_loadu_ps(&xre[i+96]);
                              zmm25 = _mm512_loadu_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(_mm512_sinh_ps(zmm24),_mm512_cos_ps(zmm25));
                              _mm512_storeu_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(_mm512_cosh_ps(zmm24),_mm512_sin_ps(zmm25));
                              _mm512_storeu_ps(&csim[i+96],zmm27);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+112],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+112],_MM_HINT_T0);
#endif 
                              zmm28 = _mm512_loadu_ps(&xre[i+112]);
                              zmm29 = _mm512_loadu_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(_mm512_sinh_ps(zmm28),_mm512_cos_ps(zmm29));
                              _mm512_storeu_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(_mm512_cosh_ps(zmm28),_mm512_sin_ps(zmm29));
                              _mm512_storeu_ps(&csim[i+112],zmm31);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
#endif                               
                              zmm0 = _mm512_loadu_ps(&xre[i+128]);
                              zmm1 = _mm512_loadu_ps(&xim[i+128]);
                              zmm2 = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+128],zmm2);
                              zmm3 = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+128],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+144],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+144],_MM_HINT_T0);
#endif 
                              zmm4 = _mm512_loadu_ps(&xre[i+144]);
                              zmm5 = _mm512_loadu_ps(&xim[i+144]);
                              zmm6 = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+144],zmm6);
                              zmm7 = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+144],zmm7);
                         }

                          for(; (i+127) < n; i += 128) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+80]);
                              zmm21 = _mm512_loadu_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_storeu_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_storeu_ps(&csim[i+80],zmm23);
                              zmm24 = _mm512_loadu_ps(&xre[i+96]);
                              zmm25 = _mm512_loadu_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(_mm512_sinh_ps(zmm24),_mm512_cos_ps(zmm25));
                              _mm512_storeu_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(_mm512_cosh_ps(zmm24),_mm512_sin_ps(zmm25));
                              _mm512_storeu_ps(&csim[i+96],zmm27);
                              zmm28 = _mm512_loadu_ps(&xre[i+112]);
                              zmm29 = _mm512_loadu_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(_mm512_sinh_ps(zmm28),_mm512_cos_ps(zmm29));
                              _mm512_storeu_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(_mm512_cosh_ps(zmm28),_mm512_sin_ps(zmm29));
                              _mm512_storeu_ps(&csim[i+112],zmm31);
                         }

                          for(; (i+79) < n; i += 80) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = std::sinh(z0)*ceph_cosf(z1);
                              csre[i]        = z2;
                              const float z3 = std::cosh(z0)*ceph_sinf(z1);
                              csim[i]        = z3;
                         }
}


                 
                 
void gms::math::csinhv_svml_zmm16r4_u10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                           const float * __restrict  __ATTR_ALIGN__(64) xim,
                                           float * __restrict __ATTR_ALIGN__(64) csre,
                                           float * __restrict __ATTR_ALIGN__(64) csim,
                                           const int32_t n)
{

                       if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                           __m512 zmm24,zmm25,zmm26,zmm27;
                           __m512 zmm28,zmm29,zmm30,zmm31;
                          int32_t i;
 
                          for(i = 0; (i+159) < n; i += 160) 
                          {
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
#endif                               
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+16],_MM_HINT_T0);
#endif 
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
#endif                               
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
#endif 
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
#endif                              
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
#endif 
                              zmm20 = _mm512_load_ps(&xre[i+80]);
                              zmm21 = _mm512_load_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_store_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_store_ps(&csim[i+80],zmm23);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
#endif                               
                              zmm24 = _mm512_load_ps(&xre[i+96]);
                              zmm25 = _mm512_load_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(_mm512_sinh_ps(zmm24),_mm512_cos_ps(zmm25));
                              _mm512_store_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(_mm512_cosh_ps(zmm24),_mm512_sin_ps(zmm25));
                              _mm512_store_ps(&csim[i+96],zmm27);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+112],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+112],_MM_HINT_T0);
#endif                              
                              zmm28 = _mm512_load_ps(&xre[i+112]);
                              zmm29 = _mm512_load_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(_mm512_sinh_ps(zmm28),_mm512_cos_ps(zmm29));
                              _mm512_store_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(_mm512_cosh_ps(zmm28),_mm512_sin_ps(zmm29));
                              _mm512_store_ps(&csim[i+112],zmm31);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
#endif                            
                              zmm0 = _mm512_load_ps(&xre[i+128]);
                              zmm1 = _mm512_load_ps(&xim[i+128]);
                              zmm2 = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+128],zmm2);
                              zmm3 = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+128],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+144],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+144],_MM_HINT_T0);
#endif                               
                              zmm4 = _mm512_load_ps(&xre[i+144]);
                              zmm5 = _mm512_load_ps(&xim[i+144]);
                              zmm6 = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+144],zmm6);
                              zmm7 = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+144],zmm7);
                         }

                          for(; (i+127) < n; i += 128) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+80]);
                              zmm21 = _mm512_load_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_store_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_store_ps(&csim[i+80],zmm23);
                              zmm24 = _mm512_load_ps(&xre[i+96]);
                              zmm25 = _mm512_load_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(_mm512_sinh_ps(zmm24),_mm512_cos_ps(zmm25));
                              _mm512_store_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(_mm512_cosh_ps(zmm24),_mm512_sin_ps(zmm25));
                              _mm512_store_ps(&csim[i+96],zmm27);
                              zmm28 = _mm512_load_ps(&xre[i+112]);
                              zmm29 = _mm512_load_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(_mm512_sinh_ps(zmm28),_mm512_cos_ps(zmm29));
                              _mm512_store_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(_mm512_cosh_ps(zmm28),_mm512_sin_ps(zmm29));
                              _mm512_store_ps(&csim[i+112],zmm31);
                         }

                          for(; (i+79) < n; i += 80) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = std::sinh(z0)*ceph_cosf(z1);
                              csre[i]        = z2;
                              const float z3 = std::cosh(z0)*ceph_sinf(z1);
                              csim[i]        = z3;
                         }
}


               
void gms::math::csinhv_svml_zmm16r4_u8x_u(const float * __restrict xre,
                                          const float * __restrict xim,
                                          float * __restrict csre,
                                          float * __restrict csim,
                                          int32_t n) 
{

                       if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                           __m512 zmm24,zmm25,zmm26,zmm27;
                           __m512 zmm28,zmm29,zmm30,zmm31;
                          int32_t i;
#if (CSINHV_SVML_ZMM16R4_PEEL_LOOP) == 1
                            while(((uintptr_t)&csre & 63  &&
                                   (uintptr_t)&csim & 63) && n) 
                            {
                                    const float z0 = *xre;
                                    const float z1 = *xim;
                                    const float z2 = std::sinh(z0)*ceph_cosf(z1);
                                    *csre          = z2;
                                    const float z3 = std::cosh(z0)*ceph_sinf(z1);
                                    *csim          = z3;
                                    xre++;
                                    xim++;
                                    csre++;
                                    csim++;
                                    n--;
                            }
#endif  
                          for(i = 0; (i+127) < n; i += 128) {
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
#endif                                
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+16],_MM_HINT_T0);
#endif                                
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
#endif                              
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
#endif  
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
#endif                              
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
#endif  
                              zmm20 = _mm512_loadu_ps(&xre[i+80]);
                              zmm21 = _mm512_loadu_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_storeu_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_storeu_ps(&csim[i+80],zmm23);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
#endif                              
                              zmm24 = _mm512_loadu_ps(&xre[i+96]);
                              zmm25 = _mm512_loadu_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(_mm512_sinh_ps(zmm24),_mm512_cos_ps(zmm25));
                              _mm512_storeu_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(_mm512_cosh_ps(zmm24),_mm512_sin_ps(zmm25));
                              _mm512_storeu_ps(&csim[i+96],zmm27);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+112],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+112],_MM_HINT_T0);
#endif                                
                              zmm28 = _mm512_loadu_ps(&xre[i+112]);
                              zmm29 = _mm512_loadu_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(_mm512_sinh_ps(zmm28),_mm512_cos_ps(zmm29));
                              _mm512_storeu_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(_mm512_cosh_ps(zmm28),_mm512_sin_ps(zmm29));
                              _mm512_storeu_ps(&csim[i+112],zmm31);
                              
                         }

                         
                          for(; (i+79) < n; i += 80) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = std::sinh(z0)*ceph_cosf(z1);
                              csre[i]        = z2;
                              const float z3 = std::cosh(z0)*ceph_sinf(z1);
                              csim[i]        = z3;
                         }
}


                  
void gms::math::csinhv_svml_zmm16r4_u8x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                          const float * __restrict  __ATTR_ALIGN__(64) xim,
                                          float * __restrict __ATTR_ALIGN__(64) csre,
                                          float * __restrict __ATTR_ALIGN__(64) csim,
                                          const int32_t n) 
{

                       if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                           __m512 zmm24,zmm25,zmm26,zmm27;
                           __m512 zmm28,zmm29,zmm30,zmm31;
                          int32_t i;
 
                          for(i = 0; (i+127) < n; i += 128) 
                          {
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
#endif                               
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
#endif  
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
#endif                              
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
#endif  
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
#endif                                
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
#endif  
                              zmm20 = _mm512_load_ps(&xre[i+80]);
                              zmm21 = _mm512_load_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_store_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_store_ps(&csim[i+80],zmm23);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
#endif                               
                              zmm24 = _mm512_load_ps(&xre[i+96]);
                              zmm25 = _mm512_load_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(_mm512_sinh_ps(zmm24),_mm512_cos_ps(zmm25));
                              _mm512_store_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(_mm512_cosh_ps(zmm24),_mm512_sin_ps(zmm25));
                              _mm512_store_ps(&csim[i+96],zmm27);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+112],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+112],_MM_HINT_T0);
#endif  
                              zmm28 = _mm512_load_ps(&xre[i+112]);
                              zmm29 = _mm512_load_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(_mm512_sinh_ps(zmm28),_mm512_cos_ps(zmm29));
                              _mm512_store_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(_mm512_cosh_ps(zmm28),_mm512_sin_ps(zmm29));
                              _mm512_store_ps(&csim[i+112],zmm31);
                           
                         }

                       
                          for(; (i+79) < n; i += 80) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = std::sinh(z0)*ceph_cosf(z1);
                              csre[i]        = z2;
                              const float z3 = std::cosh(z0)*ceph_sinf(z1);
                              csim[i]        = z3;
                         }
}


                 
void gms::math::csinhv_svml_zmm16r4_u6x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                          const float * __restrict  __ATTR_ALIGN__(64) xim,
                                          float * __restrict __ATTR_ALIGN__(64) csre,
                                          float * __restrict __ATTR_ALIGN__(64) csim,
                                          const int32_t n)
{

                       if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                          int32_t i;
 
                          for(i = 0; (i+95) < n; i += 96) 
                          {
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
#endif                             
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+16],_MM_HINT_T0);
#endif
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
#endif                              
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
#endif
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
#endif                             
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
#endif
                              zmm20 = _mm512_load_ps(&xre[i+80]);
                              zmm21 = _mm512_load_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_store_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_store_ps(&csim[i+80],zmm23);
                            
                         }

                          for(; (i+79) < n; i += 80) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) 
                          {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = std::sinh(z0)*ceph_cosf(z1);
                              csre[i]        = z2;
                              const float z3 = std::cosh(z0)*ceph_sinf(z1);
                              csim[i]        = z3;
                         }
}


                
void gms::math::csinhv_svml_zmm16r4_u6x_u(const float * __restrict xre,
                                          const float * __restrict xim,
                                          float * __restrict csre,
                                          float * __restrict csim,
                                          int32_t n)
{

                       if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           __m512 zmm20,zmm21,zmm22,zmm23;
                          int32_t i;
#if (CSINHV_SVML_ZMM16R4_PEEL_LOOP) == 1
                            while(((uintptr_t)&csre & 63  &&
                                   (uintptr_t)&csim & 63) && n) 
                            {
                                    const float z0 = *xre;
                                    const float z1 = *xim;
                                    const float z2 = std::sinh(z0)*ceph_cosf(z1);
                                    *csre          = z2;
                                    const float z3 = std::cosh(z0)*ceph_sinf(z1);
                                    *csim          = z3;
                                    xre++;
                                    xim++;
                                    csre++;
                                    csim++;
                                    n--;
                            }
#endif   
                          for(i = 0; (i+95) < n; i += 96) {
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
#endif                                
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+16],_MM_HINT_T0);
#endif                                 
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
#endif                                
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
#endif                                
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
#endif                                
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
#if (CSINHV_SVML_ZMM16R4_SOFT_PREFETCH) == 1                            
                              _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
#endif   
                              zmm20 = _mm512_loadu_ps(&xre[i+80]);
                              zmm21 = _mm512_loadu_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(_mm512_sinh_ps(zmm20),_mm512_cos_ps(zmm21));
                              _mm512_storeu_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(_mm512_cosh_ps(zmm20),_mm512_sin_ps(zmm21));
                              _mm512_storeu_ps(&csim[i+80],zmm23);
                                                           
                         }

                          for(; (i+79) < n; i += 80) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(_mm512_sinh_ps(zmm16),_mm512_cos_ps(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(_mm512_cosh_ps(zmm16),_mm512_sin_ps(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(_mm512_sinh_ps(zmm8),_mm512_cos_ps(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(_mm512_cosh_ps(zmm8),_mm512_sin_ps(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(_mm512_sinh_ps(zmm12),_mm512_cos_ps(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(_mm512_cosh_ps(zmm12),_mm512_sin_ps(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(_mm512_sinh_ps(zmm4),_mm512_cos_ps(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(_mm512_cosh_ps(zmm4),_mm512_sin_ps(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) 
                          {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(_mm512_sinh_ps(zmm0),_mm512_cos_ps(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(_mm512_cosh_ps(zmm0),_mm512_sin_ps(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) 
                          {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = std::sinh(z0)*ceph_cosf(z1);
                              csre[i]        = z2;
                              const float z3 = std::cosh(z0)*ceph_sinf(z1);
                              csim[i]        = z3;
                         }
}




                 
                  

