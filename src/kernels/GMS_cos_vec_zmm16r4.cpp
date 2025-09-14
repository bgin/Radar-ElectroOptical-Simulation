


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





#include "GMS_cos_vec_zmm16r4.h"
#if (COS_VEC_ZMM16R4_USE_SLEEF) == 1
#include "GMS_sleefsimdsp.h"
#endif 
//#include "GMS_cephes.h"



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
constexpr float FOPI  = 1.27323954473516f;
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


/*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cos_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
                 
void gms::math::cosv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                          float * __restrict __ATTR_ALIGN__(64) y,
                                          const __m512 a,
                                          const __m512 b,
                                          const __m512 c,
                                          const int32_t n) 
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_cos_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts.
            for(i = 0; (i+159) < n; i += 160) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cos_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
#endif 
                            const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cos_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
#endif 
                            const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cos_ps(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
#endif 
                            const __m512 zmm16      = _mm512_load_ps(&x[i+128]); 
                            const __m512 zmm17      = _mm512_cos_ps(zmm16);
                           _mm512_store_ps(&y[i+128], zmm17);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
#endif 
                            const __m512 zmm18      = _mm512_load_ps(&x[i+144]); 
                            const __m512 zmm19      = _mm512_cos_ps(zmm18);
                           _mm512_store_ps(&y[i+144], zmm19);
#else
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_load_ps(&x[i+112]);  
                            const __m512 zmm16     = _mm512_load_ps(&x[i+128]); 
                            const __m512 zmm18     = _mm512_load_ps(&x[i+144]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                            const __m512 zmm11     = _mm512_cos_ps(zmm10);
                            const __m512 zmm13     = _mm512_cos_ps(zmm12);
                            const __m512 zmm15     = _mm512_cos_ps(zmm14);
                            const __m512 zmm17     = _mm512_cos_ps(zmm16);
                            const __m512 zmm19     = _mm512_cos_ps(zmm18);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm512_store_ps(&y[i+144], zmm19);
#endif
            }

            for(; (i+127) < n; i += 128) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cos_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cos_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cos_ps(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                            const __m512 zmm11     = _mm512_cos_ps(zmm10);
                            const __m512 zmm13     = _mm512_cos_ps(zmm12);
                            const __m512 zmm15     = _mm512_cos_ps(zmm14);
                            const __m512 zmm17     = _mm512_cos_ps(zmm16);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
#endif
            }

            for(; (i+63) < n; i += 64) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
            }

            for(; (i+31) < n; i += 32) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
#endif
            }

            for(; (i+15) < n; i += 15) 
            {

                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);

            } 

            for(; (i+0) < n; i += 1) 
            {
                  y[i] = ceph_cosf(x[i]);
           }
}


                 
                  

                 
                  



              /*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cos_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
                  
void gms::math::cosv_zmm16r4_unroll_10x_u(const float * __restrict  x,
                                          float * __restrict  y,
                                          const __m512 a,
                                          const __m512 b,
                                          const __m512 c,
                                          const int32_t n)
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_cos_ps(_mm512_loadu_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts.
            for(i = 0; (i+159) < n; i += 160) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cos_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
#endif 
                            const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cos_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
#endif 
                            const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cos_ps(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
#endif 
                            const __m512 zmm16      = _mm512_loadu_ps(&x[i+128]); 
                            const __m512 zmm17      = _mm512_cos_ps(zmm16);
                           _mm512_storeu_ps(&y[i+128], zmm17);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
#endif 
                            const __m512 zmm18      = _mm512_loadu_ps(&x[i+144]); 
                            const __m512 zmm19      = _mm512_cos_ps(zmm18);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#else
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_loadu_ps(&x[i+112]);  
                            const __m512 zmm16     = _mm512_loadu_ps(&x[i+128]); 
                            const __m512 zmm18     = _mm512_loadu_ps(&x[i+144]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                            const __m512 zmm11     = _mm512_cos_ps(zmm10);
                            const __m512 zmm13     = _mm512_cos_ps(zmm12);
                            const __m512 zmm15     = _mm512_cos_ps(zmm14);
                            const __m512 zmm17     = _mm512_cos_ps(zmm16);
                            const __m512 zmm19     = _mm512_cos_ps(zmm18);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#endif
            }

            for(; (i+127) < n; i += 128) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cos_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cos_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cos_ps(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                            const __m512 zmm11     = _mm512_cos_ps(zmm10);
                            const __m512 zmm13     = _mm512_cos_ps(zmm12);
                            const __m512 zmm15     = _mm512_cos_ps(zmm14);
                            const __m512 zmm17     = _mm512_cos_ps(zmm16);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#endif
            }

            for(; (i+63) < n; i += 64) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
            }

            for(; (i+31) < n; i += 32)
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#endif
            }

            for(; (i+15) < n; i += 15)
            {

                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1      = _mm512_cos_ps(zmm0); 
                           _mm512_storeu_ps(&y[i+0], zmm1);

            } 

            for(; (i+0) < n; i += 1) 
            {
                  y[i] = ceph_cosf(x[i]);
            }
}

               
               /*
                    Calls non-SVML implementation of sine function
                    SLEEF version is inlined.
                */
#if (COS_VEC_ZMM16R4_USE_SLEEF) == 1

                   void gms::math::cosv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                  float * __restrict __ATTR_ALIGN__(64) y,
                                                  const int32_t n) {

                       if(__builtin_expect(0==n,0)) {return;}
                       int32_t i;
                      
                       for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                            const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                            const vfloat zmm9      = xcosf(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const vfloat zmm10      = _mm512_load_ps(&x[i+80]); 
                            const vfloat zmm11      = xcosf(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                            const vfloat zmm12      = _mm512_load_ps(&x[i+96]); 
                            const vfloat zmm13      = xcosf(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                            const vfloat zmm14      = _mm512_load_ps(&x[i+112]); 
                            const vfloat zmm15      = xcosf(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
                            const vfloat zmm16      = _mm512_load_ps(&x[i+128]); 
                            const vfloat zmm17      = xcosf(zmm16);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const vfloat zmm18      = _mm512_load_ps(&x[i+144]); 
                            const vfloat zmm19      = xcosf(zmm18);
                           _mm512_store_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                            const vfloat zmm10     = _mm512_load_ps(&x[i+80]);
                            const vfloat zmm12     = _mm512_load_ps(&x[i+96]); 
                            const vfloat zmm14     = _mm512_load_ps(&x[i+112]);  
                            const vfloat zmm16     = _mm512_load_ps(&x[i+128]); 
                            const vfloat zmm18     = _mm512_load_ps(&x[i+144]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                            const vfloat zmm9      = xcosf(zmm8);
                            const vfloat zmm11     = xcosf(zmm10);
                            const vfloat zmm13     = xcosf(zmm12);
                            const vfloat zmm15     = xcosf(zmm14);
                            const vfloat zmm17     = xcosf(zmm16);
                            const vfloat zmm19     = xcosf(zmm18);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm512_store_ps(&y[i+144], zmm19);
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                            const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                            const vfloat zmm9      = xcosf(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const vfloat zmm10      = _mm512_load_ps(&x[i+80]); 
                            const vfloat zmm11      = xcosf(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                            const vfloat zmm12      = _mm512_load_ps(&x[i+96]); 
                            const vfloat zmm13      = xcosf(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                            const vfloat zmm14      = _mm512_load_ps(&x[i+112]); 
                            const vfloat zmm15      = xcosf(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
#else
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                            const vfloat zmm10     = _mm512_load_ps(&x[i+80]);
                            const vfloat zmm12     = _mm512_load_ps(&x[i+96]); 
                            const vfloat zmm14     = _mm512_load_ps(&x[i+112]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                            const vfloat zmm9      = xcosf(zmm8);
                            const vfloat zmm11     = xcosf(zmm10);
                            const vfloat zmm13     = xcosf(zmm12);
                            const vfloat zmm15     = xcosf(zmm14);
                            const vfloat zmm17     = xcosf(zmm16);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#else
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]);
                            const vfloat zmm1      = xcosf(zmm0); 
                           _mm512_store_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = ceph_cosf(x[i]);
                     }
               } 

#endif 


#if (COS_VEC_ZMM16R4_USE_SLEEF) == 1

                   void gms::math::cosv_zmm16r4_unroll_10x_u(const float * __restrict  x,
                                                  float * __restrict  y,
                                                  const int32_t n) {

                       if(__builtin_expect(0==n,0)) {return;}
                       int32_t i;
                      
                       for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                            const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const vfloat zmm9      = xcosf(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const vfloat zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const vfloat zmm11      = xcosf(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                            const vfloat zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const vfloat zmm13      = xcosf(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                            const vfloat zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const vfloat zmm15      = xcosf(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                            const vfloat zmm16      = _mm512_loadu_ps(&x[i+128]); 
                            const vfloat zmm17      = xcosf(zmm16);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const vfloat zmm18      = _mm512_loadu_ps(&x[i+144]); 
                            const vfloat zmm19      = xcosf(zmm18);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const vfloat zmm10     = _mm512_loadu_ps(&x[i+80]);
                            const vfloat zmm12     = _mm512_loadu_ps(&x[i+96]); 
                            const vfloat zmm14     = _mm512_loadu_ps(&x[i+112]);  
                            const vfloat zmm16     = _mm512_loadu_ps(&x[i+128]); 
                            const vfloat zmm18     = _mm512_loadu_ps(&x[i+144]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                            const vfloat zmm9      = xcosf(zmm8);
                            const vfloat zmm11     = xcosf(zmm10);
                            const vfloat zmm13     = xcosf(zmm12);
                            const vfloat zmm15     = xcosf(zmm14);
                            const vfloat zmm17     = xcosf(zmm16);
                            const vfloat zmm19     = xcosf(zmm18);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                            const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const vfloat zmm9      = xcosf(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const vfloat zmm10      = _mm512_load_ps(&x[i+80]); 
                            const vfloat zmm11      = xcosf(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                            const vfloat zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const vfloat zmm13      = xcosf(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                            const vfloat zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const vfloat zmm15      = xcosf(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#else
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const vfloat zmm10     = _mm512_loadu_ps(&x[i+80]);
                            const vfloat zmm12     = _mm512_loadu_ps(&x[i+96]); 
                            const vfloat zmm14     = _mm512_loadu_ps(&x[i+112]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                            const vfloat zmm9      = xcosf(zmm8);
                            const vfloat zmm11     = xcosf(zmm10);
                            const vfloat zmm13     = xcosf(zmm12);
                            const vfloat zmm15     = xcosf(zmm14);
                            const vfloat zmm17     = xcosf(zmm16);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#else
                            const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]);
                            const vfloat zmm1      = xcosf(zmm0); 
                           _mm512_storeu_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = ceph_cosf(x[i]);
                     }
               } 

#endif 

 /*             
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cos_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/

                 
void gms::math::cosv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                         float * __restrict __ATTR_ALIGN__(64) y,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 c,
                                         const int32_t n) 
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_cos_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts here
            for(i = 0; (i+95) < n; i += 96) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cos_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
#else
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm10     = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                            const __m512 zmm11     = _mm512_cos_ps(zmm10);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
#endif
            }

            for(; (i+63) < n; i += 64) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
            }

            for(; (i+31) < n; i += 32) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3); 
#endif

            }

            for(; (i+15) < n; i += 16) 
            {

                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);

            }

            for(; (i+0) < n; i += 1) 
            {
                           y[i] = ceph_cosf(x[i]);
            }
}


               /*             
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cos_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/

                
void gms::math::cosv_zmm16r4_unroll_6x_u(const float * __restrict  x,
                                         float * __restrict  y,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 c,
                                         const int32_t n)
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_cos_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts here
            for(i = 0; (i+95) < n; i += 96) {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cos_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#else
#if (COS_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                            const __m512 zmm9      = _mm512_cos_ps(zmm8);
                            const __m512 zmm11     = _mm512_cos_ps(zmm10);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#endif
            }

            for(; (i+63) < n; i += 64) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                            const __m512 zmm5      = _mm512_cos_ps(zmm4);
                            const __m512 zmm7      = _mm512_cos_ps(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
            }

            for(; (i+31) < n; i += 32) 
            {
#if (COS_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                            const __m512 zmm3      = _mm512_cos_ps(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3); 
#endif

            }

            for(; (i+15) < n; i += 16) 
            {

                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cos_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);

            }

            for(; (i+0) < n; i += 1) 
            {
                  y[i] = ceph_cosf(x[i]);
            }
}


                /*
                    Calls non-SVML implementation of sine function
                    SLEEF version is inlined.
                */
#if (COS_VEC_ZMM16R4_USE_SLEEF) == 1

                   void gms::math::cosv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                  float * __restrict __ATTR_ALIGN__(64) y,
                                                  const int32_t n) {

                       if(__builtin_expect(0==n,0)) {return;}
                       int32_t i;
                      
                       for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                            const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                            const vfloat zmm9      = xcosf(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                            const vfloat zmm10      = _mm512_load_ps(&x[i+80]); 
                            const vfloat zmm11      = xcosf(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                          
#else
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                            const vfloat zmm10     = _mm512_load_ps(&x[i+80]);
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                            const vfloat zmm9      = xcosf(zmm8);
                            const vfloat zmm11     = xcosf(zmm10);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm5      = xcosf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                          
#else
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                            const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                            const vfloat zmm5      = xcosf(zmm4);
                            const vfloat zmm7      = xcosf(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                           
#else
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                            const vfloat zmm1      = xcosf(zmm0);
                            const vfloat zmm3      = xcosf(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                          
#endif
                      }

                     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                          
#else
                            const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                            const vfloat zmm1      = xcosf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#endif
                     }

                   
                     for(; (i+0) < n; i += 1) {
                           y[i] = ceph_cosf(x[i]);
                     }
               } 

#endif 









              

                    

  
