


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
#include "GMS_csub_vec_zmm16r4.h"



void gms::math::csubv_zmm16r4_unroll_16x_u(const float * __restrict  xre,
                                           const float * __restrict  xim,
                                           const float * __restrict  yre,
                                           const float * __restrict  yim,
                                           float * __restrict        zre,
                                           float * __restrict        zim,
                                           int32_t n)
{

                        if(__builtin_expect(0==n,0)) {return;}
#if (CSUB_VEC_ZMM16R4_ADD_PEEL_LOOP) == 1
                           while(((uintptr_t)&zre & 63 && 
                                 (uintptr_t)&zim & 63) && n) 
                           {
                                 const float xr = *xre;
                                 const float xi = *xim;
                                 const float yr = *yre;
                                 const float yi = *yim;
                                 *zre           = xr-yr;
                                 *zim           = xi-yi;
                                 xre++;
                                 xim++;
                                 yre++;
                                 yim++;
                                 zre++;
                                 zim++;
                                 n--;
                           }
#endif
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                         __m512 zmm24,zmm25,zmm26,zmm27;
                         __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) 
                        {
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+0],_MM_HINT_T0);
#endif
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+16],_MM_HINT_T0);
#endif
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+32],_MM_HINT_T0);
#endif
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
#endif
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
#endif
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
#endif                          
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
#endif
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+112],_MM_HINT_T0);
#endif                            
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
#endif
                            zmm0 = _mm512_loadu_ps(&xre[i+128]); 
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm0,zmm1));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
#endif                         
                            zmm4 = _mm512_loadu_ps(&xre[i+144]); 
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm4,zmm5));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+160],_MM_HINT_T0);
#endif                            
                            zmm8 = _mm512_loadu_ps(&xre[i+160]); 
                            zmm9 = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+160]);
                            zmm11 = _mm512_loadu_ps(&yim[i+160]);
                            _mm512_storeu_ps(&zim[i+160], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+176],_MM_HINT_T0);
#endif                            
                            zmm12 = _mm512_loadu_ps(&xre[i+176]); 
                            zmm13 = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+176]);
                            zmm15 = _mm512_loadu_ps(&yim[i+176]);
                            _mm512_storeu_ps(&zim[i+176], _mm512_sub_ps(zmm14,zmm15));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
#endif
                            zmm16 = _mm512_loadu_ps(&xre[i+192]); 
                            zmm17 = _mm512_loadu_ps(&yre[i+192]);
                            _mm512_storeu_ps(&zre[i+192], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+192]);
                            zmm19 = _mm512_loadu_ps(&yim[i+192]);
                            _mm512_storeu_ps(&zim[i+192], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+208],_MM_HINT_T0);
#endif                           
                            zmm20 = _mm512_loadu_ps(&xre[i+208]); 
                            zmm21 = _mm512_loadu_ps(&yre[i+208]);
                            _mm512_storeu_ps(&zre[i+208], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+208]);
                            zmm23 = _mm512_loadu_ps(&yim[i+208]);
                            _mm512_storeu_ps(&zim[i+208], _mm512_sub_ps(zmm22,zmm23));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+224],_MM_HINT_T0);
#endif                            
                            zmm24 = _mm512_loadu_ps(&xre[i+224]); 
                            zmm25 = _mm512_loadu_ps(&yre[i+224]);
                            _mm512_storeu_ps(&zre[i+224], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+224]);
                            zmm27 = _mm512_loadu_ps(&yim[i+224]);
                            _mm512_storeu_ps(&zim[i+224], _mm512_sub_ps(zmm26,zmm27));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+240],_MM_HINT_T0);
#endif                            
                            zmm28 = _mm512_loadu_ps(&xre[i+240]); 
                            zmm29 = _mm512_loadu_ps(&yre[i+240]);
                            _mm512_storeu_ps(&zre[i+240], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+240]);
                            zmm31 = _mm512_loadu_ps(&yim[i+240]);
                            _mm512_storeu_ps(&zim[i+240], _mm512_sub_ps(zmm30,zmm31));
                        
                        }

                       for(; (i+191) < n; i += 192) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            zmm0 = _mm512_loadu_ps(&xre[i+128]); 
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm0,zmm1));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm2,zmm3));
                            zmm4 = _mm512_loadu_ps(&xre[i+144]); 
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm4,zmm5));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm6,zmm7));
                            zmm8 = _mm512_loadu_ps(&xre[i+160]); 
                            zmm9 = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+160]);
                            zmm11 = _mm512_loadu_ps(&yim[i+160]);
                            _mm512_storeu_ps(&zim[i+160], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+176]); 
                            zmm13 = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+176]);
                            zmm15 = _mm512_loadu_ps(&yim[i+176]);
                            _mm512_storeu_ps(&zim[i+176], _mm512_sub_ps(zmm14,zmm15));

                       }

                       for(; (i+127) < n; i += 128) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));

                       }

                       for(; (i+63) < n; i += 64) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));

                       }

                       for(; (i+31) < n; i += 32) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));

                       }

                      for(; (i+15) < n; i += 16) 
                      {
                           zmm0  = _mm512_loadu_ps(&xre[i+0]);
                           zmm1  = _mm512_loadu_ps(&yre[i+0]);
                           _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                           zmm2  = _mm512_loadu_ps(&xim[i+0]);
                           zmm3  = _mm512_loadu_ps(&yim[i+0]);
                           _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));

                      }

                      for(; (i+0) < n; i += 1) {
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }
}


                  
void gms::math::csubv_zmm16r4_unroll_16x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                           const float * __restrict  __ATTR_ALIGN__(64) xim,
                                           const float * __restrict  __ATTR_ALIGN__(64) yre,
                                           const float * __restrict  __ATTR_ALIGN__(64) yim,
                                           float * __restrict        __ATTR_ALIGN__(64) zre,
                                           float * __restrict        __ATTR_ALIGN__(64) zim,
                                           const int32_t n)
{

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                         __m512 zmm24,zmm25,zmm26,zmm27;
                         __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) {
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+0],_MM_HINT_T0);
#endif                           
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+16],_MM_HINT_T0);
#endif                            
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+32],_MM_HINT_T0);
#endif
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
#endif                            
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
#endif                            
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
#endif                          
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
#endif 
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+112],_MM_HINT_T0);
#endif 
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
#endif 
                            zmm0 = _mm512_load_ps(&xre[i+128]); 
                            zmm1 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm0,zmm1));
                            zmm2 = _mm512_load_ps(&xim[i+128]);
                            zmm3 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
#endif                         
                            zmm4 = _mm512_load_ps(&xre[i+144]); 
                            zmm5 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm4,zmm5));
                            zmm6 = _mm512_load_ps(&xim[i+144]);
                            zmm7 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+160],_MM_HINT_T0);
#endif                             
                            zmm8 = _mm512_load_ps(&xre[i+160]); 
                            zmm9 = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+160]);
                            zmm11 = _mm512_load_ps(&yim[i+160]);
                            _mm512_store_ps(&zim[i+160], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+176],_MM_HINT_T0);
#endif 
                            zmm12 = _mm512_load_ps(&xre[i+176]); 
                            zmm13 = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+176]);
                            zmm15 = _mm512_load_ps(&yim[i+176]);
                            _mm512_store_ps(&zim[i+176], _mm512_sub_ps(zmm14,zmm15));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
#endif 
                            zmm16 = _mm512_load_ps(&xre[i+192]); 
                            zmm17 = _mm512_load_ps(&yre[i+192]);
                            _mm512_store_ps(&zre[i+192], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+192]);
                            zmm19 = _mm512_load_ps(&yim[i+192]);
                            _mm512_store_ps(&zim[i+192], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+208],_MM_HINT_T0);
#endif                            
                            zmm20 = _mm512_load_ps(&xre[i+208]); 
                            zmm21 = _mm512_load_ps(&yre[i+208]);
                            _mm512_store_ps(&zre[i+208], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+208]);
                            zmm23 = _mm512_load_ps(&yim[i+208]);
                            _mm512_store_ps(&zim[i+208], _mm512_sub_ps(zmm22,zmm23));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+224],_MM_HINT_T0);
#endif                             
                            zmm24 = _mm512_load_ps(&xre[i+224]); 
                            zmm25 = _mm512_load_ps(&yre[i+224]);
                            _mm512_store_ps(&zre[i+224], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+224]);
                            zmm27 = _mm512_load_ps(&yim[i+224]);
                            _mm512_store_ps(&zim[i+224], _mm512_sub_ps(zmm26,zmm27));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+240],_MM_HINT_T0);
#endif 
                            zmm28 = _mm512_load_ps(&xre[i+240]); 
                            zmm29 = _mm512_load_ps(&yre[i+240]);
                            _mm512_store_ps(&zre[i+240], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+240]);
                            zmm31 = _mm512_load_ps(&yim[i+240]);
                            _mm512_store_ps(&zim[i+240], _mm512_sub_ps(zmm30,zmm31));
                        
                        }

                       for(; (i+191) < n; i += 192) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            zmm0 = _mm512_load_ps(&xre[i+128]); 
                            zmm1 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm0,zmm1));
                            zmm2 = _mm512_load_ps(&xim[i+128]);
                            zmm3 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm2,zmm3));
                            zmm4 = _mm512_load_ps(&xre[i+144]); 
                            zmm5 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm4,zmm5));
                            zmm6 = _mm512_load_ps(&xim[i+144]);
                            zmm7 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm6,zmm7));
                            zmm8 = _mm512_load_ps(&xre[i+160]); 
                            zmm9 = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+160]);
                            zmm11 = _mm512_load_ps(&yim[i+160]);
                            _mm512_store_ps(&zim[i+160], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+176]); 
                            zmm13 = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+176]);
                            zmm15 = _mm512_load_ps(&yim[i+176]);
                            _mm512_store_ps(&zim[i+176], _mm512_sub_ps(zmm14,zmm15));
                       }

                       for(; (i+127) < n; i += 128) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));

                       }

                       for(; (i+63) < n; i += 64) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                       }

                       for(; (i+31) < n; i += 32) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));

                       }

                      for(; (i+15) < n; i += 16) 
                      {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));

                      }

                      for(; (i+0) < n; i += 1) 
                      {
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }
}

                    

                  
void gms::math::csubv_zmm16r4_unroll_10x_u(const float * __restrict  xre,
                                           const float * __restrict  xim,
                                           const float * __restrict  yre,
                                           const float * __restrict  yim,
                                           float * __restrict        zre,
                                           float * __restrict        zim,
                                           int32_t n) 
{

                        if(__builtin_expect(0==n,0)) {return;}
#if (CSUB_VEC_ZMM16R4_ADD_PEEL_LOOP) == 1
                          while(((uintptr_t)&zre & 63 &&
                                (uintptr_t)&zim &63) && n)
                         {
                                const float xr = *xre;
                                const float xi = *xim;
                                const float yr = *yre;
                                const float yi = *yim;
                                *zre           = xr-yr;
                                *zim           = xi-yi;
                                xre++;
                                xim++;
                                yre++;
                                yim++;
                                zre++;
                                zim++;
                                n--;
                         }
#endif
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
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+0],_MM_HINT_T0);
#endif                                
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+16],_MM_HINT_T0);
#endif                                  
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+32],_MM_HINT_T0);
#endif                                   
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
#endif                                   
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
#endif                                  
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
#endif       
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
#endif                                  
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+112],_MM_HINT_T0);
#endif                                   
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
#endif                                   
                            zmm0 = _mm512_loadu_ps(&xre[i+128]); 
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm0,zmm1));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
#endif                                  
                            zmm4 = _mm512_loadu_ps(&xre[i+144]); 
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm4,zmm5));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm6,zmm7));
                        }


                       for(; (i+127) < n; i += 128) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));

                       }

                       for(; (i+79) < n; i += 80) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                                                 
                       }

                       for(; (i+63) < n; i += 64) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                      }

                      for(; (i+31) < n; i += 32) 
                      {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                      }

                     for(; (i+15) < n; i += 16) 
                     {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                     }

                     for(; (i+0) < n; i += 1) 
                     {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }

}


                
void gms::math::csubv_zmm16r4_unroll_10x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                           const float * __restrict  __ATTR_ALIGN__(64) xim,
                                           const float * __restrict  __ATTR_ALIGN__(64) yre,
                                           const float * __restrict  __ATTR_ALIGN__(64) yim,
                                           float * __restrict        __ATTR_ALIGN__(64) zre,
                                           float * __restrict        __ATTR_ALIGN__(64) zim,
                                           const int32_t n) 
{

                        if(__builtin_expect(0==n,0)) {return;}
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
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+0],_MM_HINT_T0);
#endif 
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+16],_MM_HINT_T0);
#endif 
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+32],_MM_HINT_T0);
#endif 
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
#endif                           
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
#endif 
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
#endif 
                             zmm20 = _mm512_load_ps(&xre[i+80]);
                             zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                             zmm22 = _mm512_load_ps(&xim[i+80]);
                             zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
#endif                           
                             zmm24 = _mm512_load_ps(&xre[i+96]);
                             zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                             zmm26 = _mm512_load_ps(&xim[i+96]);
                             zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+112],_MM_HINT_T0);
#endif 
                             zmm28 = _mm512_load_ps(&xre[i+112]);
                             zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                             zmm30 = _mm512_load_ps(&xim[i+112]);
                             zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
#endif 
                             zmm2 = _mm512_load_ps(&xre[i+128]);
                             zmm3 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm2,zmm3));
                             zmm4 = _mm512_load_ps(&xim[i+128]);
                             zmm5 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm4,zmm5));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
#endif 
                             zmm6 = _mm512_load_ps(&xre[i+144]); 
                             zmm7 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm6,zmm7));
                             zmm8 = _mm512_load_ps(&xim[i+144]);
                             zmm9 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm8,zmm9));

                        }

                       for(; (i+127) < n; i += 128) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));

                       }

                       for(; (i+79) < n; i += 80) 
                       {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                             zmm4  = _mm512_load_ps(&xre[i+16]);
                             zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                             zmm8  = _mm512_load_ps(&xre[i+32]);
                             zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                             zmm10 = _mm512_load_ps(&xim[i+32]);
                             zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                             zmm12 = _mm512_load_ps(&xre[i+48]);
                             zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                             zmm14 = _mm512_load_ps(&xim[i+48]);
                             zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                             zmm16 = _mm512_load_ps(&xre[i+64]);
                             zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                             zmm18 = _mm512_load_ps(&xim[i+64]);
                             zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));

                       }

                       for(; (i+63) < n; i += 64) 
                       {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                             zmm4  = _mm512_load_ps(&xre[i+16]);
                             zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                             zmm8  = _mm512_load_ps(&xre[i+32]);
                             zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                             zmm10 = _mm512_load_ps(&xim[i+32]);
                             zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                             zmm12 = _mm512_load_ps(&xre[i+48]);
                             zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                             zmm14 = _mm512_load_ps(&xim[i+48]);
                             zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));

                      }

                      for(; (i+31) < n; i += 32) 
                      {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                             zmm4  = _mm512_load_ps(&xre[i+16]);
                             zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));

                      }

                     for(; (i+15) < n; i += 16) 
                     {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));

                     }

                     for(; (i+0) < n; i += 1) 
                     {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }

}
 


                  
void gms::math::csubv_zmm16r4_unroll_6x_u(const float * __restrict  xre,
                                          const float * __restrict  xim,
                                          const float * __restrict  yre,
                                          const float * __restrict  yim,
                                          float * __restrict        zre,
                                          float * __restrict        zim,
                                          int32_t n) 
{

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                        int32_t i;
#if (CSUB_VEC_ZMM16R4_ADD_PEEL_LOOP) == 1
                        while(((uintptr_t)&zre & 63  && 
                               (uintptr_t)&zim & 63) && n) 
                        {
                              const float xr = *xre;
                              const float xi = *xim;
                              const float yr = *yre;
                              const float yi = *yim;
                              *zre           = xr-yr;
                              *zim           = xi-yi;
                              xre++;
                              xim++;
                              yre++;
                              yim++;
                              zre++;
                              zim++;
                              n--;
                        }
#endif 
                        for(i = 0; (i+95) < n; i += 96) 
                        {
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+0],_MM_HINT_T0);
#endif 
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+16],_MM_HINT_T0);
#endif                             
                             zmm4  = _mm512_loadu_ps(&xre[i+16]);
                             zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+32],_MM_HINT_T0);
#endif 
                             zmm8  = _mm512_loadu_ps(&xre[i+32]);
                             zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                             zmm10 = _mm512_loadu_ps(&xim[i+32]);
                             zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
#endif 
                             zmm12 = _mm512_loadu_ps(&xre[i+48]);
                             zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                             zmm14 = _mm512_loadu_ps(&xim[i+48]);
                             zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
#endif                            
                             zmm16 = _mm512_loadu_ps(&xre[i+64]);
                             zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                             zmm18 = _mm512_loadu_ps(&xim[i+64]);
                             zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
#endif                             
                             zmm20 = _mm512_loadu_ps(&xre[i+80]);
                             zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                             zmm22 = _mm512_loadu_ps(&xim[i+80]);
                             zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                        }

                       for(; (i+63) < n; i += 64) 
                       {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                             zmm4  = _mm512_loadu_ps(&xre[i+16]);
                             zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                             zmm8  = _mm512_loadu_ps(&xre[i+32]);
                             zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                             zmm10 = _mm512_loadu_ps(&xim[i+32]);
                             zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                             zmm12 = _mm512_loadu_ps(&xre[i+48]);
                             zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                             zmm14 = _mm512_loadu_ps(&xim[i+48]);
                             zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                       }

                       for(; (i+31) < n; i += 32) 
                       {

                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                             zmm4  = _mm512_loadu_ps(&xre[i+16]);
                             zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                       }

                       for(; (i+15) < n; i += 16) 
                       {

                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));

                       }

                        for(; (i+0) < n; i += 1) 
                        {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }

}


                 
void gms::math::csubv_zmm16r4_unroll_6x_a( const float * __restrict  __ATTR_ALIGN__(64) xre,
                                           const float * __restrict  __ATTR_ALIGN__(64) xim,
                                           const float * __restrict  __ATTR_ALIGN__(64) yre,
                                           const float * __restrict  __ATTR_ALIGN__(64) yim,
                                           float * __restrict        __ATTR_ALIGN__(64) zre,
                                           float * __restrict        __ATTR_ALIGN__(64) zim,
                                           const int32_t n) 
{

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 zmm20,zmm21,zmm22,zmm23;
                        int32_t i;

                        for(i = 0; (i+95) < n; i += 96) 
                        {
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+0],_MM_HINT_T0);
#endif 
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+16],_MM_HINT_T0);
#endif 
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+32],_MM_HINT_T0);
#endif 
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
#endif 
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
#endif                           
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+80],_MM_HINT_T0);
#endif                            
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));

                        }

                       for(; (i+63) < n; i += 64) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                             zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));

                       }

                       for(; (i+31) < n; i += 32) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));

                       }

                       for(; (i+15) < n; i += 16) 
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));

                       }

                        for(; (i+0) < n; i += 1) 
                        {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                       }

}

                  
                  
void gms::math::csubv_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                           const float               s,
                                           float  * __restrict zre,
                                           const int32_t n) 
{

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm2;
                           __m512 zmm4,zmm6;
                           __m512 zmm8,zmm10;
                           __m512 zmm12,zmm14;
                           __m512 zmm16,zmm18;
                           __m512 zmmx    = _mm512_set1_ps(s);
                           int32_t i;
                         for(i = 0; (i+159) < n; i += 160) 
                         {
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
#endif 
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
#endif 
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
#endif 
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
#endif                             
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
#endif 
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
#endif 
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
#endif 
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
#endif 
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
#endif 
                              zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
#endif 
                              zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                       }

                       for(; (i+127) < n; i += 128) 
                       {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              
                       }

                       for(; (i+79) < n; i += 80) 
                       {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                             
                       }

                      for(; (i+63) < n; i += 64) 
                      {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              
                      }

                     for(; (i+31) < n; i += 32) 
                     {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                             
                     }

                    for(; (i+15) < n; i += 16) 
                    {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                               
                    }

                    for(; (i+0) < n; i += 1) 
                    {
                          
                          zre[i] = xre[i] - s;
                        
                    }

}


                  
void gms::math::csubv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                           const float            s,
                                           float  * __restrict __ATTR_ALIGN__(64) zre,
                                           const int32_t n) 
{

                         if(__builtin_expect(0==n,0)) { return;}
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           __m512 zmm16,zmm17,zmm18,zmm19;
                           const __m512 zmmx    = _mm512_set1_ps(s);
                           int32_t i;
                          for(i = 0; (i+159) < n; i += 160) 
{
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
#endif 
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
#endif 
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
#endif 
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
#endif                              
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
#endif 
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
#endif                               
                               zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
#endif 
                               zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
#endif 
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
#endif 
                               zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                              _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                             
#endif 
                              zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                             
                        }

                       for(; (i+127) < n; i += 128) 
                       {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                             
                       }

                       for(; (i+79) < n; i += 80) 
                       {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              
                       }

                      for(; (i+63) < n; i += 64) 
                      {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                               
                      }

                     for(; (i+31) < n; i += 32) 
                     {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                             
                     }

                    for(; (i+15) < n; i += 16) 
                    {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              
                    }

                    for(; (i+0) < n; i += 1) 
                    {
                         
                          zre[i] = xre[i] - s;
                         
                    }

}

                 
void gms::math::csubv_zmm16r4_unroll_16x_u(const float * __restrict  xre,
                                           const float * __restrict  yre,
                                           float * __restrict        zre,
                                           int32_t n) 
{

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm18,zmm19;
                         __m512 zmm21,zmm22;
                         __m512 zmm24,zmm25,zmm27;
                         __m512 zmm28,zmm30,zmm31;
                        int32_t i;
                        while(((uintptr_t)&zre & 63) && n)
                        {
                               const float xr = *xre;
                               const float yr = *yre;
                               *zre           = xr-yr;
                               xre++;
                               yre++;
                               zre++;
                               n--;
                        }
                        for(i = 0; (i+255) < n; i += 256) 
                        {
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
#endif 
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
#endif                             
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
#endif                             
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
#endif                            
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
#endif                            
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm12,zmm13));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
#endif                            
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm15,zmm16));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
#endif                             
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+112],_MM_HINT_T0);
#endif                             
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm21,zmm22));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
#endif                            
                            zmm24  = _mm512_loadu_ps(&xre[i+128]);
                            zmm25  = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm25,zmm24));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
#endif                             
                            zmm27  = _mm512_loadu_ps(&xre[i+144]);
                            zmm28  = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm27,zmm28));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+160],_MM_HINT_T0);
#endif                            
                            zmm30  = _mm512_loadu_ps(&xre[i+160]);
                            zmm31  = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm30,zmm31));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+176],_MM_HINT_T0);
#endif                             
                            zmm2  = _mm512_loadu_ps(&xre[i+176]);
                            zmm3  = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
#endif                                                        
                            zmm5  = _mm512_loadu_ps(&xre[i+192]);
                            zmm6  = _mm512_loadu_ps(&yre[i+192]);
                            _mm512_storeu_ps(&zre[i+192], _mm512_sub_ps(zmm5,zmm6));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+208],_MM_HINT_T0);
#endif                             
                            zmm8  = _mm512_loadu_ps(&xre[i+208]);
                            zmm9  = _mm512_loadu_ps(&yre[i+208]);
                            _mm512_storeu_ps(&zre[i+208], _mm512_sub_ps(zmm8,zmm9));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+224],_MM_HINT_T0);
#endif                             
                            zmm11 = _mm512_loadu_ps(&xre[i+224]);
                            zmm12 = _mm512_loadu_ps(&yre[i+224]);
                            _mm512_storeu_ps(&zre[i+224], _mm512_sub_ps(zmm11,zmm12));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+240],_MM_HINT_T0);
#endif                            
                            zmm14  = _mm512_loadu_ps(&xre[i+240]);
                            zmm15  = _mm512_loadu_ps(&yre[i+240]);
                            _mm512_storeu_ps(&zre[i+240], _mm512_sub_ps(zmm14,zmm15));
                            
                        }

                       for(; (i+191) < n; i += 192) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm12,zmm13));
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm15,zmm16));
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm18,zmm19));
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm21,zmm22));
                            zmm24  = _mm512_loadu_ps(&xre[i+128]);
                            zmm25  = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm25,zmm24));
                            zmm27  = _mm512_loadu_ps(&xre[i+144]);
                            zmm28  = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm27,zmm28));
                            zmm30  = _mm512_loadu_ps(&xre[i+160]);
                            zmm31  = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm30,zmm31));
                            zmm2  = _mm512_loadu_ps(&xre[i+176]);
                            zmm3  = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm2,zmm3));
                            
                       }

                       for(; (i+127) < n; i += 128)
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
                            zmm12  = _mm512_loadu_ps(&xre[i+64]);
                            zmm13  = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm12,zmm13));
                            zmm15  = _mm512_loadu_ps(&xre[i+80]);
                            zmm16  = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm15,zmm16));
                            zmm18  = _mm512_loadu_ps(&xre[i+96]);
                            zmm19  = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm18,zmm19));
                            zmm21  = _mm512_loadu_ps(&xre[i+112]);
                            zmm22  = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm21,zmm22));
                            
                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            zmm6  = _mm512_loadu_ps(&xre[i+32]);
                            zmm7  = _mm512_loadu_ps(&yre[i+32]); //yre
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
                            zmm9  = _mm512_loadu_ps(&xre[i+48]);
                            zmm10 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
                            

                       }

                       for(; (i+31) < n; i += 32) 
                       {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_loadu_ps(&xre[i+16]);
                            zmm4  = _mm512_loadu_ps(&yre[i+16]); //yre
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            

                       }

                      for(; (i+15) < n; i += 16) 
                      {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]); //yre
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            

                      }

                      for(; (i+0) < n; i += 1) 
                      {
                            zre[i] = xre[i] - yre[i];
                           
                      }
}


                  
void gms::math::csubv_zmm16r4_unroll_16x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                           const float * __restrict  __ATTR_ALIGN__(64) yre,
                                           float * __restrict        __ATTR_ALIGN__(64) zre,
                                           const int32_t n) 
{

                        if(__builtin_expect(0==n,0)) {return;}
                         __m512 zmm0,zmm1,zmm2,zmm3;
                         __m512 zmm4,zmm5,zmm6,zmm7;
                         __m512 zmm8,zmm9,zmm10,zmm11;
                         __m512 zmm12,zmm13,zmm14,zmm15;
                         __m512 zmm16,zmm18,zmm19;
                         __m512 zmm21,zmm22;
                         __m512 zmm24,zmm25,zmm27;
                         __m512 zmm28,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) 
                        {
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+0],_MM_HINT_T0);
#endif                           
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+16],_MM_HINT_T0);
#endif  
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+32],_MM_HINT_T0);
#endif                              
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
#endif                              
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
#endif                              
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm12,zmm13));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+80],_MM_HINT_T0);
#endif                              
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm15,zmm16));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
#endif                              
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm18,zmm19));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+112],_MM_HINT_T0);
#endif                              
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm21,zmm22));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
#endif                              
                            zmm24  = _mm512_load_ps(&xre[i+128]);
                            zmm25  = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm25,zmm24));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
#endif                              
                            zmm27  = _mm512_load_ps(&xre[i+144]);
                            zmm28  = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm27,zmm28));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+160],_MM_HINT_T0);
#endif                              
                            zmm30  = _mm512_load_ps(&xre[i+160]);
                            zmm31  = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm30,zmm31));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+176],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+176],_MM_HINT_T0);
#endif                             
                            zmm2  = _mm512_load_ps(&xre[i+176]);
                            zmm3  = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm2,zmm3));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
#endif                              
                            zmm5  = _mm512_load_ps(&xre[i+192]);
                            zmm6  = _mm512_load_ps(&yre[i+192]);
                            _mm512_store_ps(&zre[i+192], _mm512_sub_ps(zmm5,zmm6));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+208],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+208],_MM_HINT_T0);
#endif                              
                            zmm8  = _mm512_load_ps(&xre[i+208]);
                            zmm9  = _mm512_load_ps(&yre[i+208]);
                            _mm512_store_ps(&zre[i+208], _mm512_sub_ps(zmm8,zmm9));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+224],_MM_HINT_T0);
#endif                              
                            zmm11 = _mm512_load_ps(&xre[i+224]);
                            zmm12 = _mm512_load_ps(&yre[i+224]);
                            _mm512_store_ps(&zre[i+224], _mm512_sub_ps(zmm11,zmm12));
#if (CSUB_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                            _mm_prefetch((const char *)&xre[i+240],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+240],_MM_HINT_T0);
#endif                              
                            zmm14  = _mm512_load_ps(&xre[i+240]);
                            zmm15  = _mm512_load_ps(&yre[i+240]);
                            _mm512_store_ps(&zre[i+240], _mm512_sub_ps(zmm14,zmm15));
                            
                        }

                       for(; (i+191) < n; i += 192)
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm12,zmm13));
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm15,zmm16));
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm18,zmm19));
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm21,zmm22));
                            zmm24  = _mm512_load_ps(&xre[i+128]);
                            zmm25  = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm25,zmm24));
                            zmm27  = _mm512_load_ps(&xre[i+144]);
                            zmm28  = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm27,zmm28));
                            zmm30  = _mm512_load_ps(&xre[i+160]);
                            zmm31  = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm30,zmm31));
                            zmm2  = _mm512_load_ps(&xre[i+176]);
                            zmm3  = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm2,zmm3));
                            
                       }

                       for(; (i+127) < n; i += 128)
                       {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
                            zmm12  = _mm512_load_ps(&xre[i+64]);
                            zmm13  = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm12,zmm13));
                            zmm15  = _mm512_load_ps(&xre[i+80]);
                            zmm16  = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm15,zmm16));
                            zmm18  = _mm512_load_ps(&xre[i+96]);
                            zmm19  = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm18,zmm19));
                            zmm21  = _mm512_load_ps(&xre[i+112]);
                            zmm22  = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm21,zmm22));
                            
                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            zmm6  = _mm512_load_ps(&xre[i+32]);
                            zmm7  = _mm512_load_ps(&yre[i+32]); //yre
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm6,zmm7));
                            zmm9  = _mm512_load_ps(&xre[i+48]);
                            zmm10 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm9,zmm10));
                            

                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            zmm3  = _mm512_load_ps(&xre[i+16]);
                            zmm4  = _mm512_load_ps(&yre[i+16]); //yre
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm3,zmm4));
                            

                       }

                      for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]); //yre
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                           

                      }

                      for(; (i+0) < n; i += 1) 
                      {
                            zre[i] = xre[i] - yre[i];
                      }
}


                  
                  