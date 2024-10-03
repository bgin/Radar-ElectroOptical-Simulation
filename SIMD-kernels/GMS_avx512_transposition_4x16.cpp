



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







#include "GMS_avx512_transposition_4x16.h"


 
                     /*
                            In-place.
                            Data pointer must point to an array (__m512) of 4 elements.
                       */

                                        
		     
		     
		      void gms::math::transpose_zmm16r4_4x16(__m512 * __restrict inout) {

                           _mm_prefetch((const char*)&inout[0],_MM_HINT_T0);
                           register __m512 x0 = _mm512_unpacklo_ps(inout[0],inout[1]);
                           register __m512 x1 = _mm512_unpackhi_ps(inout[0],inout[1]);
                           _mm_prefetch((const char*)&inout[2],_MM_HINT_T0);
                           register __m512 x2 = _mm512_unpacklo_ps(inout[2],inout[3]);
                           register __m512 x3 = _mm512_unpackhi_ps(inout[2],inout[3]);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           inout[0]        = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0));
                           inout[1]        = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0));
                           inout[2]        = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1));
                           inout[3]        = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1));
                    }

                     /*
                            Out-of-place.
                            Pointers must point to an arrays (__m512) of 4 elements.
                       */

                     
		     
		      void gms::math::transpose_zmm16r4_4x16(const __m512 * __restrict in,
                                                  __m512 * __restrict out) {

                           _mm_prefetch((const char*)&in[0],_MM_HINT_T0);
                           register __m512 x0 = _mm512_unpacklo_ps(in[0],in[1]);
                           register __m512 x1 = _mm512_unpackhi_ps(in[0],in[1]);
                           _mm_prefetch((const char*)&in[2],_MM_HINT_T0);
                           register __m512 x2 = _mm512_unpacklo_ps(in[2],in[3]);
                           register __m512 x3 = _mm512_unpackhi_ps(in[2],in[3]);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           out[0]         = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0));
                           out[1]         = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0));
                           out[2]         = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1));
                           out[3]         = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1));
                    }


		    
		      void gms::math::transpose_zmm16r4_4x16_u(float * __restrict inout) {

                           _mm_prefetch((const char*)&inout[0*16],_MM_HINT_T0);
                           register __m512 zmm0 = _mm512_loadu_ps(&inout[0*16]);
                           register __m512 zmm1 = _mm512_loadu_ps(&inout[1*16]);
                           register __m512   x0 = _mm512_unpacklo_ps(zmm0,zmm1);
                           register __m512   x1 = _mm512_unpackhi_ps(zmm0,zmm1);
                           _mm_prefetch((const char*)&inout[2*16],_MM_HINT_T0);
                           register __m512 zmm2 = _mm512_loadu_ps(&inout[2*16]);
                           register __m512 zmm3 = _mm512_loadu_ps(&inout[3*16]);
                           register __m512 x2 = _mm512_unpacklo_ps(zmm2,zmm3);
                           register __m512 x3 = _mm512_unpackhi_ps(zmm2,zmm3);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           _mm512_storeu_ps(&inout[0*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_storeu_ps(&inout[1*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_storeu_ps(&inout[2*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                           _mm512_storeu_ps(&inout[3*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                    }


                   
		      
		      void gms::math::transpose_zmm16r4_4x16_a(float * __restrict __ATTR_ALIGN__(64) inout) {

                           _mm_prefetch((const char*)&inout[0*16],_MM_HINT_T0);
                           register __m512 zmm0 = _mm512_load_ps(&inout[0*16]);
                           register __m512 zmm1 = _mm512_load_ps(&inout[1*16]);
                           register __m512   x0 = _mm512_unpacklo_ps(zmm0,zmm1);
                           register __m512   x1 = _mm512_unpackhi_ps(zmm0,zmm1);
                           _mm_prefetch((const char*)&inout[2*16],_MM_HINT_T0);
                           register __m512 zmm2 = _mm512_load_ps(&inout[2*16]);
                           register __m512 zmm3 = _mm512_load_ps(&inout[3*16]);
                           register __m512 x2 = _mm512_unpacklo_ps(zmm2,zmm3);
                           register __m512 x3 = _mm512_unpackhi_ps(zmm2,zmm3);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           _mm512_store_ps(&inout[0*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_store_ps(&inout[1*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_store_ps(&inout[2*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                           _mm512_store_ps(&inout[3*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                    }


                     
		      
		      void gms::math::transpose_zmm16r4_4x16_u(float * __restrict in,
                                                    float * __restrict out) {

                           _mm_prefetch((const char*)&in[0*16],_MM_HINT_T0);
                           register __m512 zmm0 = _mm512_loadu_ps(&in[0*16]);
                           register __m512 zmm1 = _mm512_loadu_ps(&in[1*16]);
                           register __m512   x0 = _mm512_unpacklo_ps(zmm0,zmm1);
                           register __m512   x1 = _mm512_unpackhi_ps(zmm0,zmm1);
                           _mm_prefetch((const char*)&in[2*16],_MM_HINT_T0);
                           register __m512 zmm2 = _mm512_loadu_ps(&in[2*16]);
                           register __m512 zmm3 = _mm512_loadu_ps(&in[3*16]);
                           register __m512 x2 = _mm512_unpacklo_ps(zmm2,zmm3);
                           register __m512 x3 = _mm512_unpackhi_ps(zmm2,zmm3);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           _mm512_storeu_ps(&out[0*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_storeu_ps(&out[1*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_storeu_ps(&out[2*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                           _mm512_storeu_ps(&out[3*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                    }


                      
                     
		     
		      void gms::math::transpose_zmm16r4_4x16_a(float * __restrict __ATTR_ALIGN__(64) in,
                                                    float * __restrict __ATTR_ALIGN__(64) out) {

                           _mm_prefetch((const char*)&in[0*16],_MM_HINT_T0);
                           register __m512 zmm0 = _mm512_load_ps(&in[0*16]);
                           register __m512 zmm1 = _mm512_load_ps(&in[1*16]);
                           register __m512   x0 = _mm512_unpacklo_ps(zmm0,zmm1);
                           register __m512   x1 = _mm512_unpackhi_ps(zmm0,zmm1);
                           _mm_prefetch((const char*)&in[2*16],_MM_HINT_T0);
                           register __m512 zmm2 = _mm512_load_ps(&in[2*16]);
                           register __m512 zmm3 = _mm512_load_ps(&in[3*16]);
                           register __m512 x2 = _mm512_unpacklo_ps(zmm2,zmm3);
                           register __m512 x3 = _mm512_unpackhi_ps(zmm2,zmm3);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           _mm512_store_ps(&out[0*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_store_ps(&out[1*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                           _mm512_store_ps(&out[2*16] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                           _mm512_store_ps(&out[3*16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                    }


                    
		   
		      void gms::math::transpose_zmm16r4_4x16(const __m512 in0,
                                                  const __m512 in1,
                                                  const __m512 in2,
                                                  const __m512 in3,
                                                  __m512 * __restrict out0,
                                                  __m512 * __restrict out1,
                                                  __m512 * __restrict out2,
                                                  __m512 * __restrict out3) {

                           register __m512   x0 = _mm512_unpacklo_ps(in0,in1);
                           register __m512   x1 = _mm512_unpackhi_ps(in0,in1);
                           register __m512 x2 = _mm512_unpacklo_ps(in2,in3);
                           register __m512 x3 = _mm512_unpackhi_ps(in2,in3);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           *out0           = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0));
                           *out1           = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0));
                           *out2           = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1));
                           *out3           = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1));
                  }


                    

#include <cstdint>
#include "GMS_common.h"

                     
		     
		      void gms::math::transpose_zmm16r4_4x16_u8x_u(float * __restrict in,
                                                        float * __restrict out,
                                                        const int32_t n,
                                                        const int32_t m) {

                           if(__builtin_expect(0<=n,0) || 
                              __builtin_expect(0<=m,0)) { return;} 
                             
                           register __m512 zmm0,zmm1,zmm2,zmm3;
                           register __m512 x0,x1,x2,x3;
                           register __m512 y0,y1,y2,y3;
                           int32_t i,j;
                           
                           for(i = 0; (i+511) < n; i += 512) {
                               _mm_prefetch((const char*)&in[i+0],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+64],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+128],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+128]);
                               zmm1 = _mm512_loadu_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+160]);
                               zmm3 = _mm512_loadu_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+192],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+192]);
                               zmm1 = _mm512_loadu_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+224]);
                               zmm3 = _mm512_loadu_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+240],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+256]);
                               zmm1 = _mm512_loadu_ps(&in[i+272]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+288]);
                               zmm3 = _mm512_loadu_ps(&in[i+304]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+256] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+272] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+288] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+304] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+304],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+320]);
                               zmm1 = _mm512_loadu_ps(&in[i+336]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+352]);
                               zmm3 = _mm512_loadu_ps(&in[i+368]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+320] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+336] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+352] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+368] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+368],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+384]);
                               zmm1 = _mm512_loadu_ps(&in[i+400]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+416]);
                               zmm3 = _mm512_loadu_ps(&in[i+432]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+384] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+400] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+416] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+432] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+432],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+448]);
                               zmm1 = _mm512_loadu_ps(&in[i+464]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+480]);
                               zmm3 = _mm512_loadu_ps(&in[i+496]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+448] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+464] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+480] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+496] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                           }

                            for(; (i+447) < n; i += 448) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+128]);
                               zmm1 = _mm512_loadu_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+160]);
                               zmm3 = _mm512_loadu_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+192]);
                               zmm1 = _mm512_loadu_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+224]);
                               zmm3 = _mm512_loadu_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+256]);
                               zmm1 = _mm512_loadu_ps(&in[i+272]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+288]);
                               zmm3 = _mm512_loadu_ps(&in[i+304]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+256] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+272] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+288] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+304] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+320]);
                               zmm1 = _mm512_loadu_ps(&in[i+336]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+352]);
                               zmm3 = _mm512_loadu_ps(&in[i+368]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+320] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+336] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+352] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+368] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+384]);
                               zmm1 = _mm512_loadu_ps(&in[i+400]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+416]);
                               zmm3 = _mm512_loadu_ps(&in[i+432]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+384] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+400] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+416] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+432] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                           }

                            for(; (i+383) < n; i += 384) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+128]);
                               zmm1 = _mm512_loadu_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+160]);
                               zmm3 = _mm512_loadu_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+192]);
                               zmm1 = _mm512_loadu_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+224]);
                               zmm3 = _mm512_loadu_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+256]);
                               zmm1 = _mm512_loadu_ps(&in[i+272]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+288]);
                               zmm3 = _mm512_loadu_ps(&in[i+304]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+256] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+272] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+288] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+304] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+320]);
                               zmm1 = _mm512_loadu_ps(&in[i+336]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+352]);
                               zmm3 = _mm512_loadu_ps(&in[i+368]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+320] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+336] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+352] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+368] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                            for(; (i+255) < n; i += 256) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+128]);
                               zmm1 = _mm512_loadu_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+160]);
                               zmm3 = _mm512_loadu_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+192]);
                               zmm1 = _mm512_loadu_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+224]);
                               zmm3 = _mm512_loadu_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                             for(; (i+127) < n; i += 128) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                            for(; (i+63) < n; i += 64) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                            }

                          for(; (i+0) < n; i += 1) {
                              for(j = 0; j != m; ++j) {
                                  out[Ix2D(j,m,i)] = in[Ix2D(i,m,j)];
                              }
                          }
                    }


		      
		      void gms::math::transpose_zmm16r4_4x16_u4x_u(float * __restrict in,
                                                        float * __restrict out,
                                                        const int32_t n,
                                                        const int32_t m) {

                           if(__builtin_expect(0<=n,0) || 
                              __builtin_expect(0<=m,0)) { return;} 
                             
                           register __m512 zmm0,zmm1,zmm2,zmm3;
                           register __m512 x0,x1,x2,x3;
                           register __m512 y0,y1,y2,y3;
                           int32_t i,j;

                           for(i = 0; (i+255) < n; i += 256) {
                               _mm_prefetch((const char*)&in[i+0],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+64],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+128],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+128]);
                               zmm1 = _mm512_loadu_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+160]);
                               zmm3 = _mm512_loadu_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+192],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&in[i+192]);
                               zmm1 = _mm512_loadu_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+224]);
                               zmm3 = _mm512_loadu_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                            for(; (i+191) < n; i += 192) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+128]);
                               zmm1 = _mm512_loadu_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+160]);
                               zmm3 = _mm512_loadu_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                           }

                            for(; (i+127) < n; i += 128) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_loadu_ps(&in[i+64]);
                               zmm1 = _mm512_loadu_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1))); 
                          }

                            for(; (i+63) < n; i += 64) {
                               zmm0 = _mm512_loadu_ps(&in[i+0]);
                               zmm1 = _mm512_loadu_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+32]);
                               zmm3 = _mm512_loadu_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_storeu_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_storeu_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_storeu_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                         }

                             for(; (i+0) < n; i += 1) {
                                for(j = 0; j != m; ++j) {
                                     out[Ix2D(j,m,i)] = in[Ix2D(i,m,j)];
                              }
                          }

                    }


                    
		      
		      void gms::math::transpose_zmm16r4_4x16_u8x_a(float * __restrict __ATTR_ALIGN__(64) in,
                                                        float * __restrict __ATTR_ALIGN__(64) out,
                                                        const int32_t n,
                                                        const int32_t m) {

                           if(__builtin_expect(0<=n,0) || 
                              __builtin_expect(0<=m,0)) { return;} 
                             
                           register __m512 zmm0,zmm1,zmm2,zmm3;
                           register __m512 x0,x1,x2,x3;
                           register __m512 y0,y1,y2,y3;
                           int32_t i,j;
                           
                           for(i = 0; (i+511) < n; i += 512) {
                               _mm_prefetch((const char*)&in[i+0],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+64],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+96]);
                               zmm3 = _mm512_load_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+128],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+128]);
                               zmm1 = _mm512_load_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+160]);
                               zmm3 = _mm512_load_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+192],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+192]);
                               zmm1 = _mm512_load_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+224]);
                               zmm3 = _mm512_load_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+240],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+256]);
                               zmm1 = _mm512_load_ps(&in[i+272]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+288]);
                               zmm3 = _mm512_load_ps(&in[i+304]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+256] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+272] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+288] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+304] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+304],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+320]);
                               zmm1 = _mm512_load_ps(&in[i+336]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+352]);
                               zmm3 = _mm512_load_ps(&in[i+368]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+320] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+336] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+352] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+368] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+368],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+384]);
                               zmm1 = _mm512_load_ps(&in[i+400]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+416]);
                               zmm3 = _mm512_load_ps(&in[i+432]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+384] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+400] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+416] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+432] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+432],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+448]);
                               zmm1 = _mm512_load_ps(&in[i+464]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+480]);
                               zmm3 = _mm512_load_ps(&in[i+496]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+448] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+464] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+480] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+496] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                           }

                            for(; (i+447) < n; i += 448) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+96]);
                               zmm3 = _mm512_load_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+128]);
                               zmm1 = _mm512_load_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+160]);
                               zmm3 = _mm512_load_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+192]);
                               zmm1 = _mm512_load_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+224]);
                               zmm3 = _mm512_load_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+256]);
                               zmm1 = _mm512_load_ps(&in[i+272]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+288]);
                               zmm3 = _mm512_load_ps(&in[i+304]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+256] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+272] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+288] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+304] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+320]);
                               zmm1 = _mm512_load_ps(&in[i+336]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+352]);
                               zmm3 = _mm512_load_ps(&in[i+368]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+320] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+336] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+352] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+368] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+384]);
                               zmm1 = _mm512_load_ps(&in[i+400]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+416]);
                               zmm3 = _mm512_load_ps(&in[i+432]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+384] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+400] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+416] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+432] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                           }

                            for(; (i+383) < n; i += 384) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+96]);
                               zmm3 = _mm512_load_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+128]);
                               zmm1 = _mm512_load_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+160]);
                               zmm3 = _mm512_load_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+192]);
                               zmm1 = _mm512_load_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+224]);
                               zmm3 = _mm512_load_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+256]);
                               zmm1 = _mm512_load_ps(&in[i+272]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+288]);
                               zmm3 = _mm512_load_ps(&in[i+304]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+256] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+272] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+288] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+304] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+320]);
                               zmm1 = _mm512_load_ps(&in[i+336]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+352]);
                               zmm3 = _mm512_load_ps(&in[i+368]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+320] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+336] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+352] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+368] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                            for(; (i+255) < n; i += 256) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+96]);
                               zmm3 = _mm512_load_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+128]);
                               zmm1 = _mm512_load_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+160]);
                               zmm3 = _mm512_load_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+192]);
                               zmm1 = _mm512_load_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+224]);
                               zmm3 = _mm512_load_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                             for(; (i+127) < n; i += 128) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+96]);
                               zmm3 = _mm512_load_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                            for(; (i+63) < n; i += 64) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                            }

                          for(; (i+0) < n; i += 1) {
                              for(j = 0; j != m; ++j) {
                                  out[Ix2D(j,m,i)] = in[Ix2D(i,m,j)];
                              }
                          }
                    }


                    
		     
		      void gms::math::transpose_zmm16r4_4x16_u4x_a(float * __restrict __ATTR_ALIGN__(64) in,
                                                        float * __restrict __ATTR_ALIGN__(64) out,
                                                        const int32_t n,
                                                        const int32_t m) {

                           if(__builtin_expect(0<=n,0) || 
                              __builtin_expect(0<=m,0)) { return;} 
                             
                           register __m512 zmm0,zmm1,zmm2,zmm3;
                           register __m512 x0,x1,x2,x3;
                           register __m512 y0,y1,y2,y3;
                           int32_t i,j;

                           for(i = 0; (i+255) < n; i += 256) {
                               _mm_prefetch((const char*)&in[i+0],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+64],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+96]);
                               zmm3 = _mm512_load_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+128],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+128]);
                               zmm1 = _mm512_load_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+160]);
                               zmm3 = _mm512_load_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               _mm_prefetch((const char*)&in[i+192],_MM_HINT_T0);
                               zmm0 = _mm512_load_ps(&in[i+192]);
                               zmm1 = _mm512_load_ps(&in[i+208]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+224]);
                               zmm3 = _mm512_load_ps(&in[i+240]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+192] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+208] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+224] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+240] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                          }

                            for(; (i+191) < n; i += 192) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_loadu_ps(&in[i+96]);
                               zmm3 = _mm512_loadu_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+128]);
                               zmm1 = _mm512_load_ps(&in[i+144]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+160]);
                               zmm3 = _mm512_load_ps(&in[i+176]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+128] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+144] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+160] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+176] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                           }

                            for(; (i+127) < n; i += 128) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                               zmm0 = _mm512_load_ps(&in[i+64]);
                               zmm1 = _mm512_load_ps(&in[i+80]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+96]);
                               zmm3 = _mm512_load_ps(&in[i+112]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+64] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+80] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+96] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+112] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1))); 
                          }

                            for(; (i+63) < n; i += 64) {
                               zmm0 = _mm512_load_ps(&in[i+0]);
                               zmm1 = _mm512_load_ps(&in[i+16]);
                               x0   = _mm512_unpacklo_ps(zmm0,zmm1);
                               x1   = _mm512_unpackhi_ps(zmm0,zmm1);
                               zmm2 = _mm512_load_ps(&in[i+32]);
                               zmm3 = _mm512_load_ps(&in[i+48]);
                               x2   = _mm512_unpacklo_ps(zmm2,zmm3);
                               x3   = _mm512_unpackhi_ps(zmm2,zmm3);
                               y0   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                               y1   = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                               y2   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                               y3   = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                               x0   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                               x1   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                               x2   = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                               x3   = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                               _mm512_store_ps(&out[i+0] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+16] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0)));
                               _mm512_store_ps(&out[i+32] ,_mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1)));
                               _mm512_store_ps(&out[i+48] ,_mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1)));
                         }

                             for(; (i+0) < n; i += 1) {
                                for(j = 0; j != m; ++j) {
                                     out[Ix2D(j,m,i)] = in[Ix2D(i,m,j)];
                              }
                          }

                    }

                     
                                                      

