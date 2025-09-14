




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

#include "GMS_phase_pdf_sse.h"





                   __m128 gms::math::phase_pdf_xmm4r4(const __m128 a,        // signal-to-noise ratio
                                           const __m128 phi) {    // vector of phi (rad) arguments

                         const __m128 inv2pi = _mm_set1_ps(0.159154943091895335768883763373f);
                         const __m128 sq2pi  = _mm_set1_ps(2.506628274631000502415765284811f);
                         const __m128 _1     = _mm_set1_ps(1.0f);
                         const __m128 _0_5   = _mm_set1_ps(0.5f);
                         __m128 acphi,errf,expt2,hphi;
                         __m128 aa,haa,invex,t0,t1,t2,c2phi;
                         __m128 f1;
                         aa   = _mm_mul_ps(a,a);
                         hphi = _mm_mul_ps(_0_5,phi);
                         haa  = _mm_mul_ps(_0_5,aa);
                         acphi= _mm_mul_ps(a,_mm_cos_ps(phi));
                         t0   = _mm_mul_ps(sq2pi,acphi);
                         invex= _mm_div_ps(_1,_mm_exp_ps(haa));
                         t1   = _mm_fmadd_ps(sq2pi,acphi,_1);
                         errf = _mm_erf_ps(acphi);
                         t2   = _mm_cos_ps(hphi);
                         c2phi= _mm_mul_ps(t2,t2);
                         expt2= _mm_exp_ps(_mm_mul_ps(aa,c2phi));
                         t2   = _mm_mul_ps(t1,_mm_mul_ps(errf,expt2));
                         f1   = _mm_mul_ps(inv2pi,_mm_mul_ps(invex,t2));
                         return (f1);
                 }


                 
                   __m128d gms::math::phase_pdf_xmm2r8(const __m128d a,        // signal-to-noise ratio
                                            const __m128d phi) {    // vector of phi (rad) arguments

                         const __m128d inv2pi = _mm_set1_pd(0.159154943091895335768883763373);
                         const __m128d sq2pi  = _mm_set1_pd(2.506628274631000502415765284811);
                         const __m128d _1     = _mm_set1_pd(1.0f);
                         const __m128d _0_5   = _mm_set1_pd(0.5f);
                         __m128d acphi,errf,expt2,hphi;
                         __m128d aa,haa,invex,t0,t1,t2,c2phi;
                         __m128d f1;
                         aa   = _mm_mul_pd(a,a);
                         hphi = _mm_mul_pd(_0_5,phi);
                         haa  = _mm_mul_pd(_0_5,aa);
                         acphi= _mm_mul_pd(a,_mm_cos_pd(phi));
                         t0   = _mm_mul_pd(sq2pi,acphi);
                         invex= _mm_div_pd(_1,_mm_exp_pd(haa));
                         t1   = _mm_fmadd_pd(sq2pi,acphi,_1);
                         errf = _mm_erf_pd(acphi);
                         t2   = _mm_cos_pd(hphi);
                         c2phi= _mm_mul_pd(t2,t2);
                         expt2= _mm_exp_pd(_mm_mul_pd(aa,c2phi));
                         t2   = _mm_mul_pd(t1,_mm_mul_pd(errf,expt2));
                         f1   = _mm_mul_pd(inv2pi,_mm_mul_pd(invex,t2));
                         return (f1);
                 }


            
      




















