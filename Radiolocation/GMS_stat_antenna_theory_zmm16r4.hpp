

#ifndef __GMS_STAT_ANTENNA_THEORY_ZMM16R4_HPP__
#define __GMS_STAT_ANTENNA_THEORY_ZMM16R4_HPP__ 060620231414


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

namespace file_version {

    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_FULLVER =
      1000U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MAJOR+
      100U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MINOR+
      10U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MICRO;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_CREATION_DATE = "04-06-2023 14:14 PM +00200 (SUN 04 JUN 2023 GMT+2)";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_DESCRIPTION   = "AVX512 (single) optimized statistical antenna theory kernels.";

}


#include <immintrin.h>
#include <complex>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"
#include "GMS_cspint_quad.hpp"
#include "GMS_avint_quad.hpp"
#include "GMS_cephes.h"

namespace  gms {


      namespace radiolocation {
      
      
                    /* The whole content is based on:
                       Шифрин Я.С. - "Вопросы статистической теории антенн"
                    */
                    
                    
               /*
                   Work (input) arrays for integrator 'cspint'
                   Multi-threaded i.e. (two-threaded)
               */
               __ATTR_ALIGN__(64) struct CSPINT_DATA_2T {
               
                       float * __restrict  Ya1; 
                       float * __restrict  Ya2; 
                       float * __restrict  Ya3; 
                       float * __restrict  Ea;  
                       float * __restrict  WRKa; 
                       float * __restrict  Yb1;
                       float * __restrict  Yb2; 
                       float * __restrict  Yb3; 
                       float * __restrict  Eb; 
                       float * __restrict  WRKb;  
               };
               
                /*
                   Work (input) arrays for integrator 'cspint'
                   Single-threaded.
               */
                __ATTR_ALIGN__(64) struct CSPINT_DATA_1T {
               
                       float * __restrict  Y1; 
                       float * __restrict  Y2; 
                       float * __restrict  Y3; 
                       float * __restrict  E;  
                       float * __restrict  WRK; 
                       
               };
               
               
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint single data vector ZMM.
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   std::complex<float> 
                   f11_zmm16r4_cspint_16e_a(const float * __restrict __ATTR_ALIGN__(64) pAz,
                                            const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const float * __restrict __ATTR_ALIGN__(64) pz,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2) {
                                            
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         register __m512 Az    = _mm512_load_ps(&pAz[0]);
                         register __m512 phiz  = _mm512_load_pd(&pphiz[0]);
                         register __m512 z     = _mm512_load_pd(&pz[0]);
                         register __m512 stht,ear,eai,cer,cei;
                         register __m512 k,vtht;
                         std::complex<float> fth;
                         register float sumr,sumi;
                         vtht= _mm512_set1_ps(tht);
                         k   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                         stht= xsinf(tht);
                         ear = _mm512_setzero_ps();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sumr= 0.0f;
                         cer = _mm512_mul_ps(Az,cer);
                         _mm512_store_ps(&intr[0],cer);
                         sumi= 0.0f;
                         cei = _mm512_mul_ps(Az,cei);
                         _mm512_store_ps(&inti[0],cei);
                         cspint(NTAB,&pz[0],&intr[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,&pz[0],&inti[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   std::complex<float> 
                   f11_zmm16r4_cspint_16e_u(const float * __restrict  pAz,
                                            const float * __restrict  pphiz,
                                            const float * __restrict  pz,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2) {
                                            
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16]  = {}; 
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         register __m512 Az    = _mm512_loadu_ps(&pAz[0]);
                         register __m512 phiz  = _mm512_loadu_pd(&pphiz[0]);
                         register __m512 z     = _mm512_loadu_pd(&pz[0]);
                         register __m512 stht,ear,eai,cer,cei;
                         register __m512 k,vtht;
                         std::complex<float> fth;
                         register float sumr,sumi;
                         vtht= _mm512_set1_ps(tht);
                         k   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                         stht= xsinf(tht);
                         ear = _mm512_setzero_ps();
                         eai = _mm512_fmadd_pd(stht,
                                           _mm512_mul_pd(k,z),phiz);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sumr= 0.0f;
                         cer = _mm512_mul_ps(Az,cer);
                         _mm512_storeu_ps(&intr[0],cer);
                         sumi= 0.0f;
                         cei = _mm512_mul_ps(Az,cei);
                         _mm512_storeu_ps(&inti[0],cei);
                         cspint(NTAB,&pz[0],&intr[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,&pz[0],&inti[0],-L2,L2,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         fth = {sumr,sumi};
                         return (fth);
                 }
                 
                 
                   /*
                        Formula 1.1, p. 14
                        Integrator cspint n-elements data arrays.
                        Single-threaded.
                   */
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   std::complex<float>
                   f11_zmm16r4_cspint_ne_a( const float * __restrict __ATTR_ALIGN__(64) pAz,
                                            const float * __restrict __ATTR_ALIGN__(64) pphiz,
                                            const float * __restrict __ATTR_ALIGN__(64) pz,
                                            float * __restrict __ATTR_ALIGN__(64) intr,
                                            float * __restrict __ATTR_ALIGN__(64) inti,
                                            CSPINT_DATA_1T & csd,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2,
                                            const int32_t NTAB) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<float> fth = {0.0f,0.0f};
                           fth = f11_zmm16r4_cspint_16_a(pAz,pphiz,pz,tht,gam,L2);
                           return (fth);
                        }     
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vtht,stht,vk;
                        register __m512 ear,cer,cei;
                        register std::complex<float> fth;
                        register float sumr,sumi; 
                        register float sint,k;
                        int32_t i; 
                        vtht = _mm512_set1_ps(tht);
                        sint = cephes_sinf(tht);
                        vk   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                        ear  = _mm512_setzero_ps();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsinf(vtht);
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                             register __m512 t0 = _mm512_load_ps(&pz[i]);
                             register __m512 t1 = _mm512_load_ps(&pphiz[i]);
                             register __m512 eai= _mm512_fmadd_ps(stht,
                                                            _mm512_mul_ps(vk,t0),t1);
                             register __m512 t2 = _mm512_load_ps(&pAz[i]);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             cer = _mm512_mul_ps(cer,t2);
                             _mm512_store_ps(&intr[i],cer);
                             cei = _mm512_mul_ps(cei,t2);
                             _mm512_store_ps(&inti[i],cei);
                        }   
                        sumr = 0.0f;
                        sumi = 0.0f;
                        for(; i != NTAB; ++i) {
                             register float t0 = pAz[i];
                             register float t1 = pphiz[i];
                             register float t2 = pz[i];
                             register float zz = k*t0;
                             register float eai= sint*zz+t1;
                             register std::complex<float> c = std::exp({0.0f,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        cspint(NTAB,&pz[0],&intr[0],-L2,L2,&csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumr); 
                        cspint(NTAB,&pz[0],&inti[0],-L2,L2, &csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumi);  
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   std::complex<float>
                   f11_zmm16r4_cspint_ne_u( const float * __restrict  pAz,
                                            const float * __restrict  pphiz,
                                            const float * __restrict  pz,
                                            float * __restrict  intr,
                                            float * __restrict  inti,
                                            CSPINT_DATA_1T & csd,
                                            const float tht,
                                            const float gam,
                                            const int32_t L2,
                                            const int32_t NTAB) {
                                            
                        if(__builtin_expect(NTAB==16,0)) {
                           std::complex<float> fth = {0.0f,0.0f};
                           fth = f11_zmm16r4_cspint_16_u(pAz,pphiz,pz,tht,gam,L2);
                           return (fth);
                        }     
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vtht,stht,vk;
                        register __m512 ear,cer,cei;
                        register std::complex<float> fth;
                        register float sumr,sumi; 
                        register float sint,k;
                        int32_t i; 
                        vtht = _mm512_set1_ps(tht);
                        sint = cephes_sinf(tht);
                        vk   = _mm512_mul_ps(C314159265358979323846264338328,
                                                          _mm512_set1_ps(gam));
                        ear  = _mm512_setzero_ps();
                        k    = C314159265358979323846264338328*gam;
                        stht = xsinf(vtht);
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pAz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pphiz[i],_MM_HINT_T0);
                             _mm_prefetch((const char*)&pz[i],_MM_HINT_T0);
                             register __m512 t0 = _mm512_loadu_ps(&pz[i]);
                             register __m512 t1 = _mm512_loadu_ps(&pphiz[i]);
                             register __m512 eai= _mm512_fmadd_ps(stht,
                                                            _mm512_mul_ps(vk,t0),t1);
                             register __m512 t2 = _mm512_loadu_ps(&pAz[i]);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             cer = _mm512_mul_ps(cer,t2);
                             _mm512_storeu_ps(&intr[i],cer);
                             cei = _mm512_mul_ps(cei,t2);
                             _mm512_storeu_ps(&inti[i],cei);
                        }   
                        sumr = 0.0f;
                        sumi = 0.0f;
                        for(; i != NTAB; ++i) {
                             register float t0 = pAz[i];
                             register float t1 = pphiz[i];
                             register float t2 = pz[i];
                             register float zz = k*t0;
                             register float eai= sint*zz+t1;
                             register std::complex<float> c = std::exp({0.0f,eai});
                             intr[i] = c.real()*t0;
                             inti[i] = c.imag()*t0;
                        } 
                        cspint(NTAB,&pz[0],&intr[0],-L2,L2,&csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumr); 
                        cspint(NTAB,&pz[0],&inti[0],-L2,L2, &csd.Y1[0],
                               &csd.Y2[0],&csd.Y3[0],&csd.E[0],&csd.WRK[0],sumi);  
                        fth = {sumr,sumi};
                        return (fth);            
                }
                
                
                /*
                    Formula 1.3, p. 15
                    The level of an amplitude.
                */
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 f13_zmm16r4(const __m512 A0,
                                      const __m512 A) {
                                      
                          register __m512 Bx,rat;
                          rat = _mm512_div_ps(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 f13_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pA0,
                                        const float * __restrict __ATTR_ALIGN__(64) pA) {
                              
                          register __m512 A0 = _mm512_load_ps(&pA0[0]);
                          register __m512 A  = _mm512_load_ps(&pA[0]);        
                          register __m512 Bx,rat;
                          rat = _mm512_div_ps(A,A0);
                          Bx  = xlogf(rat);
                          return (Bx);                    
                 }
                 
                 
                 
      
      } // radiolocation


} // radiolocation


















#endif /*__GMS_STAT_ANTENNA_THEORY_ZMM16R4_HPP__*/
