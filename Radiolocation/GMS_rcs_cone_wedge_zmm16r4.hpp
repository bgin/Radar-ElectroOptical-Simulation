

#ifndef __GMS_RCS_CONE_WEDGE_ZMM16R4_HPP__
#define __GMS_RCS_CONE_WEDGE_ZMM16R4_HPP__ 130320231230


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

    const unsigned int GMS_RCS_CONE_WEDGE_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_CONE_WEDGE_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_CONE_WEDGE_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_CONE_WEDGE_ZMM16R4_FULLVER =
      1000U*GMS_RCS_CONE_WEDGE_ZMM16R4_MAJOR+
      100U*GMS_RCS_CONE_WEDGE_ZMM16R4_MINOR+
      10U*GMS_RCS_CONE_WEDGE_ZMM16R4_MICRO;
    const char * const GMS_RCS_CONE_WEDGE_ZMM16R4_CREATION_DATE = "13-03-2023 12:30 PM +00200 (MON 13 MAR 2023 GMT+2)";
    const char * const GMS_RCS_CONE_WEDGE_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_CONE_WEDGE_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_CONE_WEDGE_ZMM16R4_DESCRIPTION   = "AVX512 optimized Cones and Wedges Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"



namespace  gms {

 
        namespace radiolocation {


                   /*
                       Small-angle cone (alpha ~ 0).
                       Backscattered RCS.
                       Formula 6.2-12
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6212_zmm16r4(const __m512 gam0,
                                            const __m512 alp) {

                          const __m512 _4pi  = _mm512_set1_ps(12.566370614359172953850573533118f);
                          const __m512 _1    = _mm512_set1_ps(1.0f);
                          const __m512 _3    = _mm512_set1_ps(3.0f);
                          register __m512 rcs,gam2,calp,trm1,trm2,trm3,x0;
                          gam2  = _mm512_mul_ps(gam0,gam0);
                          calp  = xcosf(alp);
                          trm1  = _mm512_div_ps(gam2,_4pi);
                          x0    = _mm512_sub_ps(_1,calp);
                          trm3  = _mm512_fmadd_ps(_3,x0,_1);
                          trm2  = _mm512_mul_ps(x0,x0);
                          rcs   = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6212_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const float * __restrict __ATTR_ALIGN__(64) palp) {

                          register __m512 gam0 = _mm512_load_ps(&pgam0[0]);
                          register __m512 alp  = _mm512_load_ps(&palp[0]);
                          const __m512 _4pi  = _mm512_set1_ps(12.566370614359172953850573533118f);
                          const __m512 _1    = _mm512_set1_ps(1.0f);
                          const __m512 _3    = _mm512_set1_ps(3.0f);
                          register __m512 rcs,gam2,calp,trm1,trm2,trm3,x0;
                          gam2  = _mm512_mul_ps(gam0,gam0);
                          calp  = xcosf(alp);
                          trm1  = _mm512_div_ps(gam2,_4pi);
                          x0    = _mm512_sub_ps(_1,calp);
                          trm3  = _mm512_fmadd_ps(_3,x0,_1);
                          trm2  = _mm512_mul_ps(x0,x0);
                          rcs   = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6212_zmm16r4_u(const float * __restrict  pgam0,
                                              const float * __restrict  palp) {

                          register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);
                          register __m512 alp  = _mm512_loadu_ps(&palp[0]);
                          const __m512 _4pi  = _mm512_set1_ps(12.566370614359172953850573533118f);
                          const __m512 _1    = _mm512_set1_ps(1.0f);
                          const __m512 _3    = _mm512_set1_ps(3.0f);
                          register __m512 rcs,gam2,calp,trm1,trm2,trm3,x0;
                          gam2  = _mm512_mul_ps(gam0,gam0);
                          calp  = xcosf(alp);
                          trm1  = _mm512_div_ps(gam2,_4pi);
                          x0    = _mm512_sub_ps(_1,calp);
                          trm3  = _mm512_fmadd_ps(_3,x0,_1);
                          trm2  = _mm512_mul_ps(x0,x0);
                          rcs   = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                  }


                     /*
                           Small-angle cone (alpha ~ pi/2).
                           Backscattered RCS.
                           Formula 6.2-13
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6213_zmm16r4(const __m512 gam0,
                                            const __m512 alp) {

                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                         
                          register __m512 rcs,calp,calp2,calp4,gam2,trm1,trm2,x0;
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          calp = xcosf(alp);
                          x0   = _mm512_div_ps(gam2,_16pi);
                          calp2= _mm512_mul_ps(calp,calp);
                          trm2 = _mm512_sub_ps(_1,_mm512_add_ps(calp2,calp2));
                          calp4= _mm512_mul_ps(calp2,calp2);
                          trm1 = _mm512_mul_ps(x0,calp4);
                          rcs  = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6213_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const float * __restrict __ATTR_ALIGN__(64) palp) {

                          register __m512 gam0 = _mm512_load_ps(&pgam0[0]);
                          register __m512 alp  = _mm512_load_ps(&palp[0]);

                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                         
                          register __m512 rcs,calp,calp2,calp4,gam2,trm1,trm2,x0;
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          calp = xcosf(alp);
                          x0   = _mm512_div_ps(gam2,_16pi);
                          calp2= _mm512_mul_ps(calp,calp);
                          trm2 = _mm512_sub_ps(_1,_mm512_add_ps(calp2,calp2));
                          calp4= _mm512_mul_ps(calp2,calp2);
                          trm1 = _mm512_mul_ps(x0,calp4);
                          rcs  = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                  }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6213_zmm16r4_u(const float * __restrict  pgam0,
                                              const float * __restrict  palp) {

                          register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);
                          register __m512 alp  = _mm512_loadu_ps(&palp[0]);

                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                         
                          register __m512 rcs,calp,calp2,calp4,gam2,trm1,trm2,x0;
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          calp = xcosf(alp);
                          x0   = _mm512_div_ps(gam2,_16pi);
                          calp2= _mm512_mul_ps(calp,calp);
                          trm2 = _mm512_sub_ps(_1,_mm512_add_ps(calp2,calp2));
                          calp4= _mm512_mul_ps(calp2,calp2);
                          trm1 = _mm512_mul_ps(x0,calp4);
                          rcs  = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                  }


                   /*
                         Backscattering case.
                         E-field scattered for (phi component).
                         Formula 6.2-16
    
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6216_zmm16r4(const __m512 k0,
                                           const __m512 r,
                                           const __m512 alp,
                                           const __m512 tht,
                                           __m512 * __restrict ESr,
                                           __m512 * __restrict ESi) {
  
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _3  = _mm512_set1_ps(3.0f);
                        const __m512 _4  = _mm512_set1_ps(4.0f);
                        register __m512 k0r,inv,alph,alphs,ctht,t0r,t0i;
                        register __m512 ctht2,num,den,ear,eai,cer,cei,rat;
                        k0r  = _mm512_mul_ps(k0,r);
                        ctht = xcosf(tht);
                        alph = _mm512_mul_ps(alp,hlf);
                        ctht2= _mm512_mul_ps(ctht,ctht);
                        ear  = _mm512_setzero_ps();
                        ctht3= _mm512_mul_ps(ctht2,ctht);
                        eai  = k0r;
                        num  = _mm512_add_ps(_3,ctht2);
                        inv  = _mm512_rcp14_ps(k0r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        den  = _mm512_mul_ps(_4,ctht3);
                        t0r  = _mm512_mul_ps(cer,inv);
                        rat  = _mm512_div_ps(num,den);
                        t0i  = _mm512_mul_ps(cei,inv);
                        alphs= _mm512_mul_ps(_mm512_mul_ps(alph,alph),rat);
                        *ESr = _mm512_mul_ps(t0r,alphs);
                        *ESi = _mm512_mul_ps(t0i,alphs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6216_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) palp,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64) ESi) {
  
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 alp = _mm512_load_ps(&palp[0]);
                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _3  = _mm512_set1_ps(3.0f);
                        const __m512 _4  = _mm512_set1_ps(4.0f);
                        register __m512 k0r,inv,alph,alphs,ctht,t0r,t0i;
                        register __m512 ctht2,num,den,ear,eai,cer,cei,rat;
                        k0r  = _mm512_mul_ps(k0,r);
                        ctht = xcosf(tht);
                        alph = _mm512_mul_ps(alp,hlf);
                        ctht2= _mm512_mul_ps(ctht,ctht);
                        ear  = _mm512_setzero_ps();
                        ctht3= _mm512_mul_ps(ctht2,ctht);
                        eai  = k0r;
                        num  = _mm512_add_ps(_3,ctht2);
                        inv  = _mm512_rcp14_ps(k0r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        den  = _mm512_mul_ps(_4,ctht3);
                        t0r  = _mm512_mul_ps(cer,inv);
                        rat  = _mm512_div_ps(num,den);
                        t0i  = _mm512_mul_ps(cei,inv);
                        alphs= _mm512_mul_ps(_mm512_mul_ps(alph,alph),rat);
                        _mm512_store_ps(&ESr[0] ,_mm512_mul_ps(t0r,alphs));
                        _mm512_store_ps(&ESi[0] ,_mm512_mul_ps(t0i,alphs));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6216_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict  pr,
                                           const float * __restrict  palp,
                                           const float * __restrict  ptht,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {
  
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 alp = _mm512_loadu_ps(&palp[0]);
                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _3  = _mm512_set1_ps(3.0f);
                        const __m512 _4  = _mm512_set1_ps(4.0f);
                        register __m512 k0r,inv,alph,alphs,ctht,t0r,t0i;
                        register __m512 ctht2,num,den,ear,eai,cer,cei,rat;
                        k0r  = _mm512_mul_ps(k0,r);
                        ctht = xcosf(tht);
                        alph = _mm512_mul_ps(alp,hlf);
                        ctht2= _mm512_mul_ps(ctht,ctht);
                        ear  = _mm512_setzero_ps();
                        ctht3= _mm512_mul_ps(ctht2,ctht);
                        eai  = k0r;
                        num  = _mm512_add_ps(_3,ctht2);
                        inv  = _mm512_rcp14_ps(k0r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        den  = _mm512_mul_ps(_4,ctht3);
                        t0r  = _mm512_mul_ps(cer,inv);
                        rat  = _mm512_div_ps(num,den);
                        t0i  = _mm512_mul_ps(cei,inv);
                        alphs= _mm512_mul_ps(_mm512_mul_ps(alph,alph),rat);
                        _mm512_storeu_ps(&ESr[0] ,_mm512_mul_ps(t0r,alphs));
                        _mm512_storeu_ps(&ESi[0] ,_mm512_mul_ps(t0i,alphs));
                 }

          }


}

















#endif  /*__GMS_RCS_CONE_WEDGE_ZMM16R4_HPP__*/
