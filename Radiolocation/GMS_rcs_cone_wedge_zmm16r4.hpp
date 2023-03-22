

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


                   /*
                         Bistatic RCS case.
                         E-field scattered for (theta component).
                         Formula 6.2-14
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f6214_zmm16r4(const __m512 k0,
                                           const __m512 r,
                                           const __m512 alp,
                                           const __m512 tht1, //inc
                                           const __m512 tht2  //scat
                                           const __m512 phi1, //inc
                                           const __m512 phi2, //scat
                                           __m512 * __restrict ESr,
                                           __m512 * __restrict ESi) {

                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,htht1,htht2,phid,ear,eai,cer,cei;
                        register __m512 sphid,t0r,t0i,num1,den1,num2,den2,inv;
                        register __m512 ctht1,ctht2,sect1,sect2,x0,x1,alp2;
                        register __m512 stht1,stht2,chtht1,chtht2,rat1,rat2,cx0,cx1;
                        phid   = _mm512_sub_ps(phi1,phi2);
                        k0r    = _mm512_mul_ps(k0,r);
                        htht1  = _mm512_mul_ps(tht1,hlf);
                        sphid  = xsinf(phid);
                        htht2  = _mm512_mul_ps(tht2,hlf);
                        inv    = _mm512_rcp14_ps(k0r);
                        ear    = _mm512_setzero_ps();
                        ctht1  = xcosf(tht1);
                        cx1    = _mm512_mul_ps(ctht1);
                        x0     = _mm512_mul_ps(alp,hlf);
                        ctht2  = xcosf(tht2);
                        cx2    = _mm512_mul_ps(ctht2,ctht2);
                        alp2   = _mm512_mul_ps(x0,x0);
                        chtht1 = xcosf(htht1);
                        eai    = k0r;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        sect1  = _mm512_rcp14_ps(chtht1);
                        den1   = _mm512_add_ps(ctht1,ctht2);
                        stht1  = xsinf(tht1);
                        t0r    = _mm512_mul_ps(cer,inv);
                        chtht2 = xcosf(htht2);
                        t0i    = _mm512_mul_ps(cei,inv);
                        sect2  = _mm512_rcp14_ps(chtht2);
                        stht2  = xsinf(tht2);
                        num2   = _mm512_fmadd_ps(stht1,stht1,_mm512_mul_ps(stht2,stht2));
                        num1   = _mm512_fmadd_ps(sect1,sect1,_mm512_mul_ps(sect2,sect2));
                        x0     = _mm512_mul_ps(sphid,alp2);
                        rat1   = _mm512_div_ps(num1,den1);
                        x1     = _mm512_mul_ps(_2,_mm512_mul_ps(cx1,cx2));
                        sect1  = _mm512_add_ps(ctht1,ctht2);
                        sect2  = _mm512_mul_ps(sect1,sect1);
                        cer    = _mm512_mul_ps(t0r,x0);
                        den2   = _mm512_mul_ps(x1,sect2);
                        cei    = _mm512_mul_ps(t0i,x0);
                        rat2   = _mm512_div_ps(num2,den2);
                        *ESr   = _mm512_fmadd_ps(cer,rat1,rat2);
                        *ESi   = _mm512_fmadd_ps(cei,rat1,rat2);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f6214_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) palp,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht1,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht2,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64) ESi) {

                        register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        register __m512 r    = _mm512_load_ps(&pr[0]);
                        register __m512 alp  = _mm512_load_ps(&palp[0]);
                        register __m512 tht1 = _mm512_load_ps(&ptht1[0]);
                        register __m512 tht2 = _mm512_load_ps(&ptht2[0]);
                        register __m512 phi1 = _mm512_load_ps(&pphi1[0]);
                        register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,htht1,htht2,phid,ear,eai,cer,cei;
                        register __m512 sphid,t0r,t0i,num1,den1,num2,den2,inv;
                        register __m512 ctht1,ctht2,sect1,sect2,x0,x1,alp2;
                        register __m512 stht1,stht2,chtht1,chtht2,rat1,rat2,cx0,cx1;
                        phid   = _mm512_sub_ps(phi1,phi2);
                        k0r    = _mm512_mul_ps(k0,r);
                        htht1  = _mm512_mul_ps(tht1,hlf);
                        sphid  = xsinf(phid);
                        htht2  = _mm512_mul_ps(tht2,hlf);
                        inv    = _mm512_rcp14_ps(k0r);
                        ear    = _mm512_setzero_ps();
                        ctht1  = xcosf(tht1);
                        cx1    = _mm512_mul_ps(ctht1);
                        x0     = _mm512_mul_ps(alp,hlf);
                        ctht2  = xcosf(tht2);
                        cx2    = _mm512_mul_ps(ctht2,ctht2);
                        alp2   = _mm512_mul_ps(x0,x0);
                        chtht1 = xcosf(htht1);
                        eai    = k0r;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        sect1  = _mm512_rcp14_ps(chtht1);
                        den1   = _mm512_add_ps(ctht1,ctht2);
                        stht1  = xsinf(tht1);
                        t0r    = _mm512_mul_ps(cer,inv);
                        chtht2 = xcosf(htht2);
                        t0i    = _mm512_mul_ps(cei,inv);
                        sect2  = _mm512_rcp14_ps(chtht2);
                        stht2  = xsinf(tht2);
                        num2   = _mm512_fmadd_ps(stht1,stht1,_mm512_mul_ps(stht2,stht2));
                        num1   = _mm512_fmadd_ps(sect1,sect1,_mm512_mul_ps(sect2,sect2));
                        x0     = _mm512_mul_ps(sphid,alp2);
                        rat1   = _mm512_div_ps(num1,den1);
                        x1     = _mm512_mul_ps(_2,_mm512_mul_ps(cx1,cx2));
                        sect1  = _mm512_add_ps(ctht1,ctht2);
                        sect2  = _mm512_mul_ps(sect1,sect1);
                        cer    = _mm512_mul_ps(t0r,x0);
                        den2   = _mm512_mul_ps(x1,sect2);
                        cei    = _mm512_mul_ps(t0i,x0);
                        rat2   = _mm512_div_ps(num2,den2);
                        _mm512_store_ps(&ESr[0] ,_mm512_fmadd_ps(cer,rat1,rat2));
                        _mm512_store_ps(&ESi[0] ,_mm512_fmadd_ps(cei,rat1,rat2));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f6214_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict  pr,
                                           const float * __restrict  palp,
                                           const float * __restrict  ptht1,
                                           const float * __restrict  ptht2,
                                           const float * __restrict  pphi1,
                                           const float * __restrict  pphi2,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {

                        register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        register __m512 alp  = _mm512_loadu_ps(&palp[0]);
                        register __m512 tht1 = _mm512_loadu_ps(&ptht1[0]);
                        register __m512 tht2 = _mm512_loadu_ps(&ptht2[0]);
                        register __m512 phi1 = _mm512_loadu_ps(&pphi1[0]);
                        register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,htht1,htht2,phid,ear,eai,cer,cei;
                        register __m512 sphid,t0r,t0i,num1,den1,num2,den2,inv;
                        register __m512 ctht1,ctht2,sect1,sect2,x0,x1,alp2;
                        register __m512 stht1,stht2,chtht1,chtht2,rat1,rat2,cx0,cx1;
                        phid   = _mm512_sub_ps(phi1,phi2);
                        k0r    = _mm512_mul_ps(k0,r);
                        htht1  = _mm512_mul_ps(tht1,hlf);
                        sphid  = xsinf(phid);
                        htht2  = _mm512_mul_ps(tht2,hlf);
                        inv    = _mm512_rcp14_ps(k0r);
                        ear    = _mm512_setzero_ps();
                        ctht1  = xcosf(tht1);
                        cx1    = _mm512_mul_ps(ctht1);
                        x0     = _mm512_mul_ps(alp,hlf);
                        ctht2  = xcosf(tht2);
                        cx2    = _mm512_mul_ps(ctht2,ctht2);
                        alp2   = _mm512_mul_ps(x0,x0);
                        chtht1 = xcosf(htht1);
                        eai    = k0r;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        sect1  = _mm512_rcp14_ps(chtht1);
                        den1   = _mm512_add_ps(ctht1,ctht2);
                        stht1  = xsinf(tht1);
                        t0r    = _mm512_mul_ps(cer,inv);
                        chtht2 = xcosf(htht2);
                        t0i    = _mm512_mul_ps(cei,inv);
                        sect2  = _mm512_rcp14_ps(chtht2);
                        stht2  = xsinf(tht2);
                        num2   = _mm512_fmadd_ps(stht1,stht1,_mm512_mul_ps(stht2,stht2));
                        num1   = _mm512_fmadd_ps(sect1,sect1,_mm512_mul_ps(sect2,sect2));
                        x0     = _mm512_mul_ps(sphid,alp2);
                        rat1   = _mm512_div_ps(num1,den1);
                        x1     = _mm512_mul_ps(_2,_mm512_mul_ps(cx1,cx2));
                        sect1  = _mm512_add_ps(ctht1,ctht2);
                        sect2  = _mm512_mul_ps(sect1,sect1);
                        cer    = _mm512_mul_ps(t0r,x0);
                        den2   = _mm512_mul_ps(x1,sect2);
                        cei    = _mm512_mul_ps(t0i,x0);
                        rat2   = _mm512_div_ps(num2,den2);
                        _mm512_storeu_ps(&ESr[0] ,_mm512_fmadd_ps(cer,rat1,rat2));
                        _mm512_storeu_ps(&ESi[0] ,_mm512_fmadd_ps(cei,rat1,rat2));
                }


                    /*
                           Bistatic RCS case.
                           E-field scattered for (phi component).
                           Formula 6.2-15
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6215_zmm16r4(const __m512 k0,
                                           const __m512 r,
                                           const __m512 alp,
                                           const __m512 tht1, //inc
                                           const __m512 tht2  //scat
                                           const __m512 phi1, //inc
                                           const __m512 phi2, //scat
                                           __m512 * __restrict ESr,
                                           __m512 * __restrict ESi) {

                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 _16 = _mm512_set1_ps(16.0f);
                        register __m512 k0r,inv,alph,ear,eai,cer,cei;
                        register __m512 t0r,t0i,trm1,trm2,trm3,trm4;
                        register __m512 stht1,stht2,ctht1,ctht2;
                        register __m512 htht1,htht2,chtht1,chtht2;
                        register __m512 shtht1,shtht2,x0,x1,phid,cphid;
                        register __m512 ctht12,sect1,sect2,x2,x3,cthtp3;
                        k0r   = _mm512_mul_ps(k0,r);
                        stht1 = xsinf(tht1);
                        alph  = _mm512_mul_ps(alp,hlf);
                        ctht1 = xcosf(tht1);
                        phid  = _mm512_sub_ps(phi1,phi2);
                        stht2 = xsinf(tht2);
                        htht1 = _mm512_mul_ps(tht1,hlf);
                        ctht2 = xcosf(tht2);
                        ear   = _mm512_setzero_ps();
                        htht2 = _mm512_mul_ps(tht2);
                        eai   = k0r;
                        inv   = _mm512_rcp14_ps(k0r);
                        ctht12= _mm512_add_ps(ctht1,ctht2);
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        x0    = _mm512_fmadd_ps(stht2,stht1,stht1);
                        chtht1= xcosf(htht1);
                        cer   = _mm512_mul_ps(_mm512_mul_ps(t0r,inv),alph);
                        cthtp3= _mm512_mul_ps(ctht12,_mm512_mul_ps(ctht12,ctht12));
                        cei   = _mm512_mul_ps(_mm512_mul_ps(t0i,inv),alph);
                        chtht2= xcosf(htht2);
                        trm1  = _mm512_div_ps(x0,cthtp3);
                        x2    = _mm512_rcp14_ps(chtht1);
                        shtht1= xsinf(htht1);
                        x3    = _mm512_rcp14_ps(chtht2);
                        shtht2= xsinf(htht2);
                        x0    = _mm512_fmadd_ps(x2,x2,_mm512_mul_ps(x3,x3));
                        trm2  = _mm512_div_ps(x0,ctht12);
                        x1    = _mm512_fmadd_ps(stht1,stht1,_mm512_mul_ps(stht2,stht2));
                        x2    = _mm512_mul_ps(_2,_mm512_mul_ps(_mm512_mul_ps(chtht1,chtht1),
                                              _mm512_mul_ps(chtht2,chtht2)));
                        x3    = _mm512_mul_ps(ctht12,ctht12);
                        t0r   = _mm512_mul_ps(x2,x3);
                        trm3  = _mm512_div_ps(x1,t0r);
                        x0    = _mm512_mul_ps(_16,_mm512_mul_ps(shtht1,shtht1));
                        cphid = xcosf(phid);
                        x1    = _mm512_mul_ps(shtht2,shtht2);
                        x2    = _mm512_mul_ps(x0,x1);
                        trm4  = _mm512_div_ps(x2,cthtp3);
                        x3    = _mm512_add_ps(trm1,cphid);
                        x4    = _mm512_add_ps(trm2,_mm512_add_ps(trm3,trm4));
                        x0    = _mm512_mul_ps(x3,x4);
                        *ESr  = _mm512_mul_ps(cer,x0);
                        *ESi  = _mm512_mul_ps(cei,x0);
                                                
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6215_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) palp,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht1,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht2,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64) ESi) {

                        register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        register __m512 r    = _mm512_load_ps(&pr[0]);
                        register __m512 alp  = _mm512_load_ps(&palp[0]);
                        register __m512 tht1 = _mm512_load_ps(&ptht1[0]);
                        register __m512 tht2 = _mm512_load_ps(&ptht2[0]);
                        register __m512 phi1 = _mm512_load_ps(&pphi1[0]);
                        register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 _16 = _mm512_set1_ps(16.0f);
                        register __m512 k0r,inv,alph,ear,eai,cer,cei;
                        register __m512 t0r,t0i,trm1,trm2,trm3,trm4;
                        register __m512 stht1,stht2,ctht1,ctht2;
                        register __m512 htht1,htht2,chtht1,chtht2;
                        register __m512 shtht1,shtht2,x0,x1,phid,cphid;
                        register __m512 ctht12,sect1,sect2,x2,x3,cthtp3;
                        k0r   = _mm512_mul_ps(k0,r);
                        stht1 = xsinf(tht1);
                        alph  = _mm512_mul_ps(alp,hlf);
                        ctht1 = xcosf(tht1);
                        phid  = _mm512_sub_ps(phi1,phi2);
                        stht2 = xsinf(tht2);
                        htht1 = _mm512_mul_ps(tht1,hlf);
                        ctht2 = xcosf(tht2);
                        ear   = _mm512_setzero_ps();
                        htht2 = _mm512_mul_ps(tht2);
                        eai   = k0r;
                        inv   = _mm512_rcp14_ps(k0r);
                        ctht12= _mm512_add_ps(ctht1,ctht2);
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        x0    = _mm512_fmadd_ps(stht2,stht1,stht1);
                        chtht1= xcosf(htht1);
                        cer   = _mm512_mul_ps(_mm512_mul_ps(t0r,inv),alph);
                        cthtp3= _mm512_mul_ps(ctht12,_mm512_mul_ps(ctht12,ctht12));
                        cei   = _mm512_mul_ps(_mm512_mul_ps(t0i,inv),alph);
                        chtht2= xcosf(htht2);
                        trm1  = _mm512_div_ps(x0,cthtp3);
                        x2    = _mm512_rcp14_ps(chtht1);
                        shtht1= xsinf(htht1);
                        x3    = _mm512_rcp14_ps(chtht2);
                        shtht2= xsinf(htht2);
                        x0    = _mm512_fmadd_ps(x2,x2,_mm512_mul_ps(x3,x3));
                        trm2  = _mm512_div_ps(x0,ctht12);
                        x1    = _mm512_fmadd_ps(stht1,stht1,_mm512_mul_ps(stht2,stht2));
                        x2    = _mm512_mul_ps(_2,_mm512_mul_ps(_mm512_mul_ps(chtht1,chtht1),
                                              _mm512_mul_ps(chtht2,chtht2)));
                        x3    = _mm512_mul_ps(ctht12,ctht12);
                        t0r   = _mm512_mul_ps(x2,x3);
                        trm3  = _mm512_div_ps(x1,t0r);
                        x0    = _mm512_mul_ps(_16,_mm512_mul_ps(shtht1,shtht1));
                        cphid = xcosf(phid);
                        x1    = _mm512_mul_ps(shtht2,shtht2);
                        x2    = _mm512_mul_ps(x0,x1);
                        trm4  = _mm512_div_ps(x2,cthtp3);
                        x3    = _mm512_add_ps(trm1,cphid);
                        x4    = _mm512_add_ps(trm2,_mm512_add_ps(trm3,trm4));
                        x0    = _mm512_mul_ps(x3,x4);
                        _mm512_store_ps(&ESr[0] ,_mm512_mul_ps(cer,x0));
                        _mm512_store_ps(&ESi[0] ,_mm512_mul_ps(cei,x0));
                                                
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6215_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict  pr,
                                           const float * __restrict  palp,
                                           const float * __restrict  ptht1,
                                           const float * __restrict  ptht2,
                                           const float * __restrict  pphi1,
                                           const float * __restrict  pphi2,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {

                        register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        register __m512 alp  = _mm512_loadu_ps(&palp[0]);
                        register __m512 tht1 = _mm512_loadu_ps(&ptht1[0]);
                        register __m512 tht2 = _mm512_loadu_ps(&ptht2[0]);
                        register __m512 phi1 = _mm512_loadu_ps(&pphi1[0]);
                        register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 _16 = _mm512_set1_ps(16.0f);
                        register __m512 k0r,inv,alph,ear,eai,cer,cei;
                        register __m512 t0r,t0i,trm1,trm2,trm3,trm4;
                        register __m512 stht1,stht2,ctht1,ctht2;
                        register __m512 htht1,htht2,chtht1,chtht2;
                        register __m512 shtht1,shtht2,x0,x1,phid,cphid;
                        register __m512 ctht12,sect1,sect2,x2,x3,cthtp3;
                        k0r   = _mm512_mul_ps(k0,r);
                        stht1 = xsinf(tht1);
                        alph  = _mm512_mul_ps(alp,hlf);
                        ctht1 = xcosf(tht1);
                        phid  = _mm512_sub_ps(phi1,phi2);
                        stht2 = xsinf(tht2);
                        htht1 = _mm512_mul_ps(tht1,hlf);
                        ctht2 = xcosf(tht2);
                        ear   = _mm512_setzero_ps();
                        htht2 = _mm512_mul_ps(tht2);
                        eai   = k0r;
                        inv   = _mm512_rcp14_ps(k0r);
                        ctht12= _mm512_add_ps(ctht1,ctht2);
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        x0    = _mm512_fmadd_ps(stht2,stht1,stht1);
                        chtht1= xcosf(htht1);
                        cer   = _mm512_mul_ps(_mm512_mul_ps(t0r,inv),alph);
                        cthtp3= _mm512_mul_ps(ctht12,_mm512_mul_ps(ctht12,ctht12));
                        cei   = _mm512_mul_ps(_mm512_mul_ps(t0i,inv),alph);
                        chtht2= xcosf(htht2);
                        trm1  = _mm512_div_ps(x0,cthtp3);
                        x2    = _mm512_rcp14_ps(chtht1);
                        shtht1= xsinf(htht1);
                        x3    = _mm512_rcp14_ps(chtht2);
                        shtht2= xsinf(htht2);
                        x0    = _mm512_fmadd_ps(x2,x2,_mm512_mul_ps(x3,x3));
                        trm2  = _mm512_div_ps(x0,ctht12);
                        x1    = _mm512_fmadd_ps(stht1,stht1,_mm512_mul_ps(stht2,stht2));
                        x2    = _mm512_mul_ps(_2,_mm512_mul_ps(_mm512_mul_ps(chtht1,chtht1),
                                              _mm512_mul_ps(chtht2,chtht2)));
                        x3    = _mm512_mul_ps(ctht12,ctht12);
                        t0r   = _mm512_mul_ps(x2,x3);
                        trm3  = _mm512_div_ps(x1,t0r);
                        x0    = _mm512_mul_ps(_16,_mm512_mul_ps(shtht1,shtht1));
                        cphid = xcosf(phid);
                        x1    = _mm512_mul_ps(shtht2,shtht2);
                        x2    = _mm512_mul_ps(x0,x1);
                        trm4  = _mm512_div_ps(x2,cthtp3);
                        x3    = _mm512_add_ps(trm1,cphid);
                        x4    = _mm512_add_ps(trm2,_mm512_add_ps(trm3,trm4));
                        x0    = _mm512_mul_ps(x3,x4);
                        _mm512_storeu_ps(&ESr[0] ,_mm512_mul_ps(cer,x0));
                        _mm512_storeu_ps(&ESi[0] ,_mm512_mul_ps(cei,x0));
                                                
                }


                  /*
                           Bistatic RCS case -- Physical Optics approximated.
                           E-field scattered for (theta component).
                           Formula 6.2-17
                    */

#include "GMS_simd_utils.hpp"

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f6217_zmm16r4(const __m512 k0,
                                           const __m512 r,
                                           const __m512 alp,
                                           const __m512 tht,
                                           const __m512 phi,
                                           __m512 * __restrict ESr,
                                           __m512 * __restrict ESi) {

                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 c0  = _mm512_set1_ps(1.5f);
                         register __m512 k0r,ear,eai,cer,cei,t0r,t0i,inv,x0,chtht;
                         register __m512 L,cphi,cosa,sina,htht,argm,argp,num,den;
                         k0r   = _mm512_mul_ps(k0,r);
                         cphi  = xcosf(phi);
                         htht  = _mm512_mul_ps(tht,hlf);
                         cosa  = xcosf(alp);
                         ear   = _mm512_setzero_ps();
                         inv   = _mm512_rcp14_ps(k0r);
                         eai   = k0r;
                         argm  = _mm512_add_ps(alp,htht);
                         cexp_zmm16r4(ear,eai,&t0r,&t0i);
                         argp  = _mm512_sub_ps(alp,htht);
                         sina  = xsinf(alp);
                         x0    = _mm512_mul_ps(sina,sina);
                         chtht = xcosf(htht);
                         cer   = _mm512_mul_ps(t0r,inv);
                         num   = _mm512_mul_ps(negate_zmm16r4(x0),cosa);
                         cei   = _mm512_mul_ps(t0i,inv);
                         ear   = xcosf(argp);
                         eai   = xcosf(argm);
                         x0    = _mm512_pow_ps(_mm512_add_ps(ear,eai),c0);
                         den   = _mm512_mul_ps(_mm512_mul_ps(_4,chtht),x0);
                         L     = _mm512_div_ps(num,den);
                         t0i   = _mm512_mul_ps(L,cphi);
                         *ESr  = _mm512_mul_ps(cer,t0i);
                         *ESi  = _mm512_mul_ps(cei,t0i);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f6217_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) palp,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64) ESi) {

                         register __m512 k0   = _mm512_load_ps(&pk0[0]);
                         register __m512 r    = _mm512_load_ps(&pr[0]);
                         register __m512 alp  = _mm512_load_ps(&palp[0]);
                         register __m512 tht  = _mm512_load_ps(&ptht[0]);
                         register __m512 phi = _mm512_load_ps(&pphi[0]);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 c0  = _mm512_set1_ps(1.5f);
                         register __m512 k0r,ear,eai,cer,cei,t0r,t0i,inv,x0,chtht;
                         register __m512 L,cphi,cosa,sina,htht,argm,argp,num,den;
                         k0r   = _mm512_mul_ps(k0,r);
                         cphi  = xcosf(phi);
                         htht  = _mm512_mul_ps(tht,hlf);
                         cosa  = xcosf(alp);
                         ear   = _mm512_setzero_ps();
                         inv   = _mm512_rcp14_ps(k0r);
                         eai   = k0r;
                         argm  = _mm512_add_ps(alp,htht);
                         cexp_zmm16r4(ear,eai,&t0r,&t0i);
                         argp  = _mm512_sub_ps(alp,htht);
                         sina  = xsinf(alp);
                         x0    = _mm512_mul_ps(sina,sina);
                         chtht = xcosf(htht);
                         cer   = _mm512_mul_ps(t0r,inv);
                         num   = _mm512_mul_ps(negate_zmm16r4(x0),cosa);
                         cei   = _mm512_mul_ps(t0i,inv);
                         ear   = xcosf(argp);
                         eai   = xcosf(argm);
                         x0    = _mm512_pow_ps(_mm512_add_ps(ear,eai),c0);
                         den   = _mm512_mul_ps(_mm512_mul_ps(_4,chtht),x0);
                         L     = _mm512_div_ps(num,den);
                         t0i   = _mm512_mul_ps(L,cphi);
                         _mm512_store_ps(&ESr[0] ,_mm512_mul_ps(cer,t0i));
                         _mm512_store_ps(&ESi[0] ,_mm512_mul_ps(cei,t0i));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESth_f6217_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict  pr,
                                           const float * __restrict  palp,
                                           const float * __restrict  ptht,
                                           const float * __restrict  pphi,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {

                         register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                         register __m512 r    = _mm512_loadu_ps(&pr[0]);
                         register __m512 alp  = _mm512_loadu_ps(&palp[0]);
                         register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                         register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 c0  = _mm512_set1_ps(1.5f);
                         register __m512 k0r,ear,eai,cer,cei,t0r,t0i,inv,x0,chtht;
                         register __m512 L,cphi,cosa,sina,htht,argm,argp,num,den;
                         k0r   = _mm512_mul_ps(k0,r);
                         cphi  = xcosf(phi);
                         htht  = _mm512_mul_ps(tht,hlf);
                         cosa  = xcosf(alp);
                         ear   = _mm512_setzero_ps();
                         inv   = _mm512_rcp14_ps(k0r);
                         eai   = k0r;
                         argm  = _mm512_add_ps(alp,htht);
                         cexp_zmm16r4(ear,eai,&t0r,&t0i);
                         argp  = _mm512_sub_ps(alp,htht);
                         sina  = xsinf(alp);
                         x0    = _mm512_mul_ps(sina,sina);
                         chtht = xcosf(htht);
                         cer   = _mm512_mul_ps(t0r,inv);
                         num   = _mm512_mul_ps(negate_zmm16r4(x0),cosa);
                         cei   = _mm512_mul_ps(t0i,inv);
                         ear   = xcosf(argp);
                         eai   = xcosf(argm);
                         x0    = _mm512_pow_ps(_mm512_add_ps(ear,eai),c0);
                         den   = _mm512_mul_ps(_mm512_mul_ps(_4,chtht),x0);
                         L     = _mm512_div_ps(num,den);
                         t0i   = _mm512_mul_ps(L,cphi);
                         _mm512_storeu_ps(&ESr[0] ,_mm512_mul_ps(cer,t0i));
                         _mm512_storeu_ps(&ESi[0] ,_mm512_mul_ps(cei,t0i));
                }


                    /*
                           Bistatic RCS case -- Physical Optics approximated.
                           E-field scattered for (phi component).
                           Formula 6.2-18
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6218_zmm16r4(const __m512 k0,
                                           const __m512 r,
                                           const __m512 alp,
                                           const __m512 tht,
                                           const __m512 phi,
                                           __m512 * __restrict ESr,
                                           __m512 * __restrict ESi) {
                        
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 c0  = _mm512_set1_ps(1.5f);
                         const __m512 n1  = _mm512_set1_ps(-1.0f);
                         register __m512 k0r,ear,eai,cer,cei,t0r,t0i,inv,x0,chtht;
                         register __m512 L,sphi,cosa,sina,htht,argm,argp,num,den;
                         k0r   = _mm512_mul_ps(k0,r);
                         sphi  = xsinf(phi);
                         htht  = _mm512_mul_ps(tht,hlf);
                         cosa  = xcosf(alp);
                         ear   = _mm512_setzero_ps();
                         inv   = _mm512_rcp14_ps(k0r);
                         eai   = k0r;
                         argm  = _mm512_add_ps(alp,htht);
                         cexp_zmm16r4(ear,eai,&t0r,&t0i);
                         argp  = _mm512_sub_ps(alp,htht);
                         sina  = xsinf(alp);
                         x0    = _mm512_mul_ps(sina,sina);
                         chtht = xcosf(htht);
                         cer   = _mm512_mul_ps(t0r,inv);
                         num   = _mm512_mul_ps(negate_zmm16r4(x0),cosa);
                         cei   = _mm512_mul_ps(t0i,inv);
                         ear   = xcosf(argp);
                         eai   = xcosf(argm);
                         x0    = _mm512_pow_ps(_mm512_add_ps(ear,eai),c0);
                         den   = _mm512_mul_ps(_mm512_mul_ps(_4,chtht),x0);
                         L     = _mm512_div_ps(num,den);
                         t0i   = _mm512_mul_ps(L,cphi);
                         *ESr  = _mm512_mul_ps(n1,_mm512_mul_ps(cer,t0i)));
                         *ESi  = _mm512_mul_ps(n1,_mm512_mul_ps(cei,t0i)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6218_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) palp,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64) ESi) {

                         register __m512 k0   = _mm512_load_ps(&pk0[0]);
                         register __m512 r    = _mm512_load_ps(&pr[0]);
                         register __m512 alp  = _mm512_load_ps(&palp[0]);
                         register __m512 tht  = _mm512_load_ps(&ptht[0]);
                         register __m512 phi = _mm512_load_ps(&pphi[0]);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 c0  = _mm512_set1_ps(1.5f);
                         const __m512 n1  = _mm512_set1_ps(-1.0f);
                         register __m512 k0r,ear,eai,cer,cei,t0r,t0i,inv,x0,chtht;
                         register __m512 L,sphi,cosa,sina,htht,argm,argp,num,den;
                         k0r   = _mm512_mul_ps(k0,r);
                         sphi  = xsinf(phi);
                         htht  = _mm512_mul_ps(tht,hlf);
                         cosa  = xcosf(alp);
                         ear   = _mm512_setzero_ps();
                         inv   = _mm512_rcp14_ps(k0r);
                         eai   = k0r;
                         argm  = _mm512_add_ps(alp,htht);
                         cexp_zmm16r4(ear,eai,&t0r,&t0i);
                         argp  = _mm512_sub_ps(alp,htht);
                         sina  = xsinf(alp);
                         x0    = _mm512_mul_ps(sina,sina);
                         chtht = xcosf(htht);
                         cer   = _mm512_mul_ps(t0r,inv);
                         num   = _mm512_mul_ps(negate_zmm16r4(x0),cosa);
                         cei   = _mm512_mul_ps(t0i,inv);
                         ear   = xcosf(argp);
                         eai   = xcosf(argm);
                         x0    = _mm512_pow_ps(_mm512_add_ps(ear,eai),c0);
                         den   = _mm512_mul_ps(_mm512_mul_ps(_4,chtht),x0);
                         L     = _mm512_div_ps(num,den);
                         t0i   = _mm512_mul_ps(L,cphi);
                         _mm512_store_ps(&ESr[0]  ,_mm512_mul_ps(n1,_mm512_mul_ps(cer,t0i)));
                         _mm512_store_ps(&ESi[0]  ,_mm512_mul_ps(n1,_mm512_mul_ps(cei,t0i)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ESph_f6218_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict  pr,
                                           const float * __restrict  palp,
                                           const float * __restrict  ptht,
                                           const float * __restrict  pphi,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {

                         register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                         register __m512 r    = _mm512_loadu_ps(&pr[0]);
                         register __m512 alp  = _mm512_loadu_ps(&palp[0]);
                         register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                         register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 c0  = _mm512_set1_ps(1.5f);
                         const __m512 n1  = _mm512_set1_ps(-1.0f);
                         register __m512 k0r,ear,eai,cer,cei,t0r,t0i,inv,x0,chtht;
                         register __m512 L,sphi,cosa,sina,htht,argm,argp,num,den;
                         k0r   = _mm512_mul_ps(k0,r);
                         sphi  = xsinf(phi);
                         htht  = _mm512_mul_ps(tht,hlf);
                         cosa  = xcosf(alp);
                         ear   = _mm512_setzero_ps();
                         inv   = _mm512_rcp14_ps(k0r);
                         eai   = k0r;
                         argm  = _mm512_add_ps(alp,htht);
                         cexp_zmm16r4(ear,eai,&t0r,&t0i);
                         argp  = _mm512_sub_ps(alp,htht);
                         sina  = xsinf(alp);
                         x0    = _mm512_mul_ps(sina,sina);
                         chtht = xcosf(htht);
                         cer   = _mm512_mul_ps(t0r,inv);
                         num   = _mm512_mul_ps(negate_zmm16r4(x0),cosa);
                         cei   = _mm512_mul_ps(t0i,inv);
                         ear   = xcosf(argp);
                         eai   = xcosf(argm);
                         x0    = _mm512_pow_ps(_mm512_add_ps(ear,eai),c0);
                         den   = _mm512_mul_ps(_mm512_mul_ps(_4,chtht),x0);
                         L     = _mm512_div_ps(num,den);
                         t0i   = _mm512_mul_ps(L,cphi);
                         _mm512_storeu_ps(&ESr[0]  ,_mm512_mul_ps(n1,_mm512_mul_ps(cer,t0i)));
                         _mm512_storeu_ps(&ESi[0]  ,_mm512_mul_ps(n1,_mm512_mul_ps(cei,t0i)));
                }



                    /*
                           Physical-Optics axial-incidence bistatic RCS
                           function of theta for (0 << theta < PI-2*alpha)
                           Formula 6.2-20
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6220_zmm16r4(const __m512 gam0,
                                            const __m512 alp,
                                            const __m512 tht) {

                           const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                           const __m512 _1    = _mm512_set1_ps(1.0f);
                           register __m512 rcs,gam2,tana,tan4,x0,x1,trm1;
                           register __m512 alp2,calp2,ctht,trm2,num,den,x2;
                           gam2  = _mm512_mul_ps(gam0,gam0);
                           tana  = xtanf(alp);
                           alp2  = _mm512_add_ps(alp,alp);
                           calp2 = xcosf(alp2);
                           x0    = _mm512_mul_ps(tana,tana);
                           x1    = _mm512_add_ps(_1,calp2);
                           x2    = _mm512_mul_ps(x0,x0);
                           ctht  = xcosf(tht);
                           trm1  = _mm512_div_ps(_mm512_mul_ps(gam2,x2),_16pi);
                           x0    = _mm512_add_ps(_1,ctht);
                           x2    = _mm512_mul_ps(x1,_mm512_mul_ps(x1,x1));
                           num   = _mm512_add_ps(x2,x2);
                           x1    = _mm512_add_ps(ctht,calp2);
                           x2    = _mm512_mul_ps(x1,_mm512_mul_ps(x1,x1));
                           den   = _mm512_mul_ps(x0,x2);
                           trm2  = _mm512_div_ps(num,den);
                           rcs   = _mm512_mul_ps(trm1,trm2);
                           return (rcs);
                  }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6220_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam0,
                                            const float * __restrict __ATTR_ALIGN__(64) palp,
                                            const float * __restrict __ATTR_ALIGN__(64) ptht) {

                           register __m512  gam0  = _mm512_load_ps(&pgam0[0]);
                           register __m512  alp   = _mm512_load_ps(&palp[0]);
                           register __m512  tht   = _mm512_load_ps(&ptht[0]);
                           const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                           const __m512 _1    = _mm512_set1_ps(1.0f);
                           register __m512 rcs,gam2,tana,tan4,x0,x1,trm1;
                           register __m512 alp2,calp2,ctht,trm2,num,den,x2;
                           gam2  = _mm512_mul_ps(gam0,gam0);
                           tana  = xtanf(alp);
                           alp2  = _mm512_add_ps(alp,alp);
                           calp2 = xcosf(alp2);
                           x0    = _mm512_mul_ps(tana,tana);
                           x1    = _mm512_add_ps(_1,calp2);
                           x2    = _mm512_mul_ps(x0,x0);
                           ctht  = xcosf(tht);
                           trm1  = _mm512_div_ps(_mm512_mul_ps(gam2,x2),_16pi);
                           x0    = _mm512_add_ps(_1,ctht);
                           x2    = _mm512_mul_ps(x1,_mm512_mul_ps(x1,x1));
                           num   = _mm512_add_ps(x2,x2);
                           x1    = _mm512_add_ps(ctht,calp2);
                           x2    = _mm512_mul_ps(x1,_mm512_mul_ps(x1,x1));
                           den   = _mm512_mul_ps(x0,x2);
                           trm2  = _mm512_div_ps(num,den);
                           rcs   = _mm512_mul_ps(trm1,trm2);
                           return (rcs);
                  }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6220_zmm16r4_u(const float * __restrict  pgam0,
                                            const float * __restrict  palp,
                                            const float * __restrict  ptht) {

                           register __m512  gam0  = _mm512_loadu_ps(&pgam0[0]);
                           register __m512  alp   = _mm512_loadu_ps(&palp[0]);
                           register __m512  tht   = _mm512_loadu_ps(&ptht[0]);
                           const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                           const __m512 _1    = _mm512_set1_ps(1.0f);
                           register __m512 rcs,gam2,tana,tan4,x0,x1,trm1;
                           register __m512 alp2,calp2,ctht,trm2,num,den,x2;
                           gam2  = _mm512_mul_ps(gam0,gam0);
                           tana  = xtanf(alp);
                           alp2  = _mm512_add_ps(alp,alp);
                           calp2 = xcosf(alp2);
                           x0    = _mm512_mul_ps(tana,tana);
                           x1    = _mm512_add_ps(_1,calp2);
                           x2    = _mm512_mul_ps(x0,x0);
                           ctht  = xcosf(tht);
                           trm1  = _mm512_div_ps(_mm512_mul_ps(gam2,x2),_16pi);
                           x0    = _mm512_add_ps(_1,ctht);
                           x2    = _mm512_mul_ps(x1,_mm512_mul_ps(x1,x1));
                           num   = _mm512_add_ps(x2,x2);
                           x1    = _mm512_add_ps(ctht,calp2);
                           x2    = _mm512_mul_ps(x1,_mm512_mul_ps(x1,x1));
                           den   = _mm512_mul_ps(x0,x2);
                           trm2  = _mm512_div_ps(num,den);
                           rcs   = _mm512_mul_ps(trm1,trm2);
                           return (rcs);
                  }


                     /*
                           Axial incidence backscattering (theta = 0)
                           RCS (Spencer equation).
                           Formula 6.2-22
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6222_zmm16r4(const __m512 gam0,
                                        const __m512 alp) {

                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          register __m512 rcs,gam2,trm1,tana,x0,x1;
                          gam2 = _mm512_mul_ps(gam0,gam00;
                          tana = xtanf(alp);
                          trm1 = _mm512_div_ps(gam2,_16pi);
                          x0   = _mm512_mul_ps(tana,tana);
                          x1   = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,x1);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6222_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam0,
                                            const float * __restrict __ATTR_ALIGN__(64) palp) {

                          register __m512  gam0  = _mm512_load_ps(&pgam0[0]);
                          register __m512  alp   = _mm512_load_ps(&palp[0]);
                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          register __m512 rcs,gam2,trm1,tana,x0,x1;
                          gam2 = _mm512_mul_ps(gam0,gam00;
                          tana = xtanf(alp);
                          trm1 = _mm512_div_ps(gam2,_16pi);
                          x0   = _mm512_mul_ps(tana,tana);
                          x1   = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,x1);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6222_zmm16r4_u(const float * __restrict  pgam0,
                                            const float * __restrict  palp) {

                          register __m512  gam0  = _mm512_loadu_ps(&pgam0[0]);
                          register __m512  alp   = _mm512_loadu_ps(&palp[0]);
                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          register __m512 rcs,gam2,trm1,tana,x0,x1;
                          gam2 = _mm512_mul_ps(gam0,gam00;
                          tana = xtanf(alp);
                          trm1 = _mm512_div_ps(gam2,_16pi);
                          x0   = _mm512_mul_ps(tana,tana);
                          x1   = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,x1);
                          return (rcs);
                 }


                    /*
                           Narrow-angle cone.
                           Scattered E-field.
                           Formula 6.2-24
                       */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f6224_zmm16r4(const __m512 k0,
                                         const __m512 z,
                                         const __m512 alp,
                                         const __m512 x,
                                         const __m512 y,
                                         const __m512 z,
                                         __m512 * __restrict xre,
                                         __m512 * __restrict xim,
                                         __m512 * __restrict yre,
                                         __m512 * __restrict yim,
                                         __m512 * __restrict zre,
                                         __m512 * __restrict zim) {

                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 cr  = _mm512_set1_ps(0.0f);
                        const __m512 ci  = _mm512_set1_ps(1.0f);
                        register __m512 alph,x0,k0z,inv,ear,eai,cer,cei;
                        register __m512 t0r,t0i;
                        
                        alph = _mm512_mul_ps(alp,hlf);
                        k0z  = _mm512_mul_ps(k0,z);
                        ear  = _mm512_setzero_ps();
                        inv  = _mm512_rcp14_ps(k0z);
                        eai  = k0z;
                        x0   = _mm512_mul_ps(alph,alph);
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        cmul_zmm16r4(t0r,t0i,cr,ci,&cer,&cei);
                        cer  = _mm512_mul_ps(x0,_mm512_mul_ps(cer,inv));
                        cei  = _mm512_mul_ps(x0,_mm512_mul_ps(cei,inv));
                        *xre  = _mm512_mul_ps(x,cer);
                        *xim  = _mm512_mul_ps(x,cei);
                        *yre  = _mm512_mul_ps(y,cer);
                        *yim  = _mm512_mul_ps(y,cei);
                        *zre  = _mm512_mul_ps(z,cer);
                        *zim  = _mm512_mul_ps(z,cei);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f6224_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pz,
                                           const float * __restrict __ATTR_ALIGN__(64) palp,
                                           const float * __restrict __ATTR_ALIGN__(64) px,
                                           const float * __restrict __ATTR_ALIGN__(64) py,
                                           const float * __restrict __ATTR_ALIGN__(64) pz,
                                           float * __restrict __ATTR_ALIGN__(64) xre,
                                           float * __restrict __ATTR_ALIGN__(64) xim,
                                           float * __restrict __ATTR_ALIGN__(64) yre,
                                           float * __restrict __ATTR_ALIGN__(64) yim,
                                           float * __restrict __ATTR_ALIGN__(64) zre,
                                           float * __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 z   = _mm512_load_ps(&pz[0]);
                        register __m512 alp = _mm512_load_ps(&palp[0]);
                        register __m512 x   = _mm512_load_ps(&px[0]);
                        register __m512 y   = _mm512_load_ps(&py[0]);
                        register __m512 z   = _mm512_load_ps(&pz[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 cr  = _mm512_set1_ps(0.0f);
                        const __m512 ci  = _mm512_set1_ps(1.0f);
                        register __m512 alph,x0,k0z,inv,ear,eai,cer,cei;
                        register __m512 t0r,t0i;
                        
                        alph = _mm512_mul_ps(alp,hlf);
                        k0z  = _mm512_mul_ps(k0,z);
                        ear  = _mm512_setzero_ps();
                        inv  = _mm512_rcp14_ps(k0z);
                        eai  = k0z;
                        x0   = _mm512_mul_ps(alph,alph);
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        cmul_zmm16r4(t0r,t0i,cr,ci,&cer,&cei);
                        cer  = _mm512_mul_ps(x0,_mm512_mul_ps(cer,inv));
                        cei  = _mm512_mul_ps(x0,_mm512_mul_ps(cei,inv));
                        _mm512_store_ps(&xre[0] ,_mm512_mul_ps(x,cer));
                        _mm512_store_ps(&xim[0] ,_mm512_mul_ps(x,cei));
                        _mm512_store_ps(&yre[0] ,_mm512_mul_ps(y,cer));
                        _mm512_store_ps(&yim[0] ,_mm512_mul_ps(y,cei));
                        _mm512_store_ps(&zre[0] ,_mm512_mul_ps(z,cer));
                        _mm512_store_ps(&zim[0] ,_mm512_mul_ps(z,cei));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f6224_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict pz,
                                           const float * __restrict palp,
                                           const float * __restrict  px,
                                           const float * __restrict  py,
                                           const float * __restrict  pz,
                                           float * __restrict  xre,
                                           float * __restrict  xim,
                                           float * __restrict  yre,
                                           float * __restrict  yim,
                                           float * __restrict  zre,
                                           float * __restrict  zim) {

                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 z   = _mm512_loadu_ps(&pz[0]);
                        register __m512 alp = _mm512_loadu_ps(&palp[0]);
                        register __m512 x   = _mm512_loadu_ps(&px[0]);
                        register __m512 y   = _mm512_loadu_ps(&py[0]);
                        register __m512 z   = _mm512_loadu_ps(&pz[0]);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 cr  = _mm512_set1_ps(0.0f);
                        const __m512 ci  = _mm512_set1_ps(1.0f);
                        register __m512 alph,x0,k0z,inv,ear,eai,cer,cei;
                        register __m512 t0r,t0i;
                        
                        alph = _mm512_mul_ps(alp,hlf);
                        k0z  = _mm512_mul_ps(k0,z);
                        ear  = _mm512_setzero_ps();
                        inv  = _mm512_rcp14_ps(k0z);
                        eai  = k0z;
                        x0   = _mm512_mul_ps(alph,alph);
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        cmul_zmm16r4(t0r,t0i,cr,ci,&cer,&cei);
                        cer  = _mm512_mul_ps(x0,_mm512_mul_ps(cer,inv));
                        cei  = _mm512_mul_ps(x0,_mm512_mul_ps(cei,inv));
                        _mm512_storeu_ps(&xre[0] ,_mm512_mul_ps(x,cer));
                        _mm512_storeu_ps(&xim[0] ,_mm512_mul_ps(x,cei));
                        _mm512_storeu_ps(&yre[0] ,_mm512_mul_ps(y,cer));
                        _mm512_storeu_ps(&yim[0] ,_mm512_mul_ps(y,cei));
                        _mm512_storeu_ps(&zre[0] ,_mm512_mul_ps(z,cer));
                        _mm512_storeu_ps(&zim[0] ,_mm512_mul_ps(z,cei));
                 }


                     /*
                           Wide-angle cone.
                           Scattered E-field.
                           Formula 6.2-25
                       */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f6225_zmm16r4(const __m512 k0,
                                         const __m512 z,
                                         const __m512 alp,
                                         const __m512 x,
                                         const __m512 y,
                                         const __m512 z,
                                         __m512 * __restrict xre,
                                         __m512 * __restrict xim,
                                         __m512 * __restrict yre,
                                         __m512 * __restrict yim,
                                         __m512 * __restrict zre,
                                         __m512 * __restrict zim) {

                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        const __m512 _4  = _mm512_set1_ps(4.0f);
                        const __m512 cr  = _mm512_set1_ps(0.0f);
                        const __m512 ci  = _mm512_set1_ps(1.0f);
                        register __m512 pima,x0,pimas,k0z,inv,ear,eai,cer,cei;
                        register __m512 t0r,t0i;
                        
                        pima = _mm512_sub_ps(pi2,alp);
                        k0z  = _mm512_mul_ps(k0,z);
                        pimas= _mm512_mul_ps(pima,pima);
                        ear  = _mm512_setzero_ps();
                        x0   = _mm512_mul_ps(_mm512_mul_ps(_4,k0z),pimas);
                        inv  = _mm512_rcp14_ps(x0);
                        eai  = k0z;
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        cmul_zmm16r4(t0r,t0i,cr,ci,&cer,&cei);
                        cer  = _mm512_mul_ps(alph,_mm512_mul_ps(cer,inv));
                        cei  = _mm512_mul_ps(alph,_mm512_mul_ps(cei,inv));
                        *xre  = _mm512_mul_ps(x,cer);
                        *xim  = _mm512_mul_ps(x,cei);
                        *yre  = _mm512_mul_ps(y,cer);
                        *yim  = _mm512_mul_ps(y,cei);
                        *zre  = _mm512_mul_ps(z,cer);
                        *zim  = _mm512_mul_ps(z,cei);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f6225_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pz,
                                           const float * __restrict __ATTR_ALIGN__(64) palp,
                                           const float * __restrict __ATTR_ALIGN__(64) px,
                                           const float * __restrict __ATTR_ALIGN__(64) py,
                                           const float * __restrict __ATTR_ALIGN__(64) pz,
                                           float * __restrict __ATTR_ALIGN__(64) xre,
                                           float * __restrict __ATTR_ALIGN__(64) xim,
                                           float * __restrict __ATTR_ALIGN__(64) yre,
                                           float * __restrict __ATTR_ALIGN__(64) yim,
                                           float * __restrict __ATTR_ALIGN__(64) zre,
                                           float * __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 z   = _mm512_load_ps(&pz[0]);
                        register __m512 alp = _mm512_load_ps(&palp[0]);
                        register __m512 x   = _mm512_load_ps(&px[0]);
                        register __m512 y   = _mm512_load_ps(&py[0]);
                        register __m512 z   = _mm512_load_ps(&pz[0]);
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        const __m512 _4  = _mm512_set1_ps(4.0f);
                        const __m512 cr  = _mm512_set1_ps(0.0f);
                        const __m512 ci  = _mm512_set1_ps(1.0f);
                        register __m512 pima,x0,pimas,k0z,inv,ear,eai,cer,cei;
                        register __m512 t0r,t0i;
                        
                        pima = _mm512_sub_ps(pi2,alp);
                        k0z  = _mm512_mul_ps(k0,z);
                        pimas= _mm512_mul_ps(pima,pima);
                        ear  = _mm512_setzero_ps();
                        x0   = _mm512_mul_ps(_mm512_mul_ps(_4,k0z),pimas);
                        inv  = _mm512_rcp14_ps(x0);
                        eai  = k0z;
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        cmul_zmm16r4(t0r,t0i,cr,ci,&cer,&cei);
                        cer  = _mm512_mul_ps(alph,_mm512_mul_ps(cer,inv));
                        cei  = _mm512_mul_ps(alph,_mm512_mul_ps(cei,inv));
                        _mm512_store_ps(&xre[0] ,_mm512_mul_ps(x,cer));
                        _mm512_store_ps(&xim[0] ,_mm512_mul_ps(x,cei));
                        _mm512_store_ps(&yre[0] ,_mm512_mul_ps(y,cer));
                        _mm512_store_ps(&yim[0] ,_mm512_mul_ps(y,cei));
                        _mm512_store_ps(&zre[0] ,_mm512_mul_ps(z,cer));
                        _mm512_store_ps(&zim[0] ,_mm512_mul_ps(z,cei));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f6225_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict  pz,
                                           const float * __restrict  palp,
                                           const float * __restrict  px,
                                           const float * __restrict  py,
                                           const float * __restrict  pz,
                                           float * __restrict  xre,
                                           float * __restrict  xim,
                                           float * __restrict  yre,
                                           float * __restrict  yim,
                                           float * __restrict  zre,
                                           float * __restrict  zim) {

                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 z   = _mm512_loadu_ps(&pz[0]);
                        register __m512 alp = _mm512_loadu_ps(&palp[0]);
                        register __m512 x   = _mm512_loadu_ps(&px[0]);
                        register __m512 y   = _mm512_loadu_ps(&py[0]);
                        register __m512 z   = _mm512_loadu_ps(&pz[0]);
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        const __m512 _4  = _mm512_set1_ps(4.0f);
                        const __m512 cr  = _mm512_set1_ps(0.0f);
                        const __m512 ci  = _mm512_set1_ps(1.0f);
                        register __m512 pima,x0,pimas,k0z,inv,ear,eai,cer,cei;
                        register __m512 t0r,t0i;
                        
                        pima = _mm512_sub_ps(pi2,alp);
                        k0z  = _mm512_mul_ps(k0,z);
                        pimas= _mm512_mul_ps(pima,pima);
                        ear  = _mm512_setzero_ps();
                        x0   = _mm512_mul_ps(_mm512_mul_ps(_4,k0z),pimas);
                        inv  = _mm512_rcp14_ps(x0);
                        eai  = k0z;
                        cexp_zmm16r4(ear,eai,&t0r,&t0i);
                        cmul_zmm16r4(t0r,t0i,cr,ci,&cer,&cei);
                        cer  = _mm512_mul_ps(alph,_mm512_mul_ps(cer,inv));
                        cei  = _mm512_mul_ps(alph,_mm512_mul_ps(cei,inv));
                        _mm512_storeu_ps(&xre[0] ,_mm512_mul_ps(x,cer));
                        _mm512_storeu_ps(&xim[0] ,_mm512_mul_ps(x,cei));
                        _mm512_storeu_ps(&yre[0] ,_mm512_mul_ps(y,cer));
                        _mm512_storeu_ps(&yim[0] ,_mm512_mul_ps(y,cei));
                        _mm512_storeu_ps(&zre[0] ,_mm512_mul_ps(z,cer));
                        _mm512_storeu_ps(&zim[0] ,_mm512_mul_ps(z,cei));
                 }


                   /*
                           Wide-angle cone.
                           Radar Cross Section (angle = 0)
                           Formula 6.2-26
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6226_zmm16r4(const __m512 gam2,
                                            const __m512 alp) {

                         const __m512 pi2   = _mm512_set1_ps(1.57079632679489661923132169164f);
                         const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                         register __m512 rcs,gam2,trm1,inv,pim,x0,pim4;
                         gam2 = _mm512_mul_ps(gam0,gam0);
                         pim  = _mm512_sub_ps(pi2,alp);
                         trm1 = _mm512_div_ps(gam2,_16pi);
                         x0   = _mm512_mul_ps(pim,pim);
                         pim4 = _mm512_mul_ps(x0,x0);
                         inv  = _mm512_rcp14_ps(pim4);
                         rcs  = _mm512_mul_ps(trm1,inv);
                         return (rcs);
                 } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6226_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam2,
                                            const float * __restrict __ATTR_ALIGN__(64) palp) {

                         register __m512 gam2   = _mm512_load_ps(&pgam2[0]);
                         register __m512 alp = _mm512_load_ps(&palp[0]);
                         const __m512 pi2   = _mm512_set1_ps(1.57079632679489661923132169164f);
                         const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                         register __m512 rcs,gam2,trm1,inv,pim,x0,pim4;
                         gam2 = _mm512_mul_ps(gam0,gam0);
                         pim  = _mm512_sub_ps(pi2,alp);
                         trm1 = _mm512_div_ps(gam2,_16pi);
                         x0   = _mm512_mul_ps(pim,pim);
                         pim4 = _mm512_mul_ps(x0,x0);
                         inv  = _mm512_rcp14_ps(pim4);
                         rcs  = _mm512_mul_ps(trm1,inv);
                         return (rcs);
                 } 


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6226_zmm16r4_u(const float * __restrict  pgam2,
                                            const float * __restrict  palp) {

                         register __m512 gam2   = _mm512_loadu_ps(&pgam2[0]);
                         register __m512 alp = _mm512_loadu_ps(&palp[0]);
                         const __m512 pi2   = _mm512_set1_ps(1.57079632679489661923132169164f);
                         const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                         register __m512 rcs,gam2,trm1,inv,pim,x0,pim4;
                         gam2 = _mm512_mul_ps(gam0,gam0);
                         pim  = _mm512_sub_ps(pi2,alp);
                         trm1 = _mm512_div_ps(gam2,_16pi);
                         x0   = _mm512_mul_ps(pim,pim);
                         pim4 = _mm512_mul_ps(x0,x0);
                         inv  = _mm512_rcp14_ps(pim4);
                         rcs  = _mm512_mul_ps(trm1,inv);
                         return (rcs);
                 } 


                    /*
                          The concave-tip axial-incidence, Physical-Optics RCS.
                          Formula 6.2-29

                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6229_zmm16r4(const __m512 k0,
                                            const __m512 alp,
                                            const __m512 R) {

                           const __m512 _4pi = _mm512_set1_ps(12.566370614359172953850573533118f);
                           const __m512 _2   = _mm512_set1_ps(2.0f);
                           const __m512 _1   = _mm512_set1_ps(1.0f);
                           register __m512 rcs,trm1,ear1,eai1,ear2,eai2;
                           register __m512 cer1,cei1,cer2,cei2,inv1,inv2;
                           register __m512 t0r,t0i,x0,x1,k0R,R2,calp,x2,x3,t1r,t1i,cabs;
                           R2   = _m512_mul_ps(R,R);
                           x1   = _mm512_mul_ps(k0,k0);
                           calp = xcosf(alp);
                           k0R  = _mm512_mul_ps(k0,R);
                           x0   = _mm512_mul_ps(R2,R2);
                           ear1 = _mm512_setzero_ps();
                           trm1 = _mm512_mul_ps(_mm512_mul_ps(_4pi,x1),x0);
                           eai1 = _mm512_add_ps(k0R,k0R);
                           x0   = _mm512_mul_ps(eai1,eai1);
                           x1   = _mm512_sub_ps(eai1,_1);
                           inv1 = _mm512_rcp14_ps(x0);
                           cexp_zmm16r4(ear1,eai1,t0r,t0i);
                           x2   = _mm512_sub_ps(calp,_1);
                           cmul_zmm16r4(ear1,x1,t0r,t0i,&cer1,&cei1);
                           cer1 = _mm512_mul_ps(cer1,inv1);
                           cei1 = _mm512_mul_ps(cei1,inv1);
                           ear2 = ear1;
                           eai2 = _mm512_mul_ps(eai1,calp);
                           inv2 = _mm512_rcp14_ps(_mm512_mul_ps(eai2,eai2));
                           x3   = _mm512_mul_ps(eai1,x2); 
                           cexp_zmm16r4(ear2,eai2,t0r,t0i);
                           cmul_zmm16r4(ear2,x3,t0r,t0i,&cer2,&cei2);
                           cer2 = _mm512_mul_ps(cer2,inv2);
                           cei2 = _mm512_mul_ps(cei2,inv2);
                           t1r  = _mm512_sub_ps(cer1,cer2);
                           t1i  = _mm512_sub_ps(cei1,cei2);
                           cabs = cabs_zmm16r4(t1r,t1i);
                           rcs  = _mm512_mul_ps(trm1,cabs);
                           return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6229_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) palp,
                                              const float * __restrict __ATTR_ALIGN__(64) pR) {

                           register __m512 k0     = _mm512_load_ps(&pk0[0]);
                           register __m512 alp    = _mm512_load_ps(&palp[0]);
                           register __m512 R      = _mm512_load_ps(&pR[0]);
                           const __m512 _4pi = _mm512_set1_ps(12.566370614359172953850573533118f);
                           const __m512 _2   = _mm512_set1_ps(2.0f);
                           const __m512 _1   = _mm512_set1_ps(1.0f);
                           register __m512 rcs,trm1,ear1,eai1,ear2,eai2;
                           register __m512 cer1,cei1,cer2,cei2,inv1,inv2;
                           register __m512 t0r,t0i,x0,x1,k0R,k0R2,R2,calp,x2,x3,t1r,t1i,cabs;
                           R2   = _m512_mul_ps(R,R);
                           x1   = _mm512_mul_ps(k0,k0);
                           calp = xcosf(alp);
                           k0R  = _mm512_mul_ps(k0,R);
                           x0   = _mm512_mul_ps(R2,R2);
                           ear1 = _mm512_setzero_ps();
                           trm1 = _mm512_mul_ps(_mm512_mul_ps(_4pi,x1),x0);
                           eai1 = _mm512_add_ps(k0R,k0R);
                           x0   = _mm512_mul_ps(eai1,eai1);
                           x1   = _mm512_sub_ps(eai1,_1);
                           inv1 = _mm512_rcp14_ps(x0);
                           cexp_zmm16r4(ear1,eai1,t0r,t0i);
                           x2   = _mm512_sub_ps(calp,_1);
                           cmul_zmm16r4(ear1,x1,t0r,t0i,&cer1,&cei1);
                           cer1 = _mm512_mul_ps(cer1,inv1);
                           cei1 = _mm512_mul_ps(cei1,inv1);
                           ear2 = ear1;
                           eai2 = _mm512_mul_ps(eai1,calp);
                           inv2 = _mm512_rcp14_ps(_mm512_mul_ps(eai2,eai2));
                           x3   = _mm512_mul_ps(eai1,x2); 
                           cexp_zmm16r4(ear2,eai2,t0r,t0i);
                           cmul_zmm16r4(ear2,x3,t0r,t0i,&cer2,&cei2);
                           cer2 = _mm512_mul_ps(cer2,inv2);
                           cei2 = _mm512_mul_ps(cei2,inv2);
                           t1r  = _mm512_sub_ps(cer1,cer2);
                           t1i  = _mm512_sub_ps(cei1,cei2);
                           cabs = cabs_zmm16r4(t1r,t1i);
                           rcs  = _mm512_mul_ps(trm1,cabs);
                           return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6229_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  palp,
                                              const float * __restrict  pR) {

                           register __m512 k0     = _mm512_loadu_ps(&pk0[0]);
                           register __m512 alp    = _mm512_loadu_ps(&palp[0]);
                           register __m512 R      = _mm512_loadu_ps(&pR[0]);
                           const __m512 _4pi = _mm512_set1_ps(12.566370614359172953850573533118f);
                           const __m512 _2   = _mm512_set1_ps(2.0f);
                           const __m512 _1   = _mm512_set1_ps(1.0f);
                           register __m512 rcs,trm1,ear1,eai1,ear2,eai2;
                           register __m512 cer1,cei1,cer2,cei2,inv1,inv2;
                           register __m512 t0r,t0i,x0,x1,k0R,k0R2,R2,calp,x2,x3,t1r,t1i,cabs;
                           R2   = _m512_mul_ps(R,R);
                           x1   = _mm512_mul_ps(k0,k0);
                           calp = xcosf(alp);
                           k0R  = _mm512_mul_ps(k0,R);
                           x0   = _mm512_mul_ps(R2,R2);
                           ear1 = _mm512_setzero_ps();
                           trm1 = _mm512_mul_ps(_mm512_mul_ps(_4pi,x1),x0);
                           eai1 = _mm512_add_ps(k0R,k0R);
                           x0   = _mm512_mul_ps(eai1,eai1);
                           x1   = _mm512_sub_ps(eai1,_1);
                           inv1 = _mm512_rcp14_ps(x0);
                           cexp_zmm16r4(ear1,eai1,t0r,t0i);
                           x2   = _mm512_sub_ps(calp,_1);
                           cmul_zmm16r4(ear1,x1,t0r,t0i,&cer1,&cei1);
                           cer1 = _mm512_mul_ps(cer1,inv1);
                           cei1 = _mm512_mul_ps(cei1,inv1);
                           ear2 = ear1;
                           eai2 = _mm512_mul_ps(eai1,calp);
                           inv2 = _mm512_rcp14_ps(_mm512_mul_ps(eai2,eai2));
                           x3   = _mm512_mul_ps(eai1,x2); 
                           cexp_zmm16r4(ear2,eai2,t0r,t0i);
                           cmul_zmm16r4(ear2,x3,t0r,t0i,&cer2,&cei2);
                           cer2 = _mm512_mul_ps(cer2,inv2);
                           cei2 = _mm512_mul_ps(cei2,inv2);
                           t1r  = _mm512_sub_ps(cer1,cer2);
                           t1i  = _mm512_sub_ps(cei1,cei2);
                           cabs = cabs_zmm16r4(t1r,t1i);
                           rcs  = _mm512_mul_ps(trm1,cabs);
                           return (rcs);
                 }


                    /*
                           RCS of pointed cone and flat plate of radius b.
                           Formula 6.2-30
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6230_zmm16r4(const __m512 b,
                                            const __m512 k0) {

                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 b4,k02,k0b2,b2;
                          register __m512 trm1,trm1,rcs,x0;
                          b2  = _mm512_mul_ps(b,b);
                          k02 = _mm512_mul_ps(k0,k0);
                          b4  = _mm512_mul_ps(b2,b2);
                          k0b2= _mm512_mul_ps(k02,k02);
                          trm2= _mm512_mul_ps(k0b2,_mm512_mul_ps(pi,b2));
                          x0  = _mm512_mul_ps(_4,k02);
                          trm1= _mm512_div_ps(_mm512_mul_ps(pi,b4),x0);
                          rcs = _mm512_add_ps(trm1,trm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6230_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0) {

                          register __m512 b  = _mm512_load_ps(&pb[0]);
                          register __m512 k0 = _mm512_load_ps(&pk0[0]);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 b4,k02,k0b2,b2;
                          register __m512 trm1,trm1,rcs,x0;
                          b2  = _mm512_mul_ps(b,b);
                          k02 = _mm512_mul_ps(k0,k0);
                          b4  = _mm512_mul_ps(b2,b2);
                          k0b2= _mm512_mul_ps(k02,k02);
                          trm2= _mm512_mul_ps(k0b2,_mm512_mul_ps(pi,b2));
                          x0  = _mm512_mul_ps(_4,k02);
                          trm1= _mm512_div_ps(_mm512_mul_ps(pi,b4),x0);
                          rcs = _mm512_add_ps(trm1,trm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6230_zmm16r4_u(const float * __restrict  pb,
                                              const float * __restrict  pk0) {

                          register __m512 b  = _mm512_loadu_ps(&pb[0]);
                          register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 b4,k02,k0b2,b2;
                          register __m512 trm1,trm1,rcs,x0;
                          b2  = _mm512_mul_ps(b,b);
                          k02 = _mm512_mul_ps(k0,k0);
                          b4  = _mm512_mul_ps(b2,b2);
                          k0b2= _mm512_mul_ps(k02,k02);
                          trm2= _mm512_mul_ps(k0b2,_mm512_mul_ps(pi,b2));
                          x0  = _mm512_mul_ps(_4,k02);
                          trm1= _mm512_div_ps(_mm512_mul_ps(pi,b4),x0);
                          rcs = _mm512_add_ps(trm1,trm2);
                          return (rcs);
                }


                  /*
                       Physical-Optics approximation of rounded-tip cone,
                       for axial incidence.
                       Formula 6.2-31
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6231_zmm16r4(const __m512 alp,
                                            const __m512 k0,
                                            const __m512 b) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 b2,k0b,_2k0b,k0b2,arg;
                          register __m512 sarg,carg,salp,calp,calp2,calp4;
                          register __m512 rcs,trm1,trm2,trm3,x0,x1;
                          b2   = _mm512_mul_ps(b,b);
                          calp = xcosf(alp);
                          k0b  = _mm512_mul_ps(k0,b);
                          calp2= _mm512_mul_ps(calp,calp);
                          _2k0b= _mm512_add_ps(k0b,k0b);
                          calp4= _mm512_mul_ps(calp2,calp2);
                          k0b2 = _mm512_mul_ps(k0b,k0b);
                          salp = xsinf(alp);
                          arg  = _mm512_mul_ps(_2k0b,_mm512_sub_ps(_1,salp));
                          x0   = _mm512_mul_ps(k0b,calp2);
                          sarg = xsinf(arg);
                          x1   = _mm512_add_ps(_1,calp4);
                          trm1 = _mm512_div_ps(sarg,x0);
                          x0   = _mm512_mul_ps(_mm512_mul_ps(_4,k0b2),calp4);
                          trm2 = _mm512_div_ps(x1,x0);
                          carg = xcosf(arg);
                          x0   = _mm512_mul_ps(_mm512_add_ps(k0b2,k0b2),calp2);
                          trm3 = _mm512_div_ps(carg,x0);
                          x1   = _mm512_mul_ps(pi,b2);
                          x0   = _mm512_add_ps(_mm512_sub_ps(_1,trm1),trm2);
                          rcs  = _mm512_mul_ps(x0,_mm512_sub_ps(x0,trm3));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6231_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) palp,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0,
                                            const float * __restrict __ATTR_ALIGN__(64) pb) {

                          register __m512 alp= _mm512_load_ps(&palp[0]);
                          register __m512 b  = _mm512_load_ps(&pb[0]);
                          register __m512 k0 = _mm512_load_ps(&pk0[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 b2,k0b,_2k0b,k0b2,arg;
                          register __m512 sarg,carg,salp,calp,calp2,calp4;
                          register __m512 rcs,trm1,trm2,trm3,x0,x1;
                          b2   = _mm512_mul_ps(b,b);
                          calp = xcosf(alp);
                          k0b  = _mm512_mul_ps(k0,b);
                          calp2= _mm512_mul_ps(calp,calp);
                          _2k0b= _mm512_add_ps(k0b,k0b);
                          calp4= _mm512_mul_ps(calp2,calp2);
                          k0b2 = _mm512_mul_ps(k0b,k0b);
                          salp = xsinf(alp);
                          arg  = _mm512_mul_ps(_2k0b,_mm512_sub_ps(_1,salp));
                          x0   = _mm512_mul_ps(k0b,calp2);
                          sarg = xsinf(arg);
                          x1   = _mm512_add_ps(_1,calp4);
                          trm1 = _mm512_div_ps(sarg,x0);
                          x0   = _mm512_mul_ps(_mm512_mul_ps(_4,k0b2),calp4);
                          trm2 = _mm512_div_ps(x1,x0);
                          carg = xcosf(arg);
                          x0   = _mm512_mul_ps(_mm512_add_ps(k0b2,k0b2),calp2);
                          trm3 = _mm512_div_ps(carg,x0);
                          x1   = _mm512_mul_ps(pi,b2);
                          x0   = _mm512_add_ps(_mm512_sub_ps(_1,trm1),trm2);
                          rcs  = _mm512_mul_ps(x0,_mm512_sub_ps(x0,trm3));
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6231_zmm16r4_u(const float * __restrict  palp,
                                            const float * __restrict  pk0,
                                            const float * __restrict  pb) {

                          register __m512 alp= _mm512_loadu_ps(&palp[0]);
                          register __m512 b  = _mm512_loadu_ps(&pb[0]);
                          register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 b2,k0b,_2k0b,k0b2,arg;
                          register __m512 sarg,carg,salp,calp,calp2,calp4;
                          register __m512 rcs,trm1,trm2,trm3,x0,x1;
                          b2   = _mm512_mul_ps(b,b);
                          calp = xcosf(alp);
                          k0b  = _mm512_mul_ps(k0,b);
                          calp2= _mm512_mul_ps(calp,calp);
                          _2k0b= _mm512_add_ps(k0b,k0b);
                          calp4= _mm512_mul_ps(calp2,calp2);
                          k0b2 = _mm512_mul_ps(k0b,k0b);
                          salp = xsinf(alp);
                          arg  = _mm512_mul_ps(_2k0b,_mm512_sub_ps(_1,salp));
                          x0   = _mm512_mul_ps(k0b,calp2);
                          sarg = xsinf(arg);
                          x1   = _mm512_add_ps(_1,calp4);
                          trm1 = _mm512_div_ps(sarg,x0);
                          x0   = _mm512_mul_ps(_mm512_mul_ps(_4,k0b2),calp4);
                          trm2 = _mm512_div_ps(x1,x0);
                          carg = xcosf(arg);
                          x0   = _mm512_mul_ps(_mm512_add_ps(k0b2,k0b2),calp2);
                          trm3 = _mm512_div_ps(carg,x0);
                          x1   = _mm512_mul_ps(pi,b2);
                          x0   = _mm512_add_ps(_mm512_sub_ps(_1,trm1),trm2);
                          rcs  = _mm512_mul_ps(x0,_mm512_sub_ps(x0,trm3));
                          return (rcs);
                 }


                   /*
                          Not nose-one incidence, for theta << alpha -- 
                          PO approximation for rounded-tip cones.
                          Formula 6.2-23
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6223_zmm16r4(const __m512 tht,
                                            const __m512 k0,
                                            const __m512 b,
                                            const __m512 a,
                                            const __m512 alp) {

                          const __m512 hlf= _mm512_set1_ps(0.5f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 _3 = _mm512_set1_ps(3.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          const __m512 _6 = _mm512_set1_ps(6.0f);
                          const __m512 _8 = _mm512_set1_ps(8.0f);
                          const __m512 _13= _mm512_set1_ps(13.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 k0b,salp,arg,_2k0,b2,sarg,carg,ctht,a2,salp1,tht3;
                          register __m512 trm1,trm2,trm3;
                          register __m512 rcs,A1,A2,A3;
                          __m512 alp2,tht2,alp4,tht4,k0b2,k0b3,k0b4,a2,k0bt,k0a;
                          __m512 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
                          __m512 t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;
                          __m512 x0,x1,x2,x3,x4;
                          k0b  = _mm512_mul_ps(k0,b);
                          tht2 = _mm512_mul_ps(tht,tht);
                          _2k0 = _mm512_add_ps(k0,k0);
                          salp = xsinf(alp);
                          k0b2 = _mm512_mul_ps(k0b,k0b);
                          ctht = xcosf(tht);
                          a2   = _mm512_mul_ps(a,a);
                          salp1= _mm512_sub_ps(_1,salp);
                          k0a  = _mm512_mul_ps(k0,a);
                          alp2 = _mm512_mul_ps(alp,alp);
                          tht3 = _mm512_mul_ps(tht2,tht);
                          arg  = _mm512_mul_ps(_2k0,_mm512_mul_ps(calp,salp1));
                          tht4 = _mm512_mul_ps(tht2,tht2);
                          alp4 = _mm512_mul_ps(alp2,alp2);
                          k0b3 = _mm512_mul_ps(k0b2,k0b);
                          k0b4 = _mm512_mul_ps(k0b2,k0b2);
                          trm1 = _mm512_div_ps(_mm512_add_ps(_1,tht2),
                                               _mm512_mul_ps(_4,k0b2));
                          a2   = _mm512_mul_ps(a,a);
                          k0bt = _mm512_mul_ps(k0b,tht);
                          t0   = _mm512_add_ps(alp2,alp2); //2alp^2
                          t1   = _mm512_add_ps(tht2,tht2); //2tht^2
                          t2   = _mm512_mul_ps(alp2,tht2);//alp2*tht2
                          t3   = _mm512_mul_ps(tht4,hlf); //tht4*0.5
                          t4   = _mm512_mul_ps(_mm512_mul_ps(_2,alp4),tht2);//2*alp4*tht2
                          t5   = _mm512_mul_ps(_4,k0b2);//4*k0b2
                          t6   = _mm512_mul_ps(_mm512_mul_ps(_2,k0b2),tht2);//2*k0b2*tht2
                          t7   = _mm512_mul_ps(_mm512_mul_ps(_8,k0b3),tht2);//8*k0b3*tht2
                          t8   = _mm512_mul_ps(k0b2,tht4);//k0b2*tht4
                          t9   = _mm512_mul_ps(_mm512_mul_ps(_6,k0b2),
                                               _mm512_mul_ps(a2,tht2));//6*k0b2*a2*tht2
                          t10  = _mm512_mul_ps(_mm512_mul_ps(_8,k0b3),tht4);//8*k0b3*tht4
                          t11  = _mm512_mul_ps(_mm512_mul_ps(_13,k0b4),tht4);//13*k0b4*tht4
                          sarg = xsinf(arg);
                          x0   = _mm512_sub_ps(_mm512_add_ps(_2,t0),
                                               _mm512_add_ps(t1,alp4));
                          x1   = _mm512_add_ps(_mm512_add_ps(t2,t3),
                                               _mm512_add_ps(t4,t5));
                          x2   = _mm512_add_ps(_mm512_sub_ps(t6,t7),
                                               _mm512_add_ps(t8,t9));
                          x3   = _mm512_add_ps(t10,t11);
                          x4   = _mm512_sub_ps(x0,x1);
                          carg = xcosf(arg);
                          t9   = _mm512_mul_ps(_mm512_mul_ps(_6,k0b2),tht2);
                          A1   = _mm512_sub_ps(_mm512_sub_ps(x0,x1),
                                               _mm512_add_ps(x2,x3));
                          t7   = _mm512_mul_ps(_mm512_mul_ps(_8,k0b4),tht3);
                          x0   = _mm512_sub_ps(_mm512_sub_ps(n2,t0),
                                               _mm512_add_ps(t1,t2));
                          t8   = _mm512_mul_ps(_mm512_mul_ps(_3,k0b2),tht4);
                          x1   = _mm512_add_ps(_mm512_sub_ps(t3,t9),
                                               _mm512_add_ps(t7,t8));
                          A2   = _mm512_sub_ps(x0,x1);
                          t11  = _mm512_mul_ps(k0bt,k0bt);
                          t10  = _mm512_mul_ps(t11,k0bt);
                          t9   = _mm512_mul_ps(k0b,tht2);
                          t8   = _mm512_mul_ps(_3,_mm512_mul_ps(k0a,k0a));
                          x0   = _mm512_sub_ps(_mm512_add_ps(_1,alp2),
                                               _mm512_add_ps(t2,t8));
                          x1   = _mm512_sub_ps(k0b,_mm512_sub_ps(t9,t11));
                          A3   = _mm512_mul_ps(_mm512_mul_ps(n4,x0),
                                               _mm512_sub_ps(x1,_mm512_mul_ps(_4,t10)));
                          trm2 = _mm512_fmadd_ps(carg,A2,A1);
                          trm3 = _mm512_fmadd_ps(sarg,A3,trm2);
                          rcs  = _mm512_mul_ps(trm1,trm3);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6223_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0,
                                            const float * __restrict __ATTR_ALIGN__(64) pb,
                                            const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) palp) {

                          register __m512 tht  = _mm512_load_ps(&ptht[0]);
                          register __m512 k0   = _mm512_load_ps(&pk0[0]);
                          register __m512 b    = _mm512_load_ps(&pb[0]);
                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 alp  = _mm512_load_ps(&palp[0]);
                          const __m512 hlf= _mm512_set1_ps(0.5f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 _3 = _mm512_set1_ps(3.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          const __m512 _6 = _mm512_set1_ps(6.0f);
                          const __m512 _8 = _mm512_set1_ps(8.0f);
                          const __m512 _13= _mm512_set1_ps(13.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 k0b,salp,arg,_2k0,b2,sarg,carg,ctht,a2,salp1,tht3;
                          register __m512 trm1,trm2,trm3;
                          register __m512 rcs,A1,A2,A3;
                          __m512 alp2,tht2,alp4,tht4,k0b2,k0b3,k0b4,a2,k0bt,k0a;
                          __m512 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
                          __m512 t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;
                          __m512 x0,x1,x2,x3,x4;
                          k0b  = _mm512_mul_ps(k0,b);
                          tht2 = _mm512_mul_ps(tht,tht);
                          _2k0 = _mm512_add_ps(k0,k0);
                          salp = xsinf(alp);
                          k0b2 = _mm512_mul_ps(k0b,k0b);
                          ctht = xcosf(tht);
                          a2   = _mm512_mul_ps(a,a);
                          salp1= _mm512_sub_ps(_1,salp);
                          k0a  = _mm512_mul_ps(k0,a);
                          alp2 = _mm512_mul_ps(alp,alp);
                          tht3 = _mm512_mul_ps(tht2,tht);
                          arg  = _mm512_mul_ps(_2k0,_mm512_mul_ps(calp,salp1));
                          tht4 = _mm512_mul_ps(tht2,tht2);
                          alp4 = _mm512_mul_ps(alp2,alp2);
                          k0b3 = _mm512_mul_ps(k0b2,k0b);
                          k0b4 = _mm512_mul_ps(k0b2,k0b2);
                          trm1 = _mm512_div_ps(_mm512_add_ps(_1,tht2),
                                               _mm512_mul_ps(_4,k0b2));
                          a2   = _mm512_mul_ps(a,a);
                          k0bt = _mm512_mul_ps(k0b,tht);
                          t0   = _mm512_add_ps(alp2,alp2); //2alp^2
                          t1   = _mm512_add_ps(tht2,tht2); //2tht^2
                          t2   = _mm512_mul_ps(alp2,tht2);//alp2*tht2
                          t3   = _mm512_mul_ps(tht4,hlf); //tht4*0.5
                          t4   = _mm512_mul_ps(_mm512_mul_ps(_2,alp4),tht2);//2*alp4*tht2
                          t5   = _mm512_mul_ps(_4,k0b2);//4*k0b2
                          t6   = _mm512_mul_ps(_mm512_mul_ps(_2,k0b2),tht2);//2*k0b2*tht2
                          t7   = _mm512_mul_ps(_mm512_mul_ps(_8,k0b3),tht2);//8*k0b3*tht2
                          t8   = _mm512_mul_ps(k0b2,tht4);//k0b2*tht4
                          t9   = _mm512_mul_ps(_mm512_mul_ps(_6,k0b2),
                                               _mm512_mul_ps(a2,tht2));//6*k0b2*a2*tht2
                          t10  = _mm512_mul_ps(_mm512_mul_ps(_8,k0b3),tht4);//8*k0b3*tht4
                          t11  = _mm512_mul_ps(_mm512_mul_ps(_13,k0b4),tht4);//13*k0b4*tht4
                          sarg = xsinf(arg);
                          x0   = _mm512_sub_ps(_mm512_add_ps(_2,t0),
                                               _mm512_add_ps(t1,alp4));
                          x1   = _mm512_add_ps(_mm512_add_ps(t2,t3),
                                               _mm512_add_ps(t4,t5));
                          x2   = _mm512_add_ps(_mm512_sub_ps(t6,t7),
                                               _mm512_add_ps(t8,t9));
                          x3   = _mm512_add_ps(t10,t11);
                          x4   = _mm512_sub_ps(x0,x1);
                          carg = xcosf(arg);
                          t9   = _mm512_mul_ps(_mm512_mul_ps(_6,k0b2),tht2);
                          A1   = _mm512_sub_ps(_mm512_sub_ps(x0,x1),
                                               _mm512_add_ps(x2,x3));
                          t7   = _mm512_mul_ps(_mm512_mul_ps(_8,k0b4),tht3);
                          x0   = _mm512_sub_ps(_mm512_sub_ps(n2,t0),
                                               _mm512_add_ps(t1,t2));
                          t8   = _mm512_mul_ps(_mm512_mul_ps(_3,k0b2),tht4);
                          x1   = _mm512_add_ps(_mm512_sub_ps(t3,t9),
                                               _mm512_add_ps(t7,t8));
                          A2   = _mm512_sub_ps(x0,x1);
                          t11  = _mm512_mul_ps(k0bt,k0bt);
                          t10  = _mm512_mul_ps(t11,k0bt);
                          t9   = _mm512_mul_ps(k0b,tht2);
                          t8   = _mm512_mul_ps(_3,_mm512_mul_ps(k0a,k0a));
                          x0   = _mm512_sub_ps(_mm512_add_ps(_1,alp2),
                                               _mm512_add_ps(t2,t8));
                          x1   = _mm512_sub_ps(k0b,_mm512_sub_ps(t9,t11));
                          A3   = _mm512_mul_ps(_mm512_mul_ps(n4,x0),
                                               _mm512_sub_ps(x1,_mm512_mul_ps(_4,t10)));
                          trm2 = _mm512_fmadd_ps(carg,A2,A1);
                          trm3 = _mm512_fmadd_ps(sarg,A3,trm2);
                          rcs  = _mm512_mul_ps(trm1,trm3);
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6223_zmm16r4_u(const float * __restrict  ptht,
                                            const float * __restrict  pk0,
                                            const float * __restrict  pb,
                                            const float * __restrict  pa,
                                            const float * __restrict  palp) {

                          register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                          register __m512 b    = _mm512_loadu_ps(&pb[0]);
                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 alp  = _mm512_loadu_ps(&palp[0]);
                          const __m512 hlf= _mm512_set1_ps(0.5f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 _3 = _mm512_set1_ps(3.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          const __m512 _6 = _mm512_set1_ps(6.0f);
                          const __m512 _8 = _mm512_set1_ps(8.0f);
                          const __m512 _13= _mm512_set1_ps(13.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 k0b,salp,arg,_2k0,b2,sarg,carg,ctht,a2,salp1,tht3;
                          register __m512 trm1,trm2,trm3;
                          register __m512 rcs,A1,A2,A3;
                          __m512 alp2,tht2,alp4,tht4,k0b2,k0b3,k0b4,a2,k0bt,k0a;
                          __m512 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
                          __m512 t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;
                          __m512 x0,x1,x2,x3,x4;
                          k0b  = _mm512_mul_ps(k0,b);
                          tht2 = _mm512_mul_ps(tht,tht);
                          _2k0 = _mm512_add_ps(k0,k0);
                          salp = xsinf(alp);
                          k0b2 = _mm512_mul_ps(k0b,k0b);
                          ctht = xcosf(tht);
                          a2   = _mm512_mul_ps(a,a);
                          salp1= _mm512_sub_ps(_1,salp);
                          k0a  = _mm512_mul_ps(k0,a);
                          alp2 = _mm512_mul_ps(alp,alp);
                          tht3 = _mm512_mul_ps(tht2,tht);
                          arg  = _mm512_mul_ps(_2k0,_mm512_mul_ps(calp,salp1));
                          tht4 = _mm512_mul_ps(tht2,tht2);
                          alp4 = _mm512_mul_ps(alp2,alp2);
                          k0b3 = _mm512_mul_ps(k0b2,k0b);
                          k0b4 = _mm512_mul_ps(k0b2,k0b2);
                          trm1 = _mm512_div_ps(_mm512_add_ps(_1,tht2),
                                               _mm512_mul_ps(_4,k0b2));
                          a2   = _mm512_mul_ps(a,a);
                          k0bt = _mm512_mul_ps(k0b,tht);
                          t0   = _mm512_add_ps(alp2,alp2); //2alp^2
                          t1   = _mm512_add_ps(tht2,tht2); //2tht^2
                          t2   = _mm512_mul_ps(alp2,tht2);//alp2*tht2
                          t3   = _mm512_mul_ps(tht4,hlf); //tht4*0.5
                          t4   = _mm512_mul_ps(_mm512_mul_ps(_2,alp4),tht2);//2*alp4*tht2
                          t5   = _mm512_mul_ps(_4,k0b2);//4*k0b2
                          t6   = _mm512_mul_ps(_mm512_mul_ps(_2,k0b2),tht2);//2*k0b2*tht2
                          t7   = _mm512_mul_ps(_mm512_mul_ps(_8,k0b3),tht2);//8*k0b3*tht2
                          t8   = _mm512_mul_ps(k0b2,tht4);//k0b2*tht4
                          t9   = _mm512_mul_ps(_mm512_mul_ps(_6,k0b2),
                                               _mm512_mul_ps(a2,tht2));//6*k0b2*a2*tht2
                          t10  = _mm512_mul_ps(_mm512_mul_ps(_8,k0b3),tht4);//8*k0b3*tht4
                          t11  = _mm512_mul_ps(_mm512_mul_ps(_13,k0b4),tht4);//13*k0b4*tht4
                          sarg = xsinf(arg);
                          x0   = _mm512_sub_ps(_mm512_add_ps(_2,t0),
                                               _mm512_add_ps(t1,alp4));
                          x1   = _mm512_add_ps(_mm512_add_ps(t2,t3),
                                               _mm512_add_ps(t4,t5));
                          x2   = _mm512_add_ps(_mm512_sub_ps(t6,t7),
                                               _mm512_add_ps(t8,t9));
                          x3   = _mm512_add_ps(t10,t11);
                          x4   = _mm512_sub_ps(x0,x1);
                          carg = xcosf(arg);
                          t9   = _mm512_mul_ps(_mm512_mul_ps(_6,k0b2),tht2);
                          A1   = _mm512_sub_ps(_mm512_sub_ps(x0,x1),
                                               _mm512_add_ps(x2,x3));
                          t7   = _mm512_mul_ps(_mm512_mul_ps(_8,k0b4),tht3);
                          x0   = _mm512_sub_ps(_mm512_sub_ps(n2,t0),
                                               _mm512_add_ps(t1,t2));
                          t8   = _mm512_mul_ps(_mm512_mul_ps(_3,k0b2),tht4);
                          x1   = _mm512_add_ps(_mm512_sub_ps(t3,t9),
                                               _mm512_add_ps(t7,t8));
                          A2   = _mm512_sub_ps(x0,x1);
                          t11  = _mm512_mul_ps(k0bt,k0bt);
                          t10  = _mm512_mul_ps(t11,k0bt);
                          t9   = _mm512_mul_ps(k0b,tht2);
                          t8   = _mm512_mul_ps(_3,_mm512_mul_ps(k0a,k0a));
                          x0   = _mm512_sub_ps(_mm512_add_ps(_1,alp2),
                                               _mm512_add_ps(t2,t8));
                          x1   = _mm512_sub_ps(k0b,_mm512_sub_ps(t9,t11));
                          A3   = _mm512_mul_ps(_mm512_mul_ps(n4,x0),
                                               _mm512_sub_ps(x1,_mm512_mul_ps(_4,t10)));
                          trm2 = _mm512_fmadd_ps(carg,A2,A1);
                          trm3 = _mm512_fmadd_ps(sarg,A3,trm2);
                          rcs  = _mm512_mul_ps(trm1,trm3);
                          return (rcs);
                 }


                  /*
                         Finite cones, approximated solutions.
                         Rayleigh region.
                         RCS.
                         Formula 6.3-4
                     */

#include "GMS_simd_utils.hpp"   

                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f634_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 h) {

                          using namespace gms::math;
                          const __m512 c0 = _mm512_set1_ps(1.396263401595463661538952614791f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 k04,a4,h2,_4a,arg,earg,pih,x0,x1,nh;
                          register __m512 rcs,trm1,trm2;
                          h2  = _mm512_mul_ps(h,h);
                          x0  = _mm512_mul_ps(k0,k0);
                          x1  = _mm512_mul_ps(a,a);
                          _4a = _mm512_mul_ps(_4,a);
                          pih = _mm512_mul_ps(_mm512_set1_ps(3.14159265358979323846264338328f),h);
                          nh  = negate_zmm16r4(h);
                          k04 = _mm512_mul_ps(x0,x0);
                          arg = _mm512_div_ps(nh,_4a);
                          a4  = _mm512_mul_ps(x1,x1);
                          earg= xexpf(arg);
                          trm1= _mm512_mul_ps(_mm512_mul_ps(k04,a4),h2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(_4,earg),pih);
                          trm1= _mm512_mul_ps(c0,trm1);
                          x1  = _mm512_add_ps(_1,x0);
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f634_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ph) {

                          using namespace gms::math;
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 h   = _mm512_load_ps(&ph[0]);
                          const __m512 c0 = _mm512_set1_ps(1.396263401595463661538952614791f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 k04,a4,h2,_4a,arg,earg,pih,x0,x1,nh;
                          register __m512 rcs,trm1,trm2;
                          h2  = _mm512_mul_ps(h,h);
                          x0  = _mm512_mul_ps(k0,k0);
                          x1  = _mm512_mul_ps(a,a);
                          _4a = _mm512_mul_ps(_4,a);
                          pih = _mm512_mul_ps(_mm512_set1_ps(3.14159265358979323846264338328f),h);
                          nh  = negate_zmm16r4(h);
                          k04 = _mm512_mul_ps(x0,x0);
                          arg = _mm512_div_ps(nh,_4a);
                          a4  = _mm512_mul_ps(x1,x1);
                          earg= xexpf(arg);
                          trm1= _mm512_mul_ps(_mm512_mul_ps(k04,a4),h2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(_4,earg),pih);
                          trm1= _mm512_mul_ps(c0,trm1);
                          x1  = _mm512_add_ps(_1,x0);
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f634_zmm16r4_u(const float * __restrict pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict ph) {

                          using namespace gms::math;
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 h   = _mm512_loadu_ps(&ph[0]);
                          const __m512 c0 = _mm512_set1_ps(1.396263401595463661538952614791f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 k04,a4,h2,_4a,arg,earg,pih,x0,x1,nh;
                          register __m512 rcs,trm1,trm2;
                          h2  = _mm512_mul_ps(h,h);
                          x0  = _mm512_mul_ps(k0,k0);
                          x1  = _mm512_mul_ps(a,a);
                          _4a = _mm512_mul_ps(_4,a);
                          pih = _mm512_mul_ps(_mm512_set1_ps(3.14159265358979323846264338328f),h);
                          nh  = negate_zmm16r4(h);
                          k04 = _mm512_mul_ps(x0,x0);
                          arg = _mm512_div_ps(nh,_4a);
                          a4  = _mm512_mul_ps(x1,x1);
                          earg= xexpf(arg);
                          trm1= _mm512_mul_ps(_mm512_mul_ps(k04,a4),h2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(_4,earg),pih);
                          trm1= _mm512_mul_ps(c0,trm1);
                          x1  = _mm512_add_ps(_1,x0);
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                 }


                   /*
                          Rayleigh RCS of cone-hemispheroid.
                          Formula 6.3-5
                     */


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f635_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 b,
                                           const __m512 h) {

                          using namespace gms::math;
                          const __m512 _4pi = _mm512_set1_ps(1.27323954473516268615107010698f);
                          const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0   = _mm512_set1_ps(0.33333333333333333333333333333f);
                          const __m512 c1   = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _4   = _mm512_set1_ps(4.0f);
                          const __m512 _2   = _mm512_set1_ps(2.0f);
                          const __m512 _1   = _mm512_set1_ps(1.0f);
                          register __m512 k04,a2,_4a,h2b,nh2b,pia2,arg,earg,x0,x1,x2;
                          register __m512 rcs,trm1,trm2,trm3;
                          x0  = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          h2b = _mm512_fmadd_ps(b,_2,h);
                          k04 = _mm512_mul_ps(x0,x0);
                          _4a = _mm512_mul_ps(_4,a);
                          trm1= _mm512_mul_ps(_4pi,k04);
                          pia2= _mm512_mul_ps(pi,a2);
                          nh2b= negate_zmm16r4(h2b);
                          x0  = _mm512_mul_ps(c0,_mm512_mul_ps(pia2,h));
                          arg = _mm512_div_ps(nh2b,_4a);
                          x1  = _mm512_mul_ps(c1,_mm512_mul_ps(pia2,b));
                          earg= xexpf(arg);
                          x2  = _mm512_add_ps(x0,x1);
                          trm2= _mm512_mul_ps(x2,x2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(pi,h2b),_4a);
                          x1  = _mm512_add_ps(_1, _mm512_div_ps(earg,x0));
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f635_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) pb,
                                             const float * __restrict __ATTR_ALIGN__(64) ph) {
                                             
                          using namespace gms::math;
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 b   = _mm512_load_ps(&pb[0]);
                          register __m512 h   = _mm512_load_ps(&ph[0]);
                          const __m512 _4pi = _mm512_set1_ps(1.27323954473516268615107010698f);
                          const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0   = _mm512_set1_ps(0.33333333333333333333333333333f);
                          const __m512 c1   = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _4   = _mm512_set1_ps(4.0f);
                          const __m512 _2   = _mm512_set1_ps(2.0f);
                          const __m512 _1   = _mm512_set1_ps(1.0f);
                          register __m512 k04,a2,_4a,h2b,nh2b,pia2,arg,earg,x0,x1,x2;
                          register __m512 rcs,trm1,trm2,trm3;
                          x0  = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          h2b = _mm512_fmadd_ps(b,_2,h);
                          k04 = _mm512_mul_ps(x0,x0);
                          _4a = _mm512_mul_ps(_4,a);
                          trm1= _mm512_mul_ps(_4pi,k04);
                          pia2= _mm512_mul_ps(pi,a2);
                          nh2b= negate_zmm16r4(h2b);
                          x0  = _mm512_mul_ps(c0,_mm512_mul_ps(pia2,h));
                          arg = _mm512_div_ps(nh2b,_4a);
                          x1  = _mm512_mul_ps(c1,_mm512_mul_ps(pia2,b));
                          earg= xexpf(arg);
                          x2  = _mm512_add_ps(x0,x1);
                          trm2= _mm512_mul_ps(x2,x2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(pi,h2b),_4a);
                          x1  = _mm512_add_ps(_1, _mm512_div_ps(earg,x0));
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f635_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  pb,
                                             const float * __restrict  ph) {
                                             
                          using namespace gms::math;
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 b   = _mm512_loadu_ps(&pb[0]);
                          register __m512 h   = _mm512_loadu_ps(&ph[0]);
                          const __m512 _4pi = _mm512_set1_ps(1.27323954473516268615107010698f);
                          const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0   = _mm512_set1_ps(0.33333333333333333333333333333f);
                          const __m512 c1   = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _4   = _mm512_set1_ps(4.0f);
                          const __m512 _2   = _mm512_set1_ps(2.0f);
                          const __m512 _1   = _mm512_set1_ps(1.0f);
                          register __m512 k04,a2,_4a,h2b,nh2b,pia2,arg,earg,x0,x1,x2;
                          register __m512 rcs,trm1,trm2,trm3;
                          x0  = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          h2b = _mm512_fmadd_ps(b,_2,h);
                          k04 = _mm512_mul_ps(x0,x0);
                          _4a = _mm512_mul_ps(_4,a);
                          trm1= _mm512_mul_ps(_4pi,k04);
                          pia2= _mm512_mul_ps(pi,a2);
                          nh2b= negate_zmm16r4(h2b);
                          x0  = _mm512_mul_ps(c0,_mm512_mul_ps(pia2,h));
                          arg = _mm512_div_ps(nh2b,_4a);
                          x1  = _mm512_mul_ps(c1,_mm512_mul_ps(pia2,b));
                          earg= xexpf(arg);
                          x2  = _mm512_add_ps(x0,x1);
                          trm2= _mm512_mul_ps(x2,x2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(pi,h2b),_4a);
                          x1  = _mm512_add_ps(_1, _mm512_div_ps(earg,x0));
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                }


                  /*
                          Rayleigh RCS of cone-hemispheroid (b == a)
                          Formula 6.3-6
                      */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f636_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 h) {

                          using namespace gms::math;
                          const __m512 _4pi = _mm512_set1_ps(1.27323954473516268615107010698f);
                          const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0   = _mm512_set1_ps(0.33333333333333333333333333333f);
                          const __m512 c1   = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _4   = _mm512_set1_ps(4.0f);
                          const __m512 _2   = _mm512_set1_ps(2.0f);
                          const __m512 _1   = _mm512_set1_ps(1.0f);
                          register __m512 k04,a2,_4a,a3,h2a,nh2a,pia2,pia3,arg,earg,x0,x1,x2;
                          register __m512 rcs,trm1,trm2,trm3;
                          x0  = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          h2a = _mm512_fmadd_ps(a,_2,h);
                          a3  = _mm512_mul_ps(a2,a);
                          k04 = _mm512_mul_ps(x0,x0);
                          _4a = _mm512_mul_ps(_4,a);
                          trm1= _mm512_mul_ps(_4pi,k04);
                          pia2= _mm512_mul_ps(pi,a2);
                          nh2a= negate_zmm16r4(h2a);
                          x0  = _mm512_mul_ps(c0,_mm512_mul_ps(pia2,h));
                          pia3= _mm512_mul_ps(pi,a3);
                          arg = _mm512_div_ps(nh2a,_4a);
                          x1  = _mm512_mul_ps(c1,pia3);
                          earg= xexpf(arg);
                          x2  = _mm512_add_ps(x0,x1);
                          trm2= _mm512_mul_ps(x2,x2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(pi,h2a),_4a);
                          x1  = _mm512_add_ps(_1, _mm512_div_ps(earg,x0));
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f636_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ph) {

                          using namespace gms::math;
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 h   = _mm512_load_ps(&ph[0]);
                          const __m512 _4pi = _mm512_set1_ps(1.27323954473516268615107010698f);
                          const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0   = _mm512_set1_ps(0.33333333333333333333333333333f);
                          const __m512 c1   = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _4   = _mm512_set1_ps(4.0f);
                          const __m512 _2   = _mm512_set1_ps(2.0f);
                          const __m512 _1   = _mm512_set1_ps(1.0f);
                          register __m512 k04,a2,_4a,a3,h2a,nh2a,pia2,pia3,arg,earg,x0,x1,x2;
                          register __m512 rcs,trm1,trm2,trm3;
                          x0  = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          h2a = _mm512_fmadd_ps(a,_2,h);
                          a3  = _mm512_mul_ps(a2,a);
                          k04 = _mm512_mul_ps(x0,x0);
                          _4a = _mm512_mul_ps(_4,a);
                          trm1= _mm512_mul_ps(_4pi,k04);
                          pia2= _mm512_mul_ps(pi,a2);
                          nh2a= negate_zmm16r4(h2a);
                          x0  = _mm512_mul_ps(c0,_mm512_mul_ps(pia2,h));
                          pia3= _mm512_mul_ps(pi,a3);
                          arg = _mm512_div_ps(nh2a,_4a);
                          x1  = _mm512_mul_ps(c1,pia3);
                          earg= xexpf(arg);
                          x2  = _mm512_add_ps(x0,x1);
                          trm2= _mm512_mul_ps(x2,x2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(pi,h2a),_4a);
                          x1  = _mm512_add_ps(_1, _mm512_div_ps(earg,x0));
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f636_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  ph) {

                          using namespace gms::math;
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 h   = _mm512_loadu_ps(&ph[0]);
                          const __m512 _4pi = _mm512_set1_ps(1.27323954473516268615107010698f);
                          const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0   = _mm512_set1_ps(0.33333333333333333333333333333f);
                          const __m512 c1   = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _4   = _mm512_set1_ps(4.0f);
                          const __m512 _2   = _mm512_set1_ps(2.0f);
                          const __m512 _1   = _mm512_set1_ps(1.0f);
                          register __m512 k04,a2,_4a,a3,h2a,nh2a,pia2,pia3,arg,earg,x0,x1,x2;
                          register __m512 rcs,trm1,trm2,trm3;
                          x0  = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          h2a = _mm512_fmadd_ps(a,_2,h);
                          a3  = _mm512_mul_ps(a2,a);
                          k04 = _mm512_mul_ps(x0,x0);
                          _4a = _mm512_mul_ps(_4,a);
                          trm1= _mm512_mul_ps(_4pi,k04);
                          pia2= _mm512_mul_ps(pi,a2);
                          nh2a= negate_zmm16r4(h2a);
                          x0  = _mm512_mul_ps(c0,_mm512_mul_ps(pia2,h));
                          pia3= _mm512_mul_ps(pi,a3);
                          arg = _mm512_div_ps(nh2a,_4a);
                          x1  = _mm512_mul_ps(c1,pia3);
                          earg= xexpf(arg);
                          x2  = _mm512_add_ps(x0,x1);
                          trm2= _mm512_mul_ps(x2,x2);
                          x0  = _mm512_div_ps(_mm512_mul_ps(pi,h2a),_4a);
                          x1  = _mm512_add_ps(_1, _mm512_div_ps(earg,x0));
                          trm2= _mm512_mul_ps(x1,x1);
                          rcs = _mm512_mul_ps(trm1,_mm512_mul_ps(trm2,trm3));
                          return (rcs);
                }


                   /*
                         Flat-back cone, backscatter RCS.
                         Formula 6.3-9
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f639_zmm16r4(const __m512 gam0,
                                           const __m512 alp,
                                           const __m512 k0,
                                           const __m512 h) {

                          using namespace gms::math;
                          const __m512 c0 = _mm512_set1_ps(-0.25f);
                          const __m512 _0 = _mm512_set1_ps(-0.0f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 gam2,talp,k0h,x0,ear,eai,cer,cei;
                          register __m512 t0r,t0i,trm2r,trm2i,trm1r,trm1i,t1r,t1i;
                          gam2  = _mm512_div_ps(_mm512_mul_ps(gam0,gam0),pi);
                          talp  = xtanf(alp);
                          k0h   = _mm512_mul_ps(k0,h);
                          x0    = _mm512_mul_ps(talp,talp);
                          ear   = _mm512_setzero_ps();
                          t0r   = n1
                          eai   = _mm512_add_ps(k0h,k0h);
                          t0i   = _mm512_sub_ps(eai,_1);
                          cexp_zmm16r4(eai,ear,&cer,&cei);
                          trm1r = ear;
                          trm1i = negate_zmm16r4(x0);
                          cmul_zmm16r4(t0r,t0i,cer,cei,&trm2r,&trm2i);
                          trm2r = _mm512_add_ps(_1,trm2r);
                          trm2i = _mm512_add_ps(_1,trm2i);
                          cmul_zmm16r4(trm1r,trm1i,trm2r,&t0r,&t0i);
                          cabs = cabs_zmm16r4(t0r,t0i);
                          rcs  = _mm512_mul_ps(gam2,cabs);
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f639_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam0, 
                                             const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) palp,
                                             const float * __restrict __ATTR_ALIGN__(64) ph) {

                          using namespace gms::math;
                          register __m512 gam0= _mm512_load_ps(&pgam0[0]);
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 alp = _mm512_load_ps(&palp[0]);
                          register __m512 h   = _mm512_load_ps(&ph[0]);
                          const __m512 c0 = _mm512_set1_ps(-0.25f);
                          const __m512 _0 = _mm512_set1_ps(-0.0f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 gam2,talp,k0h,x0,ear,eai,cer,cei;
                          register __m512 t0r,t0i,trm2r,trm2i,trm1r,trm1i,t1r,t1i;
                          gam2  = _mm512_div_ps(_mm512_mul_ps(gam0,gam0),pi);
                          talp  = xtanf(alp);
                          k0h   = _mm512_mul_ps(k0,h);
                          x0    = _mm512_mul_ps(talp,talp);
                          ear   = _mm512_setzero_ps();
                          t0r   = n1
                          eai   = _mm512_add_ps(k0h,k0h);
                          t0i   = _mm512_sub_ps(eai,_1);
                          cexp_zmm16r4(eai,ear,&cer,&cei);
                          trm1r = ear;
                          trm1i = negate_zmm16r4(x0);
                          cmul_zmm16r4(t0r,t0i,cer,cei,&trm2r,&trm2i);
                          trm2r = _mm512_add_ps(_1,trm2r);
                          trm2i = _mm512_add_ps(_1,trm2i);
                          cmul_zmm16r4(trm1r,trm1i,trm2r,&t0r,&t0i);
                          cabs = cabs_zmm16r4(t0r,t0i);
                          rcs  = _mm512_mul_ps(gam2,cabs);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f639_zmm16r4_u(const float * __restrict  pgam0, 
                                             const float * __restrict  pk0,
                                             const float * __restrict  palp,
                                             const float * __restrict  ph) {

                          using namespace gms::math;
                          register __m512 gam0= _mm512_loadu_ps(&pgam0[0]);
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 alp = _mm512_loadu_ps(&palp[0]);
                          register __m512 h   = _mm512_loadu_ps(&ph[0]);
                          const __m512 c0 = _mm512_set1_ps(-0.25f);
                          const __m512 _0 = _mm512_set1_ps(-0.0f);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          register __m512 gam2,talp,k0h,x0,ear,eai,cer,cei;
                          register __m512 t0r,t0i,trm2r,trm2i,trm1r,trm1i,t1r,t1i;
                          gam2  = _mm512_div_ps(_mm512_mul_ps(gam0,gam0),pi);
                          talp  = xtanf(alp);
                          k0h   = _mm512_mul_ps(k0,h);
                          x0    = _mm512_mul_ps(talp,talp);
                          ear   = _mm512_setzero_ps();
                          t0r   = n1
                          eai   = _mm512_add_ps(k0h,k0h);
                          t0i   = _mm512_sub_ps(eai,_1);
                          cexp_zmm16r4(eai,ear,&cer,&cei);
                          trm1r = ear;
                          trm1i = negate_zmm16r4(x0);
                          cmul_zmm16r4(t0r,t0i,cer,cei,&trm2r,&trm2i);
                          trm2r = _mm512_add_ps(_1,trm2r);
                          trm2i = _mm512_add_ps(_1,trm2i);
                          cmul_zmm16r4(trm1r,trm1i,trm2r,&t0r,&t0i);
                          cabs = cabs_zmm16r4(t0r,t0i);
                          rcs  = _mm512_mul_ps(gam2,cabs);
                          return (rcs);
                 }


                   /*
                        Cone tip scattering RCS.
                        Formula 6.3-10
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6310_zmm16r4(const __m512 gam0,
                                            const __m512 alp) {

                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          register __m512 gam2,tana,x0,tan4,trm1;
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          tana = xtanf(alp);
                          trm1 = _mm512_div_ps(gam2,_16pi);
                          x0   = _mm512_mul_ps(tana,tana);
                          tan4 = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,tan4);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6310_zmm16r4_a(const float * __restrict  __ATTR_ALIGN__(64) pgam0, 
                                              const float * __restrict  __ATTR_ALIGN__(64) palp) {

                          register __m512 gam0= _mm512_load_ps(&pgam0[0]);
                          register __m512 alp = _mm512_load_ps(&palp[0]);
                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          register __m512 gam2,tana,x0,tan4,trm1;
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          tana = xtanf(alp);
                          trm1 = _mm512_div_ps(gam2,_16pi);
                          x0   = _mm512_mul_ps(tana,tana);
                          tan4 = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,tan4);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6310_zmm16r4_u(const float * __restrict  pgam0, 
                                              const float * __restrict  palp) {

                          register __m512 gam0= _mm512_loadu_ps(&pgam0[0]);
                          register __m512 alp = _mm512_loadu_ps(&palp[0]);
                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          register __m512 gam2,tana,x0,tan4,trm1;
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          tana = xtanf(alp);
                          trm1 = _mm512_div_ps(gam2,_16pi);
                          x0   = _mm512_mul_ps(tana,tana);
                          tan4 = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,tan4);
                          return (rcs);
                 }


                   /*
                         Case of flat-back cone joined by the sphere.
                         Axial-incidence backscatter RCS.
                         Formula 6.3-11
                    */

                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6311_zmm16r4(const __m512 gam0,
                                            const __m512 alp,
                                            const __m512 k0,
                                            const __m512 h,
                                            const __m512 z1) {

                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0 = _mm512_set1_ps(0.25f);
                          const __m512 c1 = _mm512_set1_ps(-0.25f);
                          register __m512 gam2,tana,tana2,cosa,seca,k0h,z1h,kz1h,i4r,ir4;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512 t0r,t0i,t1r,t1i,x0,x1,cabs;
                          k0h  = _mm512_mul_ps(k0,h);
                          tana = xtanf(alp);
                          z1h  = _mm512_add_ps(z1,h);
                          tana2= _mm512_mul_ps(tana,tana);
                          kz1h = _mm512_mul_ps(k0,z1h);
                          gam2 = _mm512_div_ps(_mm512_mul_ps(gam0,gam0),pi);
                          i4r  = _mm512_setzero_ps();
                          i4i  = c1;
                          cosa = xcosf(alp);
                          ear1 = i4r;
                          eai1 = _mm512_add_ps(k0h,k0h);
                          seca = _mm512_mul_ps(cosa,cosa);
                          t0r  = i4r;
                          x0   = _mm512_fmsub_ps(h,tana2,z1);
                          t0i  = _mm512_mul_ps(tana2,i4i);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          x1   = _mm512_sub_ps(x0,seca);
                          t1r  = i4r;
                          ear2 = i4r;
                          t1i  = _mm512_sub_ps(t0i,_mm512_mul_ps(c0,x1));
                          eai2 = _mm512_add_ps(kz1h,kz1h);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2 = i4r;
                          cmul_zmm16r4(t1r,t1i,cer1,cei1,&t0r,&t0i);
                          cei2 = _mm512_mul_ps(cei2,c0);
                          x0   = _mm512_sub_ps(t0r,cer2);
                          x1   = _mm512_sub_ps(t0i,cei2);
                          cabs = cabs_zmm16r4(x0,x1);
                          rcs  = _mm512_mul_ps(gam2,cabs);
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6311_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam0, 
                                             const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) palp,
                                             const float * __restrict __ATTR_ALIGN__(64) ph,
                                             const float * __restrict __ATTR_ALIGN__(64) pz1) {

                          register __m512 gam0= _mm512_load_ps(&pgam0[0]);
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 alp = _mm512_load_ps(&palp[0]);
                          register __m512 h   = _mm512_load_ps(&ph[0]);
                          register __m512 z1  = _mm512_load_ps(&pz1[0]);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0 = _mm512_set1_ps(0.25f);
                          const __m512 c1 = _mm512_set1_ps(-0.25f);
                          register __m512 gam2,tana,tana2,cosa,seca,k0h,z1h,kz1h,i4r,ir4;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512 t0r,t0i,t1r,t1i,x0,x1,cabs;
                          k0h  = _mm512_mul_ps(k0,h);
                          tana = xtanf(alp);
                          z1h  = _mm512_add_ps(z1,h);
                          tana2= _mm512_mul_ps(tana,tana);
                          kz1h = _mm512_mul_ps(k0,z1h);
                          gam2 = _mm512_div_ps(_mm512_mul_ps(gam0,gam0),pi);
                          i4r  = _mm512_setzero_ps();
                          i4i  = c1;
                          cosa = xcosf(alp);
                          ear1 = i4r;
                          eai1 = _mm512_add_ps(k0h,k0h);
                          seca = _mm512_mul_ps(cosa,cosa);
                          t0r  = i4r;
                          x0   = _mm512_fmsub_ps(h,tana2,z1);
                          t0i  = _mm512_mul_ps(tana2,i4i);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          x1   = _mm512_sub_ps(x0,seca);
                          t1r  = i4r;
                          ear2 = i4r;
                          t1i  = _mm512_sub_ps(t0i,_mm512_mul_ps(c0,x1));
                          eai2 = _mm512_add_ps(kz1h,kz1h);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2 = i4r;
                          cmul_zmm16r4(t1r,t1i,cer1,cei1,&t0r,&t0i);
                          cei2 = _mm512_mul_ps(cei2,c0);
                          x0   = _mm512_sub_ps(t0r,cer2);
                          x1   = _mm512_sub_ps(t0i,cei2);
                          cabs = cabs_zmm16r4(x0,x1);
                          rcs  = _mm512_mul_ps(gam2,cabs);
                          return (rcs);
                  }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6311_zmm16r4_u(const float * __restrict  pgam0, 
                                             const float * __restrict  pk0,
                                             const float * __restrict  palp,
                                             const float * __restrict  ph,
                                             const float * __restrict  pz1) {

                          register __m512 gam0= _mm512_loadu_ps(&pgam0[0]);
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 alp = _mm512_loadu_ps(&palp[0]);
                          register __m512 h   = _mm512_loadu_ps(&ph[0]);
                          register __m512 z1  = _mm512_loadu_ps(&pz1[0]);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0 = _mm512_set1_ps(0.25f);
                          const __m512 c1 = _mm512_set1_ps(-0.25f);
                          register __m512 gam2,tana,tana2,cosa,seca,k0h,z1h,kz1h,i4r,ir4;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512 t0r,t0i,t1r,t1i,x0,x1,cabs;
                          k0h  = _mm512_mul_ps(k0,h);
                          tana = xtanf(alp);
                          z1h  = _mm512_add_ps(z1,h);
                          tana2= _mm512_mul_ps(tana,tana);
                          kz1h = _mm512_mul_ps(k0,z1h);
                          gam2 = _mm512_div_ps(_mm512_mul_ps(gam0,gam0),pi);
                          i4r  = _mm512_setzero_ps();
                          i4i  = c1;
                          cosa = xcosf(alp);
                          ear1 = i4r;
                          eai1 = _mm512_add_ps(k0h,k0h);
                          seca = _mm512_mul_ps(cosa,cosa);
                          t0r  = i4r;
                          x0   = _mm512_fmsub_ps(h,tana2,z1);
                          t0i  = _mm512_mul_ps(tana2,i4i);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          x1   = _mm512_sub_ps(x0,seca);
                          t1r  = i4r;
                          ear2 = i4r;
                          t1i  = _mm512_sub_ps(t0i,_mm512_mul_ps(c0,x1));
                          eai2 = _mm512_add_ps(kz1h,kz1h);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2 = i4r;
                          cmul_zmm16r4(t1r,t1i,cer1,cei1,&t0r,&t0i);
                          cei2 = _mm512_mul_ps(cei2,c0);
                          x0   = _mm512_sub_ps(t0r,cer2);
                          x1   = _mm512_sub_ps(t0i,cei2);
                          cabs = cabs_zmm16r4(x0,x1);
                          rcs  = _mm512_mul_ps(gam2,cabs);
                          return (rcs);
                  }


                    /*
                          Flat-back cone backscattering RCS.
                          Formula 6.3-14
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6314_zmm16r4(const __m512 a,
                                            const __m512 alp) {

                          const __m512 pi3   = _mm512_set1_ps(31.006276680299820175476315067101f);
                          const __m512 _4pi2 = _mm512_set1_ps(39.478417604357434475337963999505f);
                          const __m512 _3pi  = _mm512_set1_ps(9.424777960769379715387930149839f);
                          const __m512 _3pi2 = _mm512_set1_ps(4.712388980384689857693965074919f)
                          register __m512 rcs,a2,num1,den1,arg,sarg,csc,x0,alp2;
                          a2   = _mm512_mul_ps(a,a);
                          alp2 = _mm512_add_ps(alp,alp);
                          num1 = _mm512_mul_ps(pi3,a2);
                          arg  = _mm512_div_ps(_4pi2,_mm512_add_ps(_3pi,alp2));
                          x0   = _mm512_add_ps(alp,_3pi2);
                          sarg = xsinf(arg);
                          den1 = _mm512_mul_ps(x0,x0);
                          csc  = _mm512_rcp14_ps(sarg);
                          x0   = _mm512_mul_ps(csc,csc);
                          arg  = _mm512_div_ps(num1,den1);
                          rcs  = _mm512_mul_ps(arg,x0);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6314_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) palp) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 alp = _mm512_load_ps(&palp[0]);
                          const __m512 pi3   = _mm512_set1_ps(31.006276680299820175476315067101f);
                          const __m512 _4pi2 = _mm512_set1_ps(39.478417604357434475337963999505f);
                          const __m512 _3pi  = _mm512_set1_ps(9.424777960769379715387930149839f);
                          const __m512 _3pi2 = _mm512_set1_ps(4.712388980384689857693965074919f)
                          register __m512 rcs,a2,num1,den1,arg,sarg,csc,x0,alp2;
                          a2   = _mm512_mul_ps(a,a);
                          alp2 = _mm512_add_ps(alp,alp);
                          num1 = _mm512_mul_ps(pi3,a2);
                          arg  = _mm512_div_ps(_4pi2,_mm512_add_ps(_3pi,alp2));
                          x0   = _mm512_add_ps(alp,_3pi2);
                          sarg = xsinf(arg);
                          den1 = _mm512_mul_ps(x0,x0);
                          csc  = _mm512_rcp14_ps(sarg);
                          x0   = _mm512_mul_ps(csc,csc);
                          arg  = _mm512_div_ps(num1,den1);
                          rcs  = _mm512_mul_ps(arg,x0);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6314_zmm16r4_u( const float * __restrict  pa,
                                             const float * __restrict  palp) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 alp = _mm512_loadu_ps(&palp[0]);
                          const __m512 pi3   = _mm512_set1_ps(31.006276680299820175476315067101f);
                          const __m512 _4pi2 = _mm512_set1_ps(39.478417604357434475337963999505f);
                          const __m512 _3pi  = _mm512_set1_ps(9.424777960769379715387930149839f);
                          const __m512 _3pi2 = _mm512_set1_ps(4.712388980384689857693965074919f)
                          register __m512 rcs,a2,num1,den1,arg,sarg,csc,x0,alp2;
                          a2   = _mm512_mul_ps(a,a);
                          alp2 = _mm512_add_ps(alp,alp);
                          num1 = _mm512_mul_ps(pi3,a2);
                          arg  = _mm512_div_ps(_4pi2,_mm512_add_ps(_3pi,alp2));
                          x0   = _mm512_add_ps(alp,_3pi2);
                          sarg = xsinf(arg);
                          den1 = _mm512_mul_ps(x0,x0);
                          csc  = _mm512_rcp14_ps(sarg);
                          x0   = _mm512_mul_ps(csc,csc);
                          arg  = _mm512_div_ps(num1,den1);
                          rcs  = _mm512_mul_ps(arg,x0);
                          return (rcs);
                }


                   /*
                          Scattering from a thick cylinder of length
                          being a slant height of the cone of radius
                          (4/9*a*sec(alpha)), i.e. RCS(PI/2-alpha).
                          Formula 6.3-18
                      */


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6318_zmm16r4(const __m512 gam0,
                                            const __m512 k0h,
                                            const __m512 alp) {

                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0 = _mm512_set1_ps(0.33333333333333333333333333333f);
                          register __m512 rcs,gam2,k0h3,x0,salp,calp,seca,sqr,x1,sabs;
                          x0   = _mm512_mul_ps(gam0,gam0);
                          k0h3 = _mm512_mul_ps(k0h,c0);
                          gam2 = _mm512_div_ps(x0,pi);
                          salp = xsinf(alp);
                          calp = xcosf(alp);
                          x1   = _mm512_mul_ps(k0h,salp);
                          x0   = _mm512_div_ps(x1,pi);
                          seca = _mm512_rcp14_ps(calp);
                          sqr  = _mm512_sqrt_ps(x0);
                          x1   = _mm512_mul_ps(seca,seca);
                          x0   = _mm512_mul_ps(k0h3,_mm512_mul_ps(sqr,x1));
                          sabs = _mm512_abs_ps(x0);
                          x1   = _mm512_mul_ps(sabs,sabs);
                          rcs  = _mm512_mul_ps(gam2,x1);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6318_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam0, 
                                             const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const float * __restrict __ATTR_ALIGN__(64) palp) {

                          register __m512 gam0 = _mm512_load_ps(&pgam0[0]);
                          register __m512 k0h  = _mm512_load_ps(&pk0[0]);
                          register __m512 alp  = _mm512_load_ps(&palp[0]);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0 = _mm512_set1_ps(0.33333333333333333333333333333f);
                          register __m512 rcs,gam2,k0h3,x0,salp,calp,seca,sqr,x1,sabs;
                          x0   = _mm512_mul_ps(gam0,gam0);
                          k0h3 = _mm512_mul_ps(k0h,c0);
                          gam2 = _mm512_div_ps(x0,pi);
                          salp = xsinf(alp);
                          calp = xcosf(alp);
                          x1   = _mm512_mul_ps(k0h,salp);
                          x0   = _mm512_div_ps(x1,pi);
                          seca = _mm512_rcp14_ps(calp);
                          sqr  = _mm512_sqrt_ps(x0);
                          x1   = _mm512_mul_ps(seca,seca);
                          x0   = _mm512_mul_ps(k0h3,_mm512_mul_ps(sqr,x1));
                          sabs = _mm512_abs_ps(x0);
                          x1   = _mm512_mul_ps(sabs,sabs);
                          rcs  = _mm512_mul_ps(gam2,x1);
                          return (rcs);
                 }
                                            



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f6318_zmm16r4_u(const float * __restrict  pgam0, 
                                             const float * __restrict  pk0h,
                                             const float * __restrict  palp) {

                          register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);
                          register __m512 k0h  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 alp  = _mm512_loadu_ps(&palp[0]);
                          const __m512 pi = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 c0 = _mm512_set1_ps(0.33333333333333333333333333333f);
                          register __m512 rcs,gam2,k0h3,x0,salp,calp,seca,sqr,x1,sabs;
                          x0   = _mm512_mul_ps(gam0,gam0);
                          k0h3 = _mm512_mul_ps(k0h,c0);
                          gam2 = _mm512_div_ps(x0,pi);
                          salp = xsinf(alp);
                          calp = xcosf(alp);
                          x1   = _mm512_mul_ps(k0h,salp);
                          x0   = _mm512_div_ps(x1,pi);
                          seca = _mm512_rcp14_ps(calp);
                          sqr  = _mm512_sqrt_ps(x0);
                          x1   = _mm512_mul_ps(seca,seca);
                          x0   = _mm512_mul_ps(k0h3,_mm512_mul_ps(sqr,x1));
                          sabs = _mm512_abs_ps(x0);
                          x1   = _mm512_mul_ps(sabs,sabs);
                          rcs  = _mm512_mul_ps(gam2,x1);
                          return (rcs);
                 }




                  


          }// radiolocation


}// gms

















#endif  /*__GMS_RCS_CONE_WEDGE_ZMM16R4_HPP__*/
