



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

/*
      Adapted from the DISORT "BDREF.f" file.
*/





#include "GMS_bdref_avx512.h"
#if (USE_SLEEF_LIB) == 1
#include "GMS_sleefsimddp.hpp"
#include "GMS_sleefsimdsp.hpp"
#endif


			__m512d
			gms::math::bdrf_hapke_zmm8r8(const __m512d mup,
			                  const __m512d mu,
					  const __m512d dphi,
					  const __m512d b0,
					  const __m512d hh,
					  const __m512d w,
					  const __m512d pi) {

                           const __m512d _1   = _mm512_set1_pd(1.0);
			   const __m512d _1_2 = _mm512_set1_pd(0.5);
			   const __m512d _4   = _mm512_set1_pd(4.0);
			   __m512d t0,t1,t2,t3,c0,c1,cdphi,brdf;
			   __m512d calpha,alpha,p,b,h0,gamma,h;
#if (USE_SLEEF_LIB) == 1
                           cdphi = xcos(dphi);	
#else		   
			   cdphi = _mm512_cos_pd(dphi);
#endif
			   c0    = _mm512_fmadd_pd(_2,mup,_1);
			   t0 = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                       _mm512_mul_pd(mu,mu)));
			   c1    = _mm512_fmadd_pd(_2,mu,_1); 
			   t1 = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                       _mm512_mul_pd(mup,mup)));
			   calpha = _mm512_mul_pd(_mm512_fmsub_pd(mu,mup,
			                                      _mm512_mul_pd(t0,t1)),cdphi);
			   gamma  = _mm512_sqrt_pd(_mm512_sub_pd(_1,w);
			   alpha  = _mm512_acos_pd(calpha);
			   p      = _mm512_fmadd_pd(_1_2,calpha,_1);
			   t2     = _mm512_add_pd(hh,_mm512_tan_pd(_mm512_mul_pd(alpha,_1_2)));
			   t3     = _mm512_mul_pd(b0,hh);
			   b      = _mm512_div_pd(t3,t2);
			   h0     = _mm512_div_pd(c0,_mm512_mul_pd(c0,gamma));
			   h      = _mm512_div_pd(c1,_mm512_mul_pd(c1,gamma));
			   t0     = _mm512_div_pd(w,_mm512_mul_pd(_4,pi));
			   t1     = _mm512_add_pd(mu,mup);
			   t2     = _mm512_fmadd_pd(_mm512_add_pd(_1,b),p,
			                                      _mm512_fmsub_pd(h0,h,_1));
			   brdf   = _mm512_mul_pd(_mm512_div_pd(t0,t1),t2);
			   return (brdf);
		      }


		     
			__m512
			gms::math::bdrf_hapke_zmm16r4(const __m512 mup,
			                   const __m512 mu,
					   const __m512 dphi,
					   const __m512 b0,
					   const __m512 hh,
					   const __m512 w,
					   const __m512 pi) {

                           const __m512 _1   = _mm512_set1_ps(1.0f);
			   const __m512 _1_2 = _mm512_set1_ps(0.5f);
			   const __m512 _4   = _mm512_set1_ps(4.0f);
			   __m512 t0,t1,t2,t3,c0,c1,cdphi,brdf;
			   __m512 calpha,alpha,p,b,h0,gamma,h;
			   cdphi = _mm512_cos_ps(dphi);
			   c0    = _mm512_fmadd_ps(_2,mup,_1);
			   t0 = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                       _mm512_mul_ps(mu,mu)));
			   c1    = _mm512_fmadd_ps(_2,mu,_1); 
			   t1 = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                       _mm512_mul_ps(mup,mup)));
			   calpha = _mm512_mul_ps(_mm512_fmsub_ps(mu,mup,
			                                      _mm512_mul_ps(t0,t1)),cdphi);
			   gamma  = _mm512_sqrt_ps(_mm512_sub_ps(_1,w));
			   alpha  = _mm512_acos_ps(calpha);
			   p      = _mm512_fmadd_ps(_1_2,calpha,_1);
			   t2     = _mm512_add_ps(hh,_mm512_tan_ps(_mm512_mul_ps(alpha,_1_2)));
			   t3     = _mm512_mul_ps(b0,hh);
			   b      = _mm512_div_ps(t3,t2);
			   h0     = _mm512_div_ps(c0,_mm512_mul_ps(c0,gamma));
			   h      = _mm512_div_ps(c1,_mm512_mul_ps(c1,gamma));
			   t0     = _mm512_div_ps(w,_mm512_mul_ps(_4,pi));
			   t1     = _mm512_add_ps(mu,mup);
			   t2     = _mm512_fmadd_ps(_mm512_add_ps(_1,b),p,
			                                      _mm512_fmsub_ps(h0,h,_1));
			   brdf   = _mm512_mul_pd(_mm512_div_pd(t0,t1),t2);
			   return (brdf);
		      }


                       
			__m512d
			gms::math::bdrf_rpv_zmm8r8(const __m512d mu_i,
			                const __m512d mu_r,
					const __m512d dphi,
					const __m512d rh00,
					const __m512d kappa,
					const __m512d g_hg,
					const __m512d h0) {

                           const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
			   const __m512d _1   = _mm512_set1_pd(1.0);
			   const __m512d _2   = _mm512_set1_pd(2.0);
			   const __m512d _1_5 = _mm512_set1_pd(1.5);
			   const __m512d g_hg2= _mm512_mul_pd(g_hg,g_hg);
			   __m512d brdf,cos_alpha,sin_i,sin_r,tan_i,tan_r;
			   __m512d g_sq,g,f,cdphi;
			   __m512d t0,t1,t2,t3;
			   cdphi = _mm512_cos_pd(dphi);
			   t0    = _mm512_sub_pd(_1,_mm512_mul_pd(mu_i,mu_i));
			   sin_i = _mm512_sqrt_pd(t0);
			   tan_i = sin_i/mu_i;
			   t1    = _mm512_sub_pd(_1,_mm512_mul_pd(mu_r,mu_r));
			   sin_r = _mm512_sqrt_pd(t1);
			   tan_r = sin_r/mu_r;
			   cos_alpha = _mm512_fmsub_pd(mu_i,mu_r,
			                         _mm512_mul_pd(cdphi,
						           _mm512_mul_pd(sin_i,sin_r)));
			   t2    = _mm512_fmadd_pd(tan_i,tan_i,
			                                   _mm512_mul_pd(tan_r,tan_r));
			   t3    = _mm512_mul_pd(_mm512_mul_pd(_2,tan_i),
			                         _mm512_mul_pd(tan_r,cdphi));
			   g_sq  = _mm512_add_pd(t2,t3);
			   g     = _mm512_sqrt_pd(g_sq);
			   t0    = _mm512_sub_pd(_1,g_hg2);
			   t1    = _mm512_add_pd(_mm512_add_pd(_1,g_hg2),
			                          _mm512_mul_pd(_1,
						            _mm512_mul_pd(g_hg,cos_alpha)));
			   f     = _mm512_pow_pd(_mm512_div_pd(t0,t1),_1_5);
			   t2    = _mm512_mul_pd(_mm512_mul_pd(mu_i,mu_r),
			                         _mm512_add_pd(mu_i,mu_r));
			   t3    = _mm512_add_pd(_1,
			                     _mm512_div_pd(_mm512_sub_pd(_1,h0),
					                   _mm512_add_pd(_1,g)));
			   t0    = _mm512_pow_pd(t2,_mm512_sub_pd(kappa,_1));
			   brdf  = _mm512_mul_pd(_mm512_mul_pd(rh00,t0),
			                         _mm512_mul_pd(f,t3));
			   return (brdf);
		     }


		       
			__m512
			gms::math::bdrf_rpv_zmm16r4(const __m512 mu_i,
			                const __m512 mu_r,
					const __m512 dphi,
					const __m512 rh00,
					const __m512 kappa,
					const __m512 g_hg,
					const __m512 h0) {

                           const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
			   const __m512 _1   = _mm512_set1_ps(1.0);
			   const __m512 _2   = _mm512_set1_ps(2.0);
			   const __m512 _1_5 = _mm512_set1_ps(1.5);
			   const __m512 g_hg2= _mm512_mul_ps(g_hg,g_hg);
			   __m512 brdf,cos_alpha,sin_i,sin_r,tan_i,tan_r;
			   __m512 g_sq,g,f,cdphi;
			   __m512 t0,t1,t2,t3;
			   cdphi = _mm512_cos_ps(dphi);
			   t0    = _mm512_sub_ps(_1,_mm512_mul_ps(mu_i,mu_i));
			   sin_i = _mm512_sqrt_ps(t0);
			   tan_i = sin_i/mu_i;
			   t1    = _mm512_sub_ps(_1,_mm512_mul_ps(mu_r,mu_r));
			   sin_r = _mm512_sqrt_ps(t1);
			   tan_r = sin_r/mu_r;
			   cos_alpha = _mm512_fmsub_ps(mu_i,mu_r,
			                         _mm512_mul_ps(cdphi,
						           _mm512_mul_ps(sin_i,sin_r)));
			   t2    = _mm512_fmadd_ps(tan_i,tan_i,
			                                   _mm512_mul_ps(tan_r,tan_r));
			   t3    = _mm512_mul_ps(_mm512_mul_ps(_2,tan_i),
			                         _mm512_mul_ps(tan_r,cdphi));
			   g_sq  = _mm512_add_ps(t2,t3);
			   g     = _mm512_sqrt_ps(g_sq);
			   t0    = _mm512_sub_ps(_1,g_hg2);
			   t1    = _mm512_add_ps(_mm512_add_ps(_1,g_hg2),
			                          _mm512_mul_ps(_1,
						            _mm512_mul_ps(g_hg,cos_alpha)));
			   f     = _mm512_pow_ps(_mm512_div_ps(t0,t1),_1_5);
			   t2    = _mm512_mul_ps(_mm512_mul_ps(mu_i,mu_r),
			                         _mm512_add_ps(mu_i,mu_r));
			   t3    = _mm512_add_ps(_1,
			                     _mm512_div_ps(_mm512_sub_ps(_1,h0),
					                   _mm512_add_ps(_1,g)));
			   t0    = _mm512_pow_ps(t2,_mm512_sub_ps(kappa,_1));
			   brdf  = _mm512_mul_ps(_mm512_mul_ps(rh00,t0),
			                         _mm512_mul_ps(f,t3));
			   return (brdf);
		     }


		       
			__m512d
			gms::math::bdrf_rossli_zmm8r8(const __m512d mu_i,
			                   const __m512d mu_r,
					   const __m512d dphi,
					   const __m512d k_iso,
					   const __m512d k_vol,
					   const __m512d k_geo,
					   const __m512d alpha0) {

                           /*
                                
c +--------------------------------------------------------------------
c Version 3: Ross-Li BRDF
c   Input:
c
c   MU_I:    absolute cosine of incident polar angle (positive)
c   MU_R:    absolute cosine of reflected polar angle (positive)
c   DPHI:  relative azimuth to incident vector; (pi - dphi), sun-view relative
c          azimuth sun located at phi = 180, while incident solar beam located
c          at phi = 0
c   K_ISO:   BRDF parameter, isotropic scattering kernel
c   K_VOL:   BRDF parameter, volume scattering kernel
c   K_GEO:   BRDF parameter, geometry scattering kernel
c   ALPHA0:  BRDF parameter, control hot spot (back scattering direction)
c
c   Output:
c   BRDF:  Ross-Li BRDF
c
c +--------------------------------------------------------------------
                            */

			    const __m512d sdphi= _mm512_sin_pd(dphi);
			    const __m512d _0   = _mm512_setzero_pd();
			    const __m512d _1   = _mm512_set1_pd(1.0);
			    const __m512d n1   = _mm512_set1_pd(-1.0);
			    const __m512d _2   = _mm512_set1_pd(2.0);
			    const __m512d _4   = _mm512_set1_pd(4.0);
			    const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
			    const __m512d _3pi = _mm512_set1_pd(3.0*3.14159265358979323846264338328);
			    const __m512d pih  = _mm512_set1_pd(0.5*3.14159265358979323846264338328);
			    const __m512d _43pi= _mm512_set1_pd(4.188790204786390984616857844373);
			    const __m512d _0_3 = _mm512_set1_pd(0.333333333333333333333333333333333);
			    const __m512d cdphi= _mm512_cos_pd(dphi);
			    __m512d brdf,f_vol,f_geo;
			    __m512d cos_alpha,sin_alpha,cos_alpha1,alpha;
			    __m512d sin_i,sin_r,tan_i,tan_r;
			    __m512d sin_i1,sin_r1,cos_i1,cos_r1,tan_i1,tan_r1;
			    __m512d g_sq,cos_t,t,st,ct;
			    __m512d c,t0,t1,t2,t3,c0,c1;
                            __m512d ratio_hb,ratio_br;
 
			    sin_i    = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                                       _mm512_mul_pd(mu_i,mu_i)));
			    tan_i    = _mm512_div_pd(sin_i,mu_i)
			    tan_i1   = _mm512_mul_pd(ratio_br,tan_i);
			    ratio_hb = _1;
			    sin_r    = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                                       _mm512_mul_pd(mu_r,mu_r)));
			    tan_r    = _mm512_div_pd(sin_r,mu_r);
			    ratio_br = _2;
			    tan_r1   = _mm512_mul_pd(ratio_br,tan_r);
			    c0       = _mm512_sqrt_pd(_mm512_fmadd_pd(tan_i1,tan_i1,_1));
			    cos_i1   = _mm512_div_pd(_1,c0);
			    cos_alpha= _mm512_fmsub_pd(mu_i,mu_r,
			                                   _mm512_mul_pd(cdphi,
							             _mm512_mul_pd(sin_i,sin_r)));
			    c1       = _mm512_sqrt_pd(_mm512_fmadd_pd(tan_r1,tan_r1,_1));
			    cos_r1   = _mm512_div_pd(_1,c1);
			    sin_alpha= _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                                       _mm512_mul_pd(cos_alpha,cos_alpha)));
			    alpha    = _mm512_acos_pd(cos_alpha);
			    t0       = _mm512_add_pd(_1,_mm512_div_pd(alpha,alpha0));
			    c        = _mm512_add_pd(_1,_mm512_div_pd(_1,t0));
			    sin_i1   = _mm512_div_pd(tan_i1,c0);
			    t1       = _mm512_div_pd(_1,_mm512_add_pd(mu_i,mu_r));
			    t2       = _mm512_mul_pd(c,_mm512_fmsub_pd(pih,alpha,
			                                       _mm512_mul_pd(cos_alpha,sin_alpha)));
			    f_vol    = _mm512_fmsub_pd(_43pi,_mm512_mul_pd(t1,t2),_0_3);
			    sin_r1   = _mm512_div_pd(tan_r1,c1);
			    cos_alpha1= _mm512_fmsub_pd(cosi1,cos_r1,_mm512_mul_pd(
			                                            _mm512_mul_pd(sin_i1,sin_r1),cdphi));
			    t0       = _mm512_fmadd_pd(tan_i1,tan_i1,_mm512_mul_pd(tan_r1,tan_r1));
			    t1       = _mm512_mul_pd(_mm512_mul_pd(_2,tan_i1),
			                             _mm512_mul_pd(tan_r1,cdphi));
			    g_sq     = _mm512_add_pd(t0,t1);
			    t2       = _mm512_mul_pd(ratio_hb,
			                             _mm512_div_pd(_mm512_mul_pd(cos_i1,cos_r1),
			                                    _mm512_add_pd(cos_i1,cos_r1)));
			    t3       = _mm512_mul_pd(tan_t1,_mm512_mul_pd(tan_r1,sdphi));
			    t3       = _mm512_mul_pd(t3,t3);
			    cos_t    = _mm512_mul_pd(t2,_mm512_sqrt_pd(_mm512_mul_pd(g_sq,t3)));
			    if(_mm512_cmp_pd_mask(cos_t,_1,_CMP_LE_OQ) &&
			       _mm512_cmp_pd_mask(cos_t,n1,_CMP_GE_OQ)) {
                               t = _mm512_acos_pd(cos_t);
			    }
			    else {
                               t = _0;
			    }
			    t0 = _mm512_div_pd(_mm512_add_pd(cos_i1,cos_r1),
			                       _mm512_mul_pd(pi,_mm512_mul_pd(cos_i1,cos_r1)));
			    st = _mm512_sin_pd(t);
			    ct = _mm512_cos_pd(t);
			    t1 = _mm512_mul_pd(_mm512_sub_pd(t,st),
			                       _mm512_sub_pd(ct,pi));
			    t2 = _mm512_div_pd(_mm512_add_pd(_1,cos_alpha1),
			                       _mm512_mul_pd(_2,_mm512_mul_pd(cos_i1,cos_r1)));
			    f_geo = _mm512_add_pd(_mm512_mul_pd(t0,t1),t2);
			    brdf  = _mm512_fmadd_pd(k_iso,k_geo,_mm512_fmadd_pd(f_geo,k_vol,f_vol));
			    return (brdf);
		    }


		       
			__m512
			gms::math::bdrf_rossli_zmm16r4(const __m512 mu_i,
			                   const __m512 mu_r,
					   const __m512 dphi,
					   const __m512 k_iso,
					   const __m512 k_vol,
					   const __m512 k_geo,
					   const __m512 alpha0) {

                           /*
                                
c +--------------------------------------------------------------------
c Version 3: Ross-Li BRDF
c   Input:
c
c   MU_I:    absolute cosine of incident polar angle (positive)
c   MU_R:    absolute cosine of reflected polar angle (positive)
c   DPHI:  relative azimuth to incident vector; (pi - dphi), sun-view relative
c          azimuth sun located at phi = 180, while incident solar beam located
c          at phi = 0
c   K_ISO:   BRDF parameter, isotropic scattering kernel
c   K_VOL:   BRDF parameter, volume scattering kernel
c   K_GEO:   BRDF parameter, geometry scattering kernel
c   ALPHA0:  BRDF parameter, control hot spot (back scattering direction)
c
c   Output:
c   BRDF:  Ross-Li BRDF
c
c +--------------------------------------------------------------------
                            */

			    const __m512 sdphi= _mm512_sin_ps(dphi);
			    const __m512 _0   = _mm512_setzero_ps();
			    const __m512 _1   = _mm512_set1_ps(1.0f);
			    const __m512 n1   = _mm512_set1_ps(-1.0f);
			    const __m512 _2   = _mm512_set1_ps(2.0f);
			    const __m512 _4   = _mm512_set1_ps(4.0f);
			    const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
			    const __m512 _3pi = _mm512_set1_ps(3.0f*3.14159265358979323846264338328f);
			    const __m512 pih  = _mm512_set1_ps(0.5f*3.14159265358979323846264338328f);
			    const __m512 _43pi= _mm512_set1_ps(4.188790204786390984616857844373f);
			    const __m512 _0_3 = _mm512_set1_ps(0.333333333333333333333333333333333f);
			    const __m512 cdphi= _mm512_cos_ps(dphi);
			    __m512 brdf,f_vol,f_geo;
			    __m512 cos_alpha,sin_alpha,cos_alpha1,alpha;
			    __m512 sin_i,sin_r,tan_i,tan_r;
			    __m512 sin_i1,sin_r1,cos_i1,cos_r1,tan_i1,tan_r1;
			    __m512 g_sq,cos_t,t,st,ct;
			    __m512 c,t0,t1,t2,t3,c0,c1;
                            __m512 ratio_hb,ratio_br;
 
			    sin_i    = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                                       _mm512_mul_ps(mu_i,mu_i)));
			    tan_i    = _mm512_div_ps(sin_i,mu_i)
			    tan_i1   = _mm512_mul_ps(ratio_br,tan_i);
			    ratio_hb = _1;
			    sin_r    = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                                       _mm512_mul_ps(mu_r,mu_r)));
			    tan_r    = _mm512_div_ps(sin_r,mu_r);
			    ratio_br = _2;
			    tan_r1   = _mm512_mul_ps(ratio_br,tan_r);
			    c0       = _mm512_sqrt_ps(_mm512_fmadd_ps(tan_i1,tan_i1,_1));
			    cos_i1   = _mm512_div_ps(_1,c0);
			    cos_alpha= _mm512_fmsub_ps(mu_i,mu_r,
			                                   _mm512_mul_ps(cdphi,
							             _mm512_mul_ps(sin_i,sin_r)));
			    c1       = _mm512_sqrt_ps(_mm512_fmadd_ps(tan_r1,tan_r1,_1));
			    cos_r1   = _mm512_div_ps(_1,c1);
			    sin_alpha= _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                                       _mm512_mul_ps(cos_alpha,cos_alpha)));
			    alpha    = _mm512_acos_ps(cos_alpha);
			    t0       = _mm512_add_ps(_1,_mm512_div_ps(alpha,alpha0));
			    c        = _mm512_add_ps(_1,_mm512_div_ps(_1,t0));
			    sin_i1   = _mm512_div_ps(tan_i1,c0);
			    t1       = _mm512_div_ps(_1,_mm512_add_ps(mu_i,mu_r));
			    t2       = _mm512_mul_ps(c,_mm512_fmsub_ps(pih,alpha,
			                                       _mm512_mul_ps(cos_alpha,sin_alpha)));
			    f_vol    = _mm512_fmsub_ps(_43pi,_mm512_mul_ps(t1,t2),_0_3);
			    sin_r1   = _mm512_div_ps(tan_r1,c1);
			    cos_alpha1= _mm512_fmsub_ps(cosi1,cos_r1,_mm512_mul_ps(
			                                            _mm512_mul_ps(sin_i1,sin_r1),cdphi));
			    t0       = _mm512_fmadd_ps(tan_i1,tan_i1,_mm512_mul_ps(tan_r1,tan_r1));
			    t1       = _mm512_mul_ps(_mm512_mul_ps(_2,tan_i1),
			                             _mm512_mul_ps(tan_r1,cdphi));
			    g_sq     = _mm512_add_ps(t0,t1);
			    t2       = _mm512_mul_ps(ratio_hb,
			                             _mm512_div_ps(_mm512_mul_ps(cos_i1,cos_r1),
			                                    _mm512_add_ps(cos_i1,cos_r1)));
			    t3       = _mm512_mul_ps(tan_t1,_mm512_mul_ps(tan_r1,sdphi));
			    t3       = _mm512_mul_ps(t3,t3);
			    cos_t    = _mm512_mul_ps(t2,_mm512_sqrt_ps(_mm512_mul_ps(g_sq,t3)));
			    if(_mm512_cmp_ps_mask(cos_t,_1,_CMP_LE_OQ) &&
			       _mm512_cmp_ps_mask(cos_t,n1,_CMP_GE_OQ)) {
                               t = _mm512_acos_ps(cos_t);
			    }
			    else {
                               t = _0;
			    }
			    t0 = _mm512_div_ps(_mm512_add_ps(cos_i1,cos_r1),
			                       _mm512_mul_ps(pi,_mm512_mul_ps(cos_i1,cos_r1)));
			    st = _mm512_sin_ps(t);
			    ct = _mm512_cos_ps(t);
			    t1 = _mm512_mul_ps(_mm512_sub_ps(t,st),
			                       _mm512_sub_ps(ct,pi));
			    t2 = _mm512_div_ps(_mm512_add_ps(_1,cos_alpha1),
			                       _mm512_mul_ps(_2,_mm512_mul_ps(cos_i1,cos_r1)));
			    f_geo = _mm512_add_ps(_mm512_mul_ps(t0,t1),t2);
			    brdf  = _mm512_fmadd_ps(k_iso,k_geo,_mm512_fmadd_ps(f_geo,k_vol,f_vol));
			    return (brdf);
		    }

#include "GMS_simd_utils.hpp"

		      
			__m512d
			gms::math::brdf_ocean_zmm8r8(const bool do_shadow,
			                  const __m512d refrac_idx,
					  const __m512d ws,
					  const __m512d mu_i,
					  const __m512d mu_r,
					  const __m512d dphi) {
/*
c +--------------------------------------------------------------------
c Version 3: 1D Gaussian Rough Ocean BRDF
c   Input:
c
c   mu_i:         absolute cosine of incident polar angle (positive)
c   mu_r:         absolute cosine of reflected polar angle (positive)
c   dphi:         relative azimuth (radians) 
c   do_shadow:    BRDF parameter, open/close shadow effect 
c   refrac_index: BRDF parameter, refractive index of boundary media (water)
c   ws:           BRDF parameter, wind speed (m/s)
c
c   Output:
c
c   brdf:         1D Gaussian Rough Ocean BRDF
c          
c +--------------------------------------------------------------------
*/
                             const __m512d pi     = _mm512_set1_pd(3.14159265358979323846264338328);
                             const __m512d _1     = _mm512_set1_pd(1.0);
			     const __m512d _2     = _mm512_set1_pd(2.0);
			     const __m512d _4     = _mm512_set1_pd(4.0);
			     const __m512d _1_2   = _mm512_set1_pd(0.5);
			     const __m512d c1     = _mm512_set1_pd(0.003);
			     const __m512d c2     = _mm512_set1_pd(0.00512);
			     const __m512d muimur = _mm512_add_pd(mu_i,mu_r);
			     const __m512d cdphi  = _mm512_cos_pd(dphi);
			     __m512d brdf,shadow;
			     __m512d t0,t1,t2,t3,sa1,sa2,oct;
			     __m512d sin_i,sin_r,cos_theta,sigma_sq;
			     __m512d mu_n_sq,p,n_i,n_t,cos_li,cos_lt;
			     __m512d sin_li,sin_lt,r_s,r_p,r;
                             sin_i    = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                                     _mm512_mul_pd(mu_i,mu_i)));
			     n_i      = _1;
			     sigma_sq = _mm512_fmadd_pd(c2,ws,c1);
			     n_t      = refrac_idx;
			     sun_r    = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                                     _mm512_mul_pd(mu_r,mu_r)));
			     cos_theta= _mm512_fmadd_pd(zmm8r8_negate(mu_i),mu_r,
			                                         _mm512_mul_pd(cdphi,
								          _mm512_mul_pd(sin_i,sin_r)));
			     t0       = _mm512_mul_pd(_2,_mm512_sub_pd(_1,cos_theta));
			     mu_n_sq  = _mm512_div_pd(_mm512_mul_pd(muimur,muimur),t0);
			     oct      = _mm512_sub_pd(_1,cos_theta); 
			     t1       = _mm512_div_pd(_1,_mm512_mul_pd(pi,sigma_sq));
			     sa1      = _mm512_mul_pd(_1_2,oct);
			     t2       = _mm512_div_pd(_mm512_sub_pd(_1,mu_n_sq),
			                              _mm512_mul_pd(sigma_sq,mu_n_sq));
			     sa2      = _mm512_sub_pd(_1,sa1);
			     p        = _mm512_mul_pd(t1,_mm512_exp_pd(negate_zmm8r8(t2)));
			     sin_li   = _mm512_sqrt_pd(sa2);
			     cos_li   = _mm512_sqrt_pd(sa1);
			     sin_lt   = _mm512_mul_pd(n_i,_mm512_div_pd(sin_li,n_t));
			     cos_lt   = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                                         _mm512_mul_pd(sin_lt,sin_lt)));
			     t0       = _mm512_mul_pd(n_i,cos_lt);
			     t1       = _mm512_mul_pd(n_t,cos_li);
			     t2       = _mm512_mul_pd(n_t,cos_lt);
			     t3       = _mm512_mul_pd(n_i,cos_li);
			     r_s      = _mm512_div_pd(_mm512_sub_pd(t3,t2),
			                              _mm512_add_pd(t3,t2));
			     r_p      = _mm512_div_pd(_mm512_sub_pd(t1,t0),
			                              _mm512_add_pd(t0,t1));
			     r        = _mm512_mul_pd(_1_2,_mm512_fmadd_pd(r_s,r_s,
			                              _mm512_mul_pd(r_p,r_p)));
			     t0       = _mm512_mul_pd(p,r);
			     t1       = _mm512_mul_pd(_mm512_mul_pd(mu_i,mu_r),
			                              _mm512_mul_pd(mu_n_sq,mu_n_sq));
			     brdf     = _mm512_div_pd(t0,_mm512_mul_pd(_4,t1));
			     if(do_shadow) {
                                t2    = _mm512_div_pd(_1,shadow_eta_zmm8r8(mu_i,sigma_sq,pi));
				t3    = _mm512_add_pd(shadow_eta_zmm8r8(mu_r,sigma_sq,pi),_1);
				shadow= _mm512_add_pd(t2,t3);
				brdf  = _mm512_add_pd(brdf,shadow);
				return (brdf);
			     }
			     return (brdf);                     
		      }


		     
			__m512d
			gms::math::shadow_eta_zmm8r8(const __m512d cos_theta,
			                  const __m512d sigma_sq,
					  const __m512d pi) {

/*
c +--------------------------------------------------------------------
c Version 3: shadow effect function
c            called by OCEABRDF2
c   Input:
c
c   COS_THETA     absolute cosine of incident/reflected polar angle (positive)
c   SIGMA_SQ      slope variance 
c   PI            3.141592653... constant
c
c   Output:
c
c   SHADOW_ETA:   shadow function
c +--------------------------------------------------------------------
*/

                              const __m512d  _1  = _mm512_set1_pd(1.0);
			      const __m512d _1_2 = _mm512_set1_pd(0.5);
			      const __m512d ssqpi= _mm512_div_pd(sigma_sq,pi);
			      __m512d seta;
			      __m512d sin_theta,mu,sigma_sq,term1,term2;
			      __m512d t0,t1;
			      sin_theta = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                                          _mm512_mul_pd(cos_theta,cos_theta)));
			      mu        = _mm512_div_pd(cos_theta,sin_theta);
			      t0        = _mm512_div_pd(_mm512_sqrt_pd(ssqpi),mu);
			      t1        = _mm512_mul_pd(zmm8r8_negate(mu),
			                                        _mm512_div_pd(mu,sigma_sq));
			      term1     = _mm512_mul_pd(t0,_mm512_exp_pd(t1));
			      term2     = _mm512_erfc_pd(_mm512_div_pd(mu,
			                                        _mm512_sqrt_pd(sigma_sq)));
			      seta      = _mm512_mul_pd(_1_2,
			                              _mm512_sub_pd(term1,term2));
			      return (seta);
		       }
					 


		      
			__m512
			gms::math::brdf_ocean_zmm16r4(const bool do_shadow,
			                  const __m512 refrac_idx,
					  const __m512 ws,
					  const __m512 mu_i,
					  const __m512 mu_r,
					  const __m512 dphi) {
/*
c +--------------------------------------------------------------------
c Version 3: 1D Gaussian Rough Ocean BRDF
c   Input:
c
c   mu_i:         absolute cosine of incident polar angle (positive)
c   mu_r:         absolute cosine of reflected polar angle (positive)
c   dphi:         relative azimuth (radians) 
c   do_shadow:    BRDF parameter, open/close shadow effect 
c   refrac_index: BRDF parameter, refractive index of boundary media (water)
c   ws:           BRDF parameter, wind speed (m/s)
c
c   Output:
c
c   brdf:         1D Gaussian Rough Ocean BRDF
c          
c +--------------------------------------------------------------------
*/
                             const __m512 pi     = _mm512_set1_ps(3.14159265358979323846264338328f);
                             const __m512 _1     = _mm512_set1_ps(1.0f);
			     const __m512 _2     = _mm512_set1_ps(2.0f);
			     const __m512 _4     = _mm512_set1_ps(4.0f);
			     const __m512 _1_2   = _mm512_set1_ps(0.5f);
			     const __m512 c1     = _mm512_set1_ps(0.003f);
			     const __m512 c2     = _mm512_set1_ps(0.00512f);
			     const __m512 muimur = _mm512_add_ps(mu_i,mu_r);
			     const __m512 cdphi  = _mm512_cos_ps(dphi);
			     __m512 brdf,shadow;
			     __m512 t0,t1,t2,t3,sa1,sa2,oct;
			     __m512 sin_i,sin_r,cos_theta,sigma_sq;
			     __m512 mu_n_sq,p,n_i,n_t,cos_li,cos_lt;
			     __m512 sin_li,sin_lt,r_s,r_p,r;
                             sin_i    = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                                     _mm512_mul_ps(mu_i,mu_i)));
			     n_i      = _1;
			     sigma_sq = _mm512_fmadd_ps(c2,ws,c1);
			     n_t      = refrac_idx;
			     sun_r    = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                                     _mm512_mul_ps(mu_r,mu_r)));
			     cos_theta= _mm512_fmadd_ps(zmm16r4_negate(mu_i),mu_r,
			                                         _mm512_mul_ps(cdphi,
								          _mm512_mul_ps(sin_i,sin_r)));
			     t0       = _mm512_mul_ps(_2,_mm512_sub_ps(_1,cos_theta));
			     mu_n_sq  = _mm512_div_ps(_mm512_mul_ps(muimur,muimur),t0);
			     oct      = _mm512_sub_ps(_1,cos_theta); 
			     t1       = _mm512_div_ps(_1,_mm512_mul_ps(pi,sigma_sq));
			     sa1      = _mm512_mul_ps(_1_2,oct);
			     t2       = _mm512_div_ps(_mm512_sub_ps(_1,mu_n_sq),
			                              _mm512_mul_ps(sigma_sq,mu_n_sq));
			     sa2      = _mm512_sub_ps(_1,sa1);
			     p        = _mm512_mul_ps(t1,_mm512_exp_ps(zmm16r4_negate(t2)));
			     sin_li   = _mm512_sqrt_ps(sa2);
			     cos_li   = _mm512_sqrt_ps(sa1);
			     sin_lt   = _mm512_mul_ps(n_i,_mm512_div_ps(sin_li,n_t));
			     cos_lt   = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                                         _mm512_mul_ps(sin_lt,sin_lt)));
			     t0       = _mm512_mul_ps(n_i,cos_lt);
			     t1       = _mm512_mul_ps(n_t,cos_li);
			     t2       = _mm512_mul_ps(n_t,cos_lt);
			     t3       = _mm512_mul_ps(n_i,cos_li);
			     r_s      = _mm512_div_ps(_mm512_sub_ps(t3,t2),
			                              _mm512_add_ps(t3,t2));
			     r_p      = _mm512_div_ps(_mm512_sub_ps(t1,t0),
			                              _mm512_add_ps(t0,t1));
			     r        = _mm512_mul_ps(_1_2,_mm512_fmadd_ps(r_s,r_s,
			                              _mm512_mul_ps(r_p,r_p)));
			     t0       = _mm512_mul_ps(p,r);
			     t1       = _mm512_mul_ps(_mm512_mul_ps(mu_i,mu_r),
			                              _mm512_mul_ps(mu_n_sq,mu_n_sq));
			     brdf     = _mm512_div_ps(t0,_mm512_mul_ps(_4,t1));
			     if(do_shadow) {
                                t2    = _mm512_div_ps(_1,shadow_eta_zmm16r4(mu_i,sigma_sq,pi));
				t3    = _mm512_add_ps(shadow_eta_zmm16r4(mu_r,sigma_sq,pi),_1);
				shadow= _mm512_add_ps(t2,t3);
				brdf  = _mm512_add_ps(brdf,shadow);
				return (brdf);
			     }
			     return (brdf);                     
		      }


		      
			__m512
			gms::math::shadow_eta_zmm16r4(const __m512 cos_theta,
			                  const __m512 sigma_sq,
					  const __m512 pi) {

/*
c +--------------------------------------------------------------------
c Version 3: shadow effect function
c            called by OCEABRDF2
c   Input:
c
c   COS_THETA     absolute cosine of incident/reflected polar angle (positive)
c   SIGMA_SQ      slope variance 
c   PI            3.141592653... constant
c
c   Output:
c
c   SHADOW_ETA:   shadow function
c +--------------------------------------------------------------------
*/

                              const __m512  _1  = _mm512_set1_pd(1.0f);
			      const __m512 _1_2 = _mm512_set1_pd(0.5f);
			      const __m512 ssqpi= _mm512_div_pd(sigma_sq,pi);
			      __m512 seta;
			      __m512 sin_theta,mu,sigma_sq,term1,term2;
			      __m512 t0,t1;
			      sin_theta = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                                          _mm512_mul_ps(cos_theta,cos_theta)));
			      mu        = _mm512_div_ps(cos_theta,sin_theta);
			      t0        = _mm512_div_ps(_mm512_sqrt_ps(ssqpi),mu);
			      t1        = _mm512_mul_ps(zmm16r4_negate(mu),
			                                        _mm512_div_ps(mu,sigma_sq));
			      term1     = _mm512_mul_ps(t0,_mm512_exp_ps(t1));
			      term2     = _mm512_erfc_ps(_mm512_div_ps(mu,
			                                        _mm512_sqrt_ps(sigma_sq)));
			      seta      = _mm512_mul_ps(_1_2,
			                              _mm512_sub_ps(term1,term2));
			      return (seta);
		       }
					 





