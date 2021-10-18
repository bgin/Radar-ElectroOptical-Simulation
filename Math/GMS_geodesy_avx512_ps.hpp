
#ifndef __GMS_GEODESY_AVX512_PS_HPP__
#define __GMS_GEODESY_AVX512_PS_HPP__ 181020210942


namespace file_version {

    const unsigned int gGMS_GEODESY_AVX512_PS_MAJOR = 1U;
    const unsigned int gGMS_GEODESY_AVX512_PS_MINOR = 0U;
    const unsigned int gGMS_GEODESY_AVX512_PS_MICRO = 0U;
    const unsigned int gGMS_GEODESY_AVX512_PS_FULLVER =
      1000U*gGMS_GEODESY_AVX512_PS_MAJOR+
      100U*gGMS_GEODESY_AVX512_PS_MINOR+
      10U*gGMS_GEODESY_AVX512_PS_MICRO;
    const char * const pgGMS_GEODESY_AVX512_PS_CREATION_DATE = "18-10-2021 09:42 AM +00200 (MON 18 OCT 2021 GMT+2)";
    const char * const pgGMS_GEODESY_AVX512_PS_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_GEODESY_AVX512_PS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_GEODESY_AVX512_PS_DESCRIPTION   = "Vectorized (AVX512 packed single-precision) geodesic computation implementation."

}


#include <immintrin.h>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include "GMS_config.h"


namespace  gms {

          namespace math {


	                 /*
                              Cartesian to geodetic conversion (kernel).
                              
                          */
	                __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                static inline
			void
			cart_to_geodetic_zmm16r4(    const __m512 pos_x, //input position x [km]
			                             const __m512 pos_y, //input position y [km]
						     const __m512 pos_z, //input position z [km]
						     const __m512 a, // semi-minor axis [km]
						     const __m512 b, // semi-major axis [km]
						     __m512 &alt, //output altitude [km]
						     __m512 &lon, //output longtitude [rad]
						     __m512 &lat) { // output latitiude [km]
                            const __m512 _0   = _mm512_setzero_ps();
			    const __m512 _1   = _mm512_set1_ps(1.0F);
			    const __m512 _2   = _mm512_set1_ps(1.0F);
			    const __m512 _54  = _mm512_set1_ps(54.0F);
			    const __m512 _3   = _mm512_set1_ps(3.0F);
			    const __m512 _0_5 = _mm512_set1_ps(0.5F);
			    const __m512 _0_3 = _mm512_set1_ps(0.3333333333333333333333333333333F);
			    __m512 vaa = _0;
			    __m512 vbb = _0;
			    __m512 vf  = _0;
			    __m512 ve2 = _0;
			    __m512 vep = _0;
			    __m512 vzz = _0;
			    __m512 vr  = _0;
                            __m512 vrr = _0;
			    __m512 vee = _0;
			    __m512 vff = _0;
			    __m512 vg  = _0;
			    __m512 vc  = _0;
			    __m512 vs  = _0;
			    __m512 vpp = _0;
			   
			    __m512 vq  = _0;
			    __m512 vr0 = _0;
			    __m512 vt0 = _0;
			    __m512 vt1 = _0;
			    __m512 vt2 = _0;
			    __m512 vt3 = _0;
			    __m512 vt4 = _0;
			    __m512 vt5 = _0;
			    __m512 vt6 = _0;
			    __m512 vt7 = _0;
			    __m512 vt8 = _0;
			    __m512 vt9 = _0;
			    __m512 vt10 = _0;
			    __m512 vt11 = _0;
			    __m512 vt12 = _0;
			    __m512 vt13 = _0;
			    __m512 vt14 = _0;
			    __m512 vt15 = _0;
			    __m512 vu  = _0;
			    __m512 vv  = _0;
			    __m512 vz0 = _0;
			   			    
			    vf  = _mm512_div_ps(_mm512_sub_ps(a,b),a);
			    vaa = _mm512_mul_ps(a,a);
			    ve2 = _mm512_sub_ps(_mm512_mul_ps(_2,vf),
			                        _mm512_mul_ps(vf,vf));
			    vbb = _mm512_mul_ps(b,b);
			    vzz = _mm512_mul_ps(pos_z,pos_z);
			    ve2 = _mm512_sub_ps(aa,bb);
			    vr  = _mm512_sqrt_ps(_mm512_mul_ps(pos_x,pos_x),
			                         _mm512_mul_ps(pos_y,pos_y));
			    vrr = _mm512_mul_ps(vr,vr);
			    vff = _mm512_mul_ps(_54,
			                        _mm512_mul_ps(vbb,vzz));
			    vep = _mm512_sqrt_ps(_mm512_sub_ps(
			                         _mm512_div_ps(vaa,vbb),_1));
			    vt0 = _mm512_add_ps(vrr,_mm512_mul_ps(
			                        _mm512_sub_ps(_1,ve2),vzz));
			    vg  = _mm512_sub_ps(vt0,_mm512_mul_pd(ve2,vee));
			    const __m512 vg2 = _mm512_mul_ps(vg,vg);
			    const __m512 vg3  = _mm512_mul_ps(vg2,vg);
			    const __m512 ve22 = _mm512_mul_ps(ve2,ve2);
			    vc = _mm512_div_ps(_mm512_mul_ps(ve22,
			                       _mm512_mul_ps(vff,vrr)),vg3);
			    vt1 = _mm512_add_ps(_1,vc);
			    vt2 = _mm512_mul_ps(vc,vc);
			    vs  = _mm512_pow_ps(_mm512_add_ps(vt1,
			                     _mm512_sqrt_ps(
					     _mm512_fmadd_ps(_2,vc,vt2))),_0_3);
			    vt3 = _mm512_add_ps(vs,_mm512_add_ps(
			                           _mm512_div_ps(_1,vs),_1));
                            vt4 = _mm512_mul_ps(_mm512_mul_ps(_3,
			                        _mm512_mul_ps(vt3,vt3)),vg2);
			    vpp = _mm512_div_ps(ff,vt4);
			    vt5 = _mm512_mul_ps(ve22,vpp);
			    vq  = _mm512_sqrt_ps(_mm512_fmadd_pd(_2,vt5,_1));
			    vt6 = _mm512_sub_ps(_0,vpp);
			    vt7 = _mm512_div_ps(_mm512_mul_ps(vt6,
			                        _mm512_mul_ps(ve2,vr)),
						_mm512_add_ps(_1,vq));
			    vt8 = _mm512_mul_ps(_0_5,vaa);
			    vt9 = _mm512_mul_ps(vt8,_mm512_add_ps(_1,
			                                          _mm512_div_ps(_1,vq)));
			    vt10 = _mm512_mul_ps(vpp,_mm512_mul_ps(
			                                      _mm512_sub_ps(_1,ve2),vzz));
			    vt11 = _mm512_mul_ps(_1,_mm512_add_ps(_1,vq));
			    vt12 = _mm512_mul_ps(vpp,_mm512_mul_ps(
			                             _mm512_sub_ps(_1,ve2),vzz));
			    vt13 = _mm512_div_ps(vt12,_mm512_mul_ps(vq,
			                              _mm512_add_ps(_1,vq)));
			    vt14 = _mm512_mul_ps(_0_5,_mm512_mul_ps(vpp,vrr));
			    vt15 = _mm512_sub_ps(vt9,_mm512_sub_ps(vt13,vt14))
			    vr0  = _mm512_add_ps(vt7,_mm512_sqrt_ps(_mm512_max_ps(_0,vt15)));
			    vt15 = _mm512_sub_ps(vr,_mm512_mul_ps(ve2,vr0));
			    const __m512d vt16 = _mm512_mul_ps(vt15,vt15);
			    vu = _mm512_sqrt_ps(_mm512_add_ps(vt16,vzz));
			    vv = _mm512_sqrt_ps(_mm512_add_ps(vt16,
			                        _mm512_mul_ps(_mm512_sub_ps(_1,ve2),vzz)));
			    vt14 = _mm512_mul_ps(_mm512_mul_ps(vbb,vbb),vz);
			    vt13 = _mm512_mul_ps(a,vv);
			    vz0  = _mm512_div_ps(vt14,vt13);
			    vt12 = _mm512_sub_ps(_1,_mm512_div_ps(vbb,vt13));
			    alt = _mm512_mul_ps(vu,vt12);
			    lat = _mm512_atan2_ps(_mm512_fmadd_ps(
			                          _mm512_mul_ps(vep,vep),vz0,pos_z),vr);
			    lon = _mm512_atan2_ps(pos_y,pos_x);
			}


			
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			cart_to_geodetic_u_zmm16r4_looped(const float * __restrict pos_x,
			                                  const float * __restrict pos_y,
							  const float * __restrict pos_z,
							  const float a,
							  const float b,
							  float * __restrict alt,
							  float * __restrict lon,
							  float * __restrict lat,
							  const int32_t n) {

			      if(__builtin_expect(n<=0,0)) {return;}
			         const __m512 _0   = _mm512_setzero_ps();
			         const __m512 _1   = _mm512_set1_ps(1.0F);
			         const __m512 _2   = _mm512_set1_ps(1.0F);
			         const __m512 _54  = _mm512_set1_ps(54.0F);
			         const __m512 _3   = _mm512_set1_ps(3.0F);
			         const __m512 _0_5 = _mm512_set1_ps(0.5F);
			         const __m512 _0_3 = _mm512_set1_ps(0.3333333333333333333333333333333F);
				 const __m512 va   = _mm512_set1_ps(a);
				 const __m512 vb   = _mm512_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			         for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
				 
                                      __m512  vf  = _mm512_div_ps(_mm512_sub_ps(va,vb),va);
			              __m512  vaa = _mm512_mul_ps(va,va);
			              __m512  ve2 = _mm512_sub_ps(_mm512_mul_ps(_2,vf),
			                                           _mm512_mul_ps(vf,vf));
			              __m512  vbb = _mm512_mul_ps(vb,vb);
				       _mm_prefetch((const char*)&pos_z[i+16],_MM_HINT_T0);
				      const __m512 vz = _mm512_loadu_ps(&pos_z[i]);
			              __m512  vzz = _mm512_mul_ps(vz,vz);
			              __m512  ve2 = _mm512_sub_ps(aa,bb);
				      _mm_prefetch((const char*)&pos_x[i+16],_MM_HINT_T0);
				      _mm_prefetch((const char*)&pos_y[i+16],_MM_HINT_T0);
				      const __m512 vx = _mm512_loadu_ps(&pos_x[i]);
				      const __m512 vy = _mm512_loadu_ps(&pos_y[i]);
			              __m512 vr  = _mm512_sqrt_ps(_mm512_mul_ps(vx,vx),
			                                            _mm512_mul_ps(vy,vy));
			              __m512  vrr = _mm512_mul_ps(vr,vr);
			              __m512  vff = _mm512_mul_ps(_54,
			                                   _mm512_mul_ps(vbb,vzz));
			              __m512  vep = _mm512_sqrt_ps(_mm512_sub_ps(
			                                   _mm512_div_ps(vaa,vbb),_1));
			              __m512  vt0 = _mm512_add_ps(vrr,_mm512_mul_ps(
			                                    _mm512_sub_ps(_1,ve2),vzz));
			              __m512   vg  = _mm512_sub_ps(vt0,_mm512_mul_ps(ve2,vee));
			              const __m512 vg2 = _mm512_mul_ps(vg,vg);
			              const __m512 vg3  = _mm512_mul_ps(vg2,vg);
			              const __m512 ve22 = _mm512_mul_ps(ve2,ve2);
			              __m512   vc = _mm512_div_ps(_mm512_mul_pd(ve22,
			                                     _mm512_mul_ps(vff,vrr)),vg3);
			              __m512   vt1 = _mm512_add_ps(_1,vc);
			              __m512   vt2 = _mm512_mul_ps(vc,vc);
			              __m512   vs  = _mm512_pow_ps(_mm512_add_ps(vt1,
			                                        _mm512_sqrt_ps(
					                           _mm512_fmadd_ps(_2,vc,vt2))),_0_3);
			              __m512   vt3 = _mm512_add_ps(vs,_mm512_add_ps(
			                                         _mm512_div_ps(_1,vs),_1));
                                      __m512   vt4 = _mm512_mul_ps(_mm512_mul_ps(_3,
			                                    _mm512_mul_ps(vt3,vt3)),vg2);
			              __m512   vpp = _mm512_div_ps(ff,vt4);
			              __m512   vt5 = _mm512_mul_ps(ve22,vpp);
			              __m512   vq  = _mm512_sqrt_ps(_mm512_fmadd_ps(_2,vt5,_1));
			              __m512   vt6 = _mm512_sub_ps(_0,vpp);
			              __m512   vt7 = _mm512_div_ps(_mm512_mul_ps(vt6,
			                                              _mm512_mul_ps(ve2,vr)),
						                      _mm512_add_ps(_1,vq));
			              __m512   vt8 = _mm512_mul_ps(_0_5,vaa);
			              __m512   vt9 = _mm512_mul_ps(vt8,_mm512_add_ps(_1,
			                                          _mm512_div_ps(_1,vq)));
			              __m512   vt10 = _mm512_mul_ps(vpp,_mm512_mul_ps(
			                                          _mm512_sub_ps(_1,ve2),vzz));
			              __m512   vt11 = _mm512_mul_ps(_1,_mm512_add_ps(_1,vq));
			              __m512   vt12 = _mm512_mul_ps(vpp,_mm512_mul_ps(
			                                         _mm512_sub_ps(_1,ve2),vzz));
			              __m512   vt13 = _mm512_div_ps(vt12,_mm512_mul_ps(vq,
			                                         _mm512_add_ps(_1,vq)));
			              __m512   vt14 = _mm512_mul_ps(_0_5,_mm512_mul_ps(vpp,vrr));
			              __m512   vt15 = _mm512_sub_ps(vt9,_mm512_sub_ps(vt13,vt14))
			              __m512   vr0  = _mm512_add_ps(vt7,_mm512_sqrt_ps(_mm512_max_ps(_0,vt15)));
			                               vt15 = _mm512_sub_ps(vr,_mm512_mul_ps(ve2,vr0));
			              const __m512 vt16 = _mm512_mul_ps(vt15,vt15);
			              __m512   vu = _mm512_sqrt_ps(_mm512_add_ps(vt16,vzz));
			              __m512   vv = _mm512_sqrt_ps(_mm512_add_ps(vt16,
			                                       _mm512_mul_ps(_mm512_sub_ps(_1,ve2),vzz)));
			                        vt14 = _mm512_mul_ps(_mm512_mul_ps(vbb,vbb),vz);
			                        vt13 = _mm512_mul_ps(va,vv);
			              __m512   vz0  = _mm512_div_ps(vt14,vt13);
			                        vt12 = _mm512_sub_ps(_1,_mm512_div_ps(vbb,vt13));
				       __mm512_storeu_ps(&alt[i], _mm512_mul_ps(vu,vt12));
			               __m512_storeu_ps(&lat[i],   _mm512_atan2_ps(_mm512_fmadd_ps(
			                          _mm512_mul_ps(vep,vep),vz0,vz),vr));
			               __m512_storeu_ps(&lon[i], _mm512_atan2_ps(vy,vx));
				 }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#endif
                                    for(; i != n; ++i) {

				          float aa = a*a;
					  float bb = b*b;
					  float f  = (a-b)/a;
					  float e2 = (2.0F*f-f*f);
					  float ep = std::sqrtf(aa/bb-1.0F);
					  const float z  = pos_z[i];
					  float zz = z*z;
					  const float x  = pos_x[i];
					  const float y  = pos_y[i];
					  float r  = std::sqrtf(x*x+y*y);
					  float rr = r*r;
					  float ee = aa-bb;
					  float ff = 54.0*bb*zz;
					  float g  = rr+(1.0-e2)*zz-e2*ee;
					  float c  = e2*e2*ff*rr/(g*g*g);
					  float t0 = 1.0+c;
					  float s  = std::powf(t0+std::sqrtf(c*c+2.0*c),0.3333333333333333333333333333333F);
					  float t1 = (s+1.0F/s+1.0F)*(s+1.0F/s+1.0F);
					  float pp = ff/(3.0F*t1*g*g);
					  float t2 = e2*e2*pp;
					  float q  = std::sqrtf(1.0F+2.0F*t2);
					  float t3 = -pp*e2*r/(1.0F+q);
					  float t4 = 0.5*aa*(1.0F+1.0F/q);
					  float t5 = pp*(1.0F-e2)*zz;
					  float t6 = q*(1.0F+q);
					  float t7 = 0.5*pp*rr;
					  float r0 = t3+std::sqrtf(std::max(0.0,t4-(t5/t6)-t7));
					  const float t8 = (r-e2*r0)*(r-e2*r0);
					  float u  = std::sqrtf(t0+zz);
					  float v  = std::sqrtf(t0+(1.0F-e2)*zz);
					  const float t9 = a*v;
					  float z0 = b*b*z/t9;
					  lat[i] = std::atan2((z+ep*ep*z0),r);
					  alt[i] = u*(1.0F-b2/t9);
					  lon[i] = std::atan2(y,x);
				    }
				 
			}



			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			cart_to_geodetic_a_zmm1r4_looped(float * __restrict __ATTR_ALIGN__(64) pos_x,
			                                 float * __restrict __ATTR_ALIGN__(64) pos_y,
							 float * __restrict __ATTR_ALIGN__(64) pos_z,
							 const float a,
							 const float b,
							 float * __restrict __ATTR_ALIGN__(64) alt,
							 float * __restrict __ATTR_ALIGN__(64) lon,
							 float * __restrict __ATTR_ALIGN__(64) lat,
							 const int32_t n) {

			      if(__builtin_expect(n<=0,0)) {return;}
			         const __m512 _0   = _mm512_setzero_ps();
			         const __m512 _1   = _mm512_set1_ps(1.0F);
			         const __m512 _2   = _mm512_set1_ps(1.0F);
			         const __m512 _54  = _mm512_set1_ps(54.0F);
			         const __m512 _3   = _mm512_set1_ps(3.0F);
			         const __m512 _0_5 = _mm512_set1_ps(0.5F);
			         const __m512 _0_3 = _mm512_set1_ps(0.3333333333333333333333333333333F);
				 const __m512 va   = _mm512_set1_ps(a);
				 const __m512 vb   = _mm512_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                                 __assume_aligned(pos_x,64);
				 __assume_aligned(pos_y,64);
				 __assume_aligned(pos_z,64);
				 __assume_aligned(alt,64);
				 __assume_aligned(lon,64);
				 __assume_aligned(lat,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                                 pos_x = (float*)__builtin_assume_aligned(pos_x,64);
				 pos_y = (float*)__builtin_assume_aligned(pos_y,64);
				 pos_z = (float*)__builtin_assume_aligned(pos_z,64);
				 alt   = (float*)__builtin_assume_aligned(alt,64);
				 lon   = (float*)__builtin_assume_aligned(lon,64);
				 lat   = (float*)__builtin_assume_aligned(lat,64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			         for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
				 
                                      __m512  vf  = _mm512_div_ps(_mm512_sub_ps(va,vb),va);
			              __m512  vaa = _mm512_mul_ps(va,va);
			              __m512  ve2 = _mm512_sub_ps(_mm512_mul_ps(_2,vf),
			                                           _mm512_mul_ps(vf,vf));
			              __m512  vbb = _mm512_mul_ps(vb,vb);
				       _mm_prefetch((const char*)&pos_z[i+16],_MM_HINT_T0);
				      const __m512 vz = _mm512_load_ps(&pos_z[i]);
			              __m512  vzz = _mm512_mul_ps(vz,vz);
			              __m512  ve2 = _mm512_sub_ps(aa,bb);
				      _mm_prefetch((const char*)&pos_x[i+16],_MM_HINT_T0);
				      _mm_prefetch((const char*)&pos_y[i+16],_MM_HINT_T0);
				      const __m512 vx = _mm512_load_ps(&pos_x[i]);
				      const __m512 vy = _mm512_load_ps(&pos_y[i]);
			              __m512 vr  = _mm512_sqrt_ps(_mm512_mul_ps(vx,vx),
			                                            _mm512_mul_ps(vy,vy));
			              __m512  vrr = _mm512_mul_ps(vr,vr);
			              __m512  vff = _mm512_mul_ps(_54,
			                                   _mm512_mul_ps(vbb,vzz));
			              __m512  vep = _mm512_sqrt_ps(_mm512_sub_ps(
			                                   _mm512_div_ps(vaa,vbb),_1));
			              __m512  vt0 = _mm512_add_ps(vrr,_mm512_mul_ps(
			                                    _mm512_sub_ps(_1,ve2),vzz));
			              __m512   vg  = _mm512_sub_ps(vt0,_mm512_mul_ps(ve2,vee));
			              const __m512 vg2 = _mm512_mul_ps(vg,vg);
			              const __m512 vg3  = _mm512_mul_ps(vg2,vg);
			              const __m512 ve22 = _mm512_mul_ps(ve2,ve2);
			              __m512   vc = _mm512_div_ps(_mm512_mul_pd(ve22,
			                                     _mm512_mul_ps(vff,vrr)),vg3);
			              __m512   vt1 = _mm512_add_ps(_1,vc);
			              __m512   vt2 = _mm512_mul_ps(vc,vc);
			              __m512   vs  = _mm512_pow_ps(_mm512_add_ps(vt1,
			                                        _mm512_sqrt_ps(
					                           _mm512_fmadd_ps(_2,vc,vt2))),_0_3);
			              __m512   vt3 = _mm512_add_ps(vs,_mm512_add_ps(
			                                         _mm512_div_ps(_1,vs),_1));
                                      __m512   vt4 = _mm512_mul_ps(_mm512_mul_ps(_3,
			                                    _mm512_mul_ps(vt3,vt3)),vg2);
			              __m512   vpp = _mm512_div_ps(ff,vt4);
			              __m512   vt5 = _mm512_mul_ps(ve22,vpp);
			              __m512   vq  = _mm512_sqrt_ps(_mm512_fmadd_ps(_2,vt5,_1));
			              __m512   vt6 = _mm512_sub_ps(_0,vpp);
			              __m512   vt7 = _mm512_div_ps(_mm512_mul_ps(vt6,
			                                              _mm512_mul_ps(ve2,vr)),
						                      _mm512_add_ps(_1,vq));
			              __m512   vt8 = _mm512_mul_ps(_0_5,vaa);
			              __m512   vt9 = _mm512_mul_ps(vt8,_mm512_add_ps(_1,
			                                          _mm512_div_ps(_1,vq)));
			              __m512   vt10 = _mm512_mul_ps(vpp,_mm512_mul_ps(
			                                          _mm512_sub_ps(_1,ve2),vzz));
			              __m512   vt11 = _mm512_mul_ps(_1,_mm512_add_ps(_1,vq));
			              __m512   vt12 = _mm512_mul_ps(vpp,_mm512_mul_ps(
			                                         _mm512_sub_ps(_1,ve2),vzz));
			              __m512   vt13 = _mm512_div_ps(vt12,_mm512_mul_ps(vq,
			                                         _mm512_add_ps(_1,vq)));
			              __m512   vt14 = _mm512_mul_ps(_0_5,_mm512_mul_ps(vpp,vrr));
			              __m512   vt15 = _mm512_sub_ps(vt9,_mm512_sub_ps(vt13,vt14))
			              __m512   vr0  = _mm512_add_ps(vt7,_mm512_sqrt_ps(_mm512_max_ps(_0,vt15)));
			                               vt15 = _mm512_sub_ps(vr,_mm512_mul_ps(ve2,vr0));
			              const __m512 vt16 = _mm512_mul_ps(vt15,vt15);
			              __m512   vu = _mm512_sqrt_ps(_mm512_add_ps(vt16,vzz));
			              __m512   vv = _mm512_sqrt_ps(_mm512_add_ps(vt16,
			                                       _mm512_mul_ps(_mm512_sub_ps(_1,ve2),vzz)));
			                        vt14 = _mm512_mul_ps(_mm512_mul_ps(vbb,vbb),vz);
			                        vt13 = _mm512_mul_ps(va,vv);
			              __m512   vz0  = _mm512_div_ps(vt14,vt13);
			                        vt12 = _mm512_sub_ps(_1,_mm512_div_ps(vbb,vt13));
				       __mm512_store_ps(&alt[i], _mm512_mul_ps(vu,vt12));
			               __m512_store_ps(&lat[i],   _mm512_atan2_ps(_mm512_fmadd_ps(
			                          _mm512_mul_ps(vep,vep),vz0,vz),vr));
			               __m512_store_ps(&lon[i], _mm512_atan2_ps(vy,vx));
			         
			            
				 }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#endif
                                    for(; i != n; ++i) {

				          float aa = a*a;
					  float bb = b*b;
					  float f  = (a-b)/a;
					  float e2 = (2.0F*f-f*f);
					  float ep = std::sqrtf(aa/bb-1.0F);
					  const float z  = pos_z[i];
					  float zz = z*z;
					  const float x  = pos_x[i];
					  const float y  = pos_y[i];
					  float r  = std::sqrtf(x*x+y*y);
					  float rr = r*r;
					  float ee = aa-bb;
					  float ff = 54.0*bb*zz;
					  float g  = rr+(1.0-e2)*zz-e2*ee;
					  float c  = e2*e2*ff*rr/(g*g*g);
					  float t0 = 1.0+c;
					  float s  = std::powf(t0+std::sqrtf(c*c+2.0*c),0.3333333333333333333333333333333F);
					  float t1 = (s+1.0F/s+1.0F)*(s+1.0F/s+1.0F);
					  float pp = ff/(3.0F*t1*g*g);
					  float t2 = e2*e2*pp;
					  float q  = std::sqrtf(1.0F+2.0F*t2);
					  float t3 = -pp*e2*r/(1.0F+q);
					  float t4 = 0.5*aa*(1.0F+1.0F/q);
					  float t5 = pp*(1.0F-e2)*zz;
					  float t6 = q*(1.0F+q);
					  float t7 = 0.5*pp*rr;
					  float r0 = t3+std::sqrtf(std::max(0.0,t4-(t5/t6)-t7));
					  const float t8 = (r-e2*r0)*(r-e2*r0);
					  float u  = std::sqrtf(t0+zz);
					  float v  = std::sqrtf(t0+(1.0F-e2)*zz);
					  const float t9 = a*v;
					  float z0 = b*b*z/t9;
					  lat[i] = std::atan2((z+ep*ep*z0),r);
					  alt[i] = u*(1.0F-b2/t9);
					  lon[i] = std::atan2(y,x);
				    }
				 
			}


                    	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                static inline
			void geodetic_to_cart_zmm16r4(const __m512 a,
			                              const __m512 b,
						      const __m512 lat,
						      const __m512 lon,
						      const __m512 alt,
						      __m512 &pos_x,
						      __m512 &pos_y,
						      __m512 &pos_z) {
 
                              register const __m512 _0 = _mm512_setzero_pd();
			      register const __m512 _1 = _mm512_set1_pd(1.0F);
			      register __m512 t0   = _0;
			      register __m512 t1   = _0;
			      register __m512 t2   = _0;
			      register __m512 t3   = _0;
			      register __m512 t4   = _0;
			      register __m512 zmm0 = _0; //sin(lat)
			      register __m512 zmm1 = _0; //cos(lat)
			      register __m512 zmm2 = _0; //tan(lat)
			      register __m512 zmm3 = _0; //sin(lon)
			      register __m512 zmm4 = _0; //cos(lon)
			      register __m512 ve2  = _0;
			      register __m512 vom2 = _0;
			      register __m512 vd   = _0;
			      register __m512 vq   = _0;
			      register __m512 vad  = _0;
			      register __m512 vaa  = _0;
			      register __m512 vbb  = _0;
			      zmm0  = _mm512_sin_ps(lat);
			      vaa   = _mm512_mul_ps(a,a);
			      vbb   = _mm512_mul_ps(b,b);
			      zmm1  = _mm512_cos_ps(lat);
			      ve2   = _mm512_sub_ps(_1,_mm512_div_ps(vbb,vaa));
			      vom2  = _mm512_sub_ps(_1,ve2);
			      zmm2  = _mm512_tan_ps(lat);
			      t0    = _mm512_mul_ps(zmm2,zmm2);
			      vd    = _mm512_sqrt_ps(_mm512_fmadd_ps(vom2,t0,_1));
			      t1    = _mm512_mul_ps(ve2,_mm512_mul_ps(zmm0,zmm0));
			      vq    = _mm512_sqrt_ps(_mm512_sub_ps(_1,t1));
			      zmm3  = _mm512_sin_ps(lon);
			      vad   = _mm512_div_ps(a,vd);
			      zmm4  = _mm512_cos_ps(lon);
			      t2    = _mm512_mul_ps(zmm1,zmm4);
			      pos_x = _mm512_fmadd_ps(alt,t2,_mm512_mul_ps(vad,zmm4));
			      t3    = _mm512_mul_ps(zmm3,zmm1);
			      pos_y = _mm512_fmadd_ps(alt,t3,_mm512_mul_ps(vad,zmm3));
			      t4    = _mm512_div_ps(zmm0,vq);
			      t3    = _mm512_mul_ps(a,_mm512_mul_ps(vom2,t4));
			      pos_z = _mm512_fmadd_ps(alt,zmm0,t3);
			      
			}


                    	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			geodetic_to_cart_u_zmm16r4_looped(const float a,
			                                  const float b,
							  const float * __restrict lat,
							  const float * __restrict lon,
							  const float * __restrict alt,
							  float * __restrict pos_x,
							  float * __restrict pos_y,
							  float * __restrict pos_z,
							  const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			         register const __m512 _0 = _mm512_setzero_ps();
			         register const __m512 _1 = _mm512_set1_ps(1.0F);
				 register const __m512 va   = _mm512_set1_ps(a);
				 register const __m512 vb   = _mm512_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
			            _mm_prefetch((const char*)&lat[i+16],_MM_HINT_T0);
			            register const __m512 vlat = _mm512_loadu_ps(&lat[i]);
                                    register __m512 zmm0  = _mm512_sin_ps(vlat);
			            register __m512 vaa   = _mm512_mul_ps(va,va);
			            register __m512 vbb   = _mm512_mul_ps(vb,vb);
				    register __m512 zmm1  = _mm512_cos_ps(vlat);
			            register __m512 ve2   = _mm512_sub_ps(_1,_mm512_div_ps(vbb,vaa));
			            register __m512 vom2  = _mm512_sub_ps(_1,ve2);
			            register __m512 zmm2  = _mm512_tan_ps(vlat);
			            register __m512 t0    = _mm512_mul_ps(zmm2,zmm2);
			            register __m512 vd    = _mm512_sqrt_ps(_mm512_fmadd_ps(vom2,t0,_1));
			            register __m512 t1    = _mm512_mul_ps(ve2,_mm512_mul_ps(zmm0,zmm0));
			            register __m512 vq    = _mm512_sqrt_ps(_mm512_sub_ps(_1,t1));
			            _mm_prefetch((const char*)&lon[i+16],_MM_HINT_T0);
				    register const __m512 vlon = _mm512_loadu_ps(&lon[i]);
			            register __m512 zmm3  = _mm512_sin_ps(vlon);
			            register __m512 vad   = _mm512_div_ps(va,vd);
			            register __m512 zmm4  = _mm512_cos_ps(vlon);
			            register __m512 t2    = _mm512_mul_ps(zmm1,zmm4);
				    _mm_prefetch((const char*)&alt[i+16],_MM_HINT_T0);
				    register const __m512 valt = _mm512_loadu_ps(&alt[i]);
				    _mm512_storeu_ps(&pos_x[i],
				                     _mm512_fmadd_ps(valt,t2,_mm512_mul_ps(vad,zmm4)));
			            t3    = _mm512_mul_ps(zmm3,zmm1);
			            _mm512_storeu_ps(&pos_y[i],
				                     _mm512_fmadd_ps(valt,t3,_mm512_mul_ps(vad,zmm3)));
			            t4    = _mm512_div_ps(zmm0,vq);
			            t3    = _mm512_mul_ps(va,_mm512_mul_ps(vom2,t4));
			            _mm512_storeu_ps(&pos_z[i],
				                     _mm512_fmadd_ps(valt,zmm0,t3));
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#endif
                                    for(; i != n; ++i) {
				       
                                        register const  float glat = lat[i];
					register const  float s0 = std::sin(glat);
					register const  float glon = lon[i];
					register const  float s1 = std::cos(glat);
					register float ee  = 1.0F-(b*b)/(a*a);
					register const float s3 = std::sin(glon);
					register float om2 = 1.0F-ee;
					register float q  = std::sqrtf(1.0-ee*s0*s0);
					register const float salt = alt[i];
					register const float s4 = std::cos(glon);
					register const float s2 = std::tan(glat);
					register float d = std::sqrtf(1.0F-ee*s2*s2);
					register float ad = a/d;
					pos_x[i] = ad*s4+salt*s3*s1;
					pos_y[i] = ad*s3+salt*s4*1;
					pos_z[i] = a*om2*s0/q+salt*s0;
				    }
		        }


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		        static inline
			void
			geodetic_to_cart_a_zmm16r4_looped(const float a,
			                                  const float b,
							  float * __restrict __ATTR_ALIGN__(64) lat,
							  float * __restrict __ATTR_ALIGN__(64) lon,
							  float * __restrict __ATTR_ALIGN__(64) alt,
							  float * __restrict __ATTR_ALIGN__(64) pos_x,
							  float * __restrict __ATTR_ALIGN__(64) pos_y,
							  float * __restrict __ATTR_ALIGN__(64) pos_z,
							  const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			         register const __m512 _0 = _mm512_setzero_ps();
			         register const __m512 _1 = _mm512_set1_pd(1.0F);
				 register const __m512 va   = _mm512_set1_ps(a);
				 register const __m512 vb   = _mm512_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                                 __assume_aligned(pos_x,64);
				 __assume_aligned(pos_y,64);
				 __assume_aligned(pos_z,64);
				 __assume_aligned(alt,64);
				 __assume_aligned(lon,64);
				 __assume_aligned(lat,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                                 pos_x = (float*)__builtin_assume_aligned(pos_x,64);
				 pos_y = (float*)__builtin_assume_aligned(pos_y,64);
				 pos_z = (float*)__builtin_assume_aligned(pos_z,64);
				 alt   = (float*)__builtin_assume_aligned(alt,64);
				 lon   = (float*)__builtin_assume_aligned(lon,64);
				 lat   = (float*)__builtin_assume_aligned(lat,64);
#endif				 
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
			            _mm_prefetch((const char*)&lat[i+16],_MM_HINT_T0);
			            register const __m512 vlat = _mm512_load_ps(&lat[i]);
                                    register __m512 zmm0  = _mm512_sin_ps(vlat);
			            register __m512 vaa   = _mm512_mul_ps(va,va);
			            register __m512 vbb   = _mm512_mul_ps(vb,vb);
				    register __m512 zmm1  = _mm512_cos_ps(vlat);
			            register __m512 ve2   = _mm512_sub_ps(_1,_mm512_div_ps(vbb,vaa));
			            register __m512 vom2  = _mm512_sub_ps(_1,ve2);
			            register __m512 zmm2  = _mm512_tan_ps(vlat);
			            register __m512 t0    = _mm512_mul_ps(zmm2,zmm2);
			            register __m512 vd    = _mm512_sqrt_ps(_mm512_fmadd_ps(vom2,t0,_1));
			            register __m512 t1    = _mm512_mul_ps(ve2,_mm512_mul_ps(zmm0,zmm0));
			            register __m512 vq    = _mm512_sqrt_ps(_mm512_sub_ps(_1,t1));
			            _mm_prefetch((const char*)&lon[i+16],_MM_HINT_T0);
				    register const __m512 vlon = _mm512_load_ps(&lon[i]);
			            register __m512 zmm3  = _mm512_sin_ps(vlon);
			            register __m512 vad   = _mm512_div_ps(va,vd);
			            register __m512 zmm4  = _mm512_cos_ps(vlon);
			            register __m512 t2    = _mm512_mul_ps(zmm1,zmm4);
				    _mm_prefetch((const char*)&alt[i+16],_MM_HINT_T0);
				    register const __m512 valt = _mm512_load_ps(&alt[i]);
				    _mm512_store_ps(&pos_x[i],
				                     _mm512_fmadd_ps(valt,t2,_mm512_mul_ps(vad,zmm4)));
			            t3    = _mm512_mul_ps(zmm3,zmm1);
			            _mm512_store_ps(&pos_y[i],
				                     _mm512_fmadd_ps(valt,t3,_mm512_mul_ps(vad,zmm3)));
			            t4    = _mm512_div_ps(zmm0,vq);
			            t3    = _mm512_mul_ps(va,_mm512_mul_ps(vom2,t4));
			            _mm512_store_ps(&pos_z[i],
				                     _mm512_fmadd_ps(valt,zmm0,t3));
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#endif
                                    for(; i != n; ++i) {
				        register const  float glat = lat[i];
					register const  float s0 = std::sin(glat);
					register const  float glon = lon[i];
					register const  float s1 = std::cos(glat);
					register float ee  = 1.0F-(b*b)/(a*a);
					register const float s3 = std::sin(glon);
					register float om2 = 1.0F-ee;
					register float q  = std::sqrtf(1.0-ee*s0*s0);
					register const float salt = alt[i];
					register const float s4 = std::cos(glon);
					register const float s2 = std::tan(glat);
					register float d = std::sqrtf(1.0F-ee*s2*s2);
					register float ad = a/d;
					pos_x[i] = ad*s4+salt*s3*s1;
					pos_y[i] = ad*s3+salt*s4*1;
					pos_z[i] = a*om2*s0/q+salt*s0;
                                       
				    }
		        }


			


						     
     }

}














#endif /*__GMS_GEODESY_AVX512_PS_HPP__*/
