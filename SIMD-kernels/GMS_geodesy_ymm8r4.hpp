
#ifndef __GMS_GEODESY_YMM8R4_HPP__
#define __GMS_GEODESY_YMM8R4_HPP__ 201020210942


namespace file_version {

    const unsigned int gGMS_GEODESY_YMM8R4_MAJOR = 1U;
    const unsigned int gGMS_GEODESY_YMM8R4_MINOR = 0U;
    const unsigned int gGMS_GEODESY_YMM8R4_MICRO = 0U;
    const unsigned int gGMS_GEODESY_YMM8R4_FULLVER =
      1000U*gGMS_GEODESY_YMM8R4_MAJOR+
      100U*gGMS_GEODESY_YMM8R4_MINOR+
      10U*gGMS_GEODESY_YMM8R4_MICRO;
    const char * const pgGMS_GEODESY_YMM8R4_CREATION_DATE = "20-10-2021 09:42 AM +00200 (WED 20 OCT 2021 GMT+2)";
    const char * const pgGMS_GEODESY_YMM8R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_GEODESY_YMM8R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_GEODESY_YMM8R4_DESCRIPTION   = "Vectorized (AVX packed single-precision) geodesic computation implementation."

}


#include <immintrin.h>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include "GMS_config.h"
#include "GMS_cephes.h"


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
			cart_to_geodetic_ymm8r4(     const __m256 pos_x, //input position x [km]
			                             const __m256 pos_y, //input position y [km]
						     const __m256 pos_z, //input position z [km]
						     const __m256 a, // semi-minor axis [km]
						     const __m256 b, // semi-major axis [km]
						     __m256 &alt, //output altitude [km]
						     __m256 &lon, //output longtitude [rad]
						     __m256 &lat) { // output latitiude [km]
                            const __m256 _0   = _mm256_setzero_ps();
			    const __m256 _1   = _mm256_set1_ps(1.0F);
			    const __m256 _2   = _mm256_set1_ps(1.0F);
			    const __m256 _54  = _mm256_set1_ps(54.0F);
			    const __m256 _3   = _mm256_set1_ps(3.0F);
			    const __m256 _0_5 = _mm256_set1_ps(0.5F);
			    const __m256 _0_3 = _mm256_set1_ps(0.3333333333333333333333333333333F);
			    __m256 vaa = _0;
			    __m256 vbb = _0;
			    __m256 vf  = _0;
			    __m256 ve2 = _0;
			    __m256 vep = _0;
			    __m256 vzz = _0;
			    __m256 vr  = _0;
                            __m256 vrr = _0;
			    __m256 vee = _0;
			    __m256 vff = _0;
			    __m256 vg  = _0;
			    __m256 vc  = _0;
			    __m256 vs  = _0;
			    __m256 vpp = _0;
			   
			    __m256 vq  = _0;
			    __m256 vr0 = _0;
			    __m256 vt0 = _0;
			    __m256 vt1 = _0;
			    __m256 vt2 = _0;
			    __m256 vt3 = _0;
			    __m256 vt4 = _0;
			    __m256 vt5 = _0;
			    __m256 vt6 = _0;
			    __m256 vt7 = _0;
			    __m256 vt8 = _0;
			    __m256 vt9 = _0;
			    __m256 vt10 = _0;
			    __m256 vt11 = _0;
			    __m256 vt12 = _0;
			    __m256 vt13 = _0;
			    __m256 vt14 = _0;
			    __m256 vt15 = _0;
			    __m256 vu  = _0;
			    __m256 vv  = _0;
			    __m256 vz0 = _0;
			   			    
			    vf  = _mm256_div_ps(_mm256_sub_ps(a,b),a);
			    vaa = _mm256_mul_ps(a,a);
			    ve2 = _mm256_sub_ps(_mm256_mul_ps(_2,vf),
			                        _mm256_mul_ps(vf,vf));
			    vbb = _mm256_mul_ps(b,b);
			    vzz = _mm256_mul_ps(pos_z,pos_z);
			    ve2 = _mm256_sub_ps(aa,bb);
			    vr  = _mm256_sqrt_ps(_mm256_mul_ps(pos_x,pos_x),
			                         _mm256_mul_ps(pos_y,pos_y));
			    vrr = _mm256_mul_ps(vr,vr);
			    vff = _mm256_mul_ps(_54,
			                        _mm256_mul_ps(vbb,vzz));
			    vep = _mm256_sqrt_ps(_mm256_sub_ps(
			                         _mm256_div_ps(vaa,vbb),_1));
			    vt0 = _mm256_add_ps(vrr,_mm256_mul_ps(
			                        _mm256_sub_ps(_1,ve2),vzz));
			    vg  = _mm256_sub_ps(vt0,_mm256_mul_ps(ve2,vee));
			    const __m256 vg2 = _mm256_mul_ps(vg,vg);
			    const __m256 vg3  = _mm256_mul_ps(vg2,vg);
			    const __m256 ve22 = _mm256_mul_ps(ve2,ve2);
			    vc = _mm256_div_ps(_mm256_mul_ps(ve22,
			                       _mm256_mul_ps(vff,vrr)),vg3);
			    vt1 = _mm256_add_ps(_1,vc);
			    vt2 = _mm256_mul_ps(vc,vc);
			    vs  = _mm256_pow_ps(_mm256_add_ps(vt1,
			                     _mm256_sqrt_ps(
					     _mm256_fmadd_ps(_2,vc,vt2))),_0_3);
			    vt3 = _mm256_add_ps(vs,_mm256_add_ps(
			                           _mm256_div_ps(_1,vs),_1));
                            vt4 = _mm256_mul_ps(_mm256_mul_ps(_3,
			                        _mm256_mul_ps(vt3,vt3)),vg2);
			    vpp = _mm256_div_ps(ff,vt4);
			    vt5 = _mm256_mul_ps(ve22,vpp);
			    vq  = _mm256_sqrt_ps(_mm256_fmadd_ps(_2,vt5,_1));
			    vt6 = _mm256_sub_ps(_0,vpp);
			    vt7 = _mm256_div_ps(_mm256_mul_ps(vt6,
			                        _mm256_mul_ps(ve2,vr)),
						_mm256_add_ps(_1,vq));
			    vt8 = _mm256_mul_ps(_0_5,vaa);
			    vt9 = _mm256_mul_ps(vt8,_mm256_add_ps(_1,
			                                          _mm256_div_ps(_1,vq)));
			    vt10 = _mm256_mul_ps(vpp,_mm256_mul_ps(
			                                      _mm256_sub_ps(_1,ve2),vzz));
			    vt11 = _mm256_mul_ps(_1,_mm256_add_ps(_1,vq));
			    vt12 = _mm256_mul_ps(vpp,_mm256_mul_ps(
			                             _mm256_sub_ps(_1,ve2),vzz));
			    vt13 = _mm256_div_ps(vt12,_mm256_mul_ps(vq,
			                              _mm256_add_ps(_1,vq)));
			    vt14 = _mm256_mul_ps(_0_5,_mm256_mul_ps(vpp,vrr));
			    vt15 = _mm256_sub_ps(vt9,_mm256_sub_ps(vt13,vt14))
			    vr0  = _mm256_add_ps(vt7,_mm256_sqrt_ps(_mm256_max_ps(_0,vt15)));
			    vt15 = _mm256_sub_ps(vr,_mm256_mul_ps(ve2,vr0));
			    const __m256d vt16 = _mm256_mul_ps(vt15,vt15);
			    vu = _mm256_sqrt_ps(_mm256_add_ps(vt16,vzz));
			    vv = _mm256_sqrt_ps(_mm256_add_ps(vt16,
			                        _mm256_mul_ps(_mm256_sub_ps(_1,ve2),vzz)));
			    vt14 = _mm256_mul_ps(_mm256_mul_ps(vbb,vbb),vz);
			    vt13 = _mm256_mul_ps(a,vv);
			    vz0  = _mm256_div_ps(vt14,vt13);
			    vt12 = _mm256_sub_ps(_1,_mm256_div_ps(vbb,vt13));
			    alt = _mm256_mul_ps(vu,vt12);
			    lat = _mm256_atan2_ps(_mm256_fmadd_ps(
			                          _mm256_mul_ps(vep,vep),vz0,pos_z),vr);
			    lon = _mm256_atan2_ps(pos_y,pos_x);
			}


			
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			cart_to_geodetic_u_ymm8r4_looped( const float * __restrict pos_x,
			                                  const float * __restrict pos_y,
							  const float * __restrict pos_z,
							  const float a,
							  const float b,
							  float * __restrict alt,
							  float * __restrict lon,
							  float * __restrict lat,
							  const int32_t n) {

			      if(__builtin_expect(n<=0,0)) {return;}
			         const __m256 _0   = _mm256_setzero_ps();
			         const __m256 _1   = _mm256_set1_ps(1.0F);
			         const __m256 _2   = _mm256_set1_ps(1.0F);
			         const __m256 _54  = _mm256_set1_ps(54.0F);
			         const __m256 _3   = _mm256_set1_ps(3.0F);
			         const __m256 _0_5 = _mm256_set1_ps(0.5F);
			         const __m256 _0_3 = _mm256_set1_ps(0.3333333333333333333333333333333F);
				 const __m256 va   = _mm256_set1_ps(a);
				 const __m256 vb   = _mm256_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			         for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
				 
                                      __m256  vf  = _mm256_div_ps(_mm256_sub_ps(va,vb),va);
			              __m256  vaa = _mm256_mul_ps(va,va);
			              __m256  ve2 = _mm256_sub_ps(_mm256_mul_ps(_2,vf),
			                                           _mm256_mul_ps(vf,vf));
			              __m256  vbb = _mm256_mul_ps(vb,vb);
				       _mm_prefetch((const char*)&pos_z[i+32],_MM_HINT_T0);
				      const __m256 vz = _mm256_loadu_ps(&pos_z[i]);
			              __m256  vzz = _mm256_mul_ps(vz,vz);
			              __m256  ve2 = _mm256_sub_ps(aa,bb);
				      _mm_prefetch((const char*)&pos_x[i+32],_MM_HINT_T0);
				      _mm_prefetch((const char*)&pos_y[i+32],_MM_HINT_T0);
				      const __m256 vx = _mm256_loadu_ps(&pos_x[i]);
				      const __m256 vy = _mm256_loadu_ps(&pos_y[i]);
			              __m256 vr  = _mm256_sqrt_ps(_mm256_mul_ps(vx,vx),
			                                            _mm256_mul_ps(vy,vy));
			              __m256  vrr = _mm256_mul_ps(vr,vr);
			              __m256  vff = _mm256_mul_ps(_54,
			                                   _mm256_mul_ps(vbb,vzz));
			              __m256  vep = _mm256_sqrt_ps(_mm256_sub_ps(
			                                   _mm256_div_ps(vaa,vbb),_1));
			              __m256  vt0 = _mm256_add_ps(vrr,_mm256_mul_ps(
			                                    _mm256_sub_ps(_1,ve2),vzz));
			              __m256   vg  = _mm256_sub_ps(vt0,_mm256_mul_ps(ve2,vee));
			              const __m256 vg2 = _mm256_mul_ps(vg,vg);
			              const __m256 vg3  = _mm256_mul_ps(vg2,vg);
			              const __m256 ve22 = _mm256_mul_ps(ve2,ve2);
			              __m256   vc = _mm256_div_ps(_mm256_mul_ps(ve22,
			                                     _mm256_mul_ps(vff,vrr)),vg3);
			              __m256   vt1 = _mm256_add_ps(_1,vc);
			              __m256   vt2 = _mm256_mul_ps(vc,vc);
			              __m256   vs  = _mm256_pow_ps(_mm256_add_ps(vt1,
			                                        _mm256_sqrt_ps(
					                           _mm256_fmadd_ps(_2,vc,vt2))),_0_3);
			              __m256   vt3 = _mm256_add_ps(vs,_mm256_add_ps(
			                                         _mm256_div_ps(_1,vs),_1));
                                      __m256   vt4 = _mm256_mul_ps(_mm256_mul_ps(_3,
			                                    _mm256_mul_ps(vt3,vt3)),vg2);
			              __m256   vpp = _mm256_div_ps(ff,vt4);
			              __m256   vt5 = _mm256_mul_ps(ve22,vpp);
			              __m256   vq  = _mm256_sqrt_ps(_mm256_fmadd_ps(_2,vt5,_1));
			              __m256   vt6 = _mm256_sub_ps(_0,vpp);
			              __m256   vt7 = _mm256_div_ps(_mm256_mul_ps(vt6,
			                                              _mm256_mul_ps(ve2,vr)),
						                      _mm256_add_ps(_1,vq));
			              __m256   vt8 = _mm256_mul_ps(_0_5,vaa);
			              __m256   vt9 = _mm256_mul_ps(vt8,_mm256_add_ps(_1,
			                                          _mm256_div_ps(_1,vq)));
			              __m256   vt10 = _mm256_mul_ps(vpp,_mm256_mul_ps(
			                                          _mm256_sub_ps(_1,ve2),vzz));
			              __m256   vt11 = _mm256_mul_ps(_1,_mm256_add_ps(_1,vq));
			              __m256   vt12 = _mm256_mul_ps(vpp,_mm256_mul_ps(
			                                         _mm256_sub_ps(_1,ve2),vzz));
			              __m256   vt13 = _mm256_div_ps(vt12,_mm256_mul_ps(vq,
			                                         _mm256_add_ps(_1,vq)));
			              __m256   vt14 = _mm256_mul_ps(_0_5,_mm256_mul_ps(vpp,vrr));
			              __m256   vt15 = _mm256_sub_ps(vt9,_mm256_sub_ps(vt13,vt14))
			              __m256   vr0  = _mm256_add_ps(vt7,_mm256_sqrt_ps(_mm256_max_ps(_0,vt15)));
			                               vt15 = _mm256_sub_ps(vr,_mm256_mul_ps(ve2,vr0));
			              const __m256 vt16 = _mm256_mul_ps(vt15,vt15);
			              __m256   vu = _mm256_sqrt_ps(_mm256_add_ps(vt16,vzz));
			              __m256   vv = _mm256_sqrt_ps(_mm256_add_ps(vt16,
			                                       _mm256_mul_ps(_mm256_sub_ps(_1,ve2),vzz)));
			                        vt14 = _mm256_mul_ps(_mm256_mul_ps(vbb,vbb),vz);
			                        vt13 = _mm256_mul_ps(va,vv);
			              __m256   vz0  = _mm256_div_ps(vt14,vt13);
			                        vt12 = _mm256_sub_ps(_1,_mm256_div_ps(vbb,vt13));
				       __mm256_storeu_ps(&alt[i], _mm256_mul_ps(vu,vt12));
			               __m256_storeu_ps(&lat[i],   _mm256_atan2_ps(_mm256_fmadd_ps(
			                          _mm256_mul_ps(vep,vep),vz0,vz),vr));
			               __m256_storeu_ps(&lon[i], _mm256_atan2_ps(vy,vx));
				 }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {

				          float aa = a*a;
					  float bb = b*b;
					  float f  = (a-b)/a;
					  float e2 = (2.0F*f-f*f);
					  float ep = cephes_sqrtf(aa/bb-1.0F);
					  const float z  = pos_z[i];
					  float zz = z*z;
					  const float x  = pos_x[i];
					  const float y  = pos_y[i];
					  float r  = cephes_sqrtf(x*x+y*y);
					  float rr = r*r;
					  float ee = aa-bb;
					  float ff = 54.0*bb*zz;
					  float g  = rr+(1.0-e2)*zz-e2*ee;
					  float c  = e2*e2*ff*rr/(g*g*g);
					  float t0 = 1.0+c;
					  float s  = cephes_powf(t0+cephes_sqrtf(c*c+2.0*c),0.3333333333333333333333333333333F);
					  float t1 = (s+1.0F/s+1.0F)*(s+1.0F/s+1.0F);
					  float pp = ff/(3.0F*t1*g*g);
					  float t2 = e2*e2*pp;
					  float q  = cephes_sqrtf(1.0F+2.0F*t2);
					  float t3 = -pp*e2*r/(1.0F+q);
					  float t4 = 0.5*aa*(1.0F+1.0F/q);
					  float t5 = pp*(1.0F-e2)*zz;
					  float t6 = q*(1.0F+q);
					  float t7 = 0.5*pp*rr;
					  float r0 = t3+cephes_sqrtf(std::max(0.0,t4-(t5/t6)-t7));
					  const float t8 = (r-e2*r0)*(r-e2*r0);
					  float u  = cephes_sqrtf(t0+zz);
					  float v  = cephes_sqrtf(t0+(1.0F-e2)*zz);
					  const float t9 = a*v;
					  float z0 = b*b*z/t9;
					  lat[i] = cephes_atan2f((z+ep*ep*z0),r);
					  alt[i] = u*(1.0F-b2/t9);
					  lon[i] = cephes_atan2f(y,x);
				    }
				 
			}



			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			cart_to_geodetic_a_ymm8r4_looped(float * __restrict __ATTR_ALIGN__(32) pos_x,
			                                 float * __restrict __ATTR_ALIGN__(32) pos_y,
							 float * __restrict __ATTR_ALIGN__(32) pos_z,
							 const float a,
							 const float b,
							 float * __restrict __ATTR_ALIGN__(32) alt,
							 float * __restrict __ATTR_ALIGN__(32) lon,
							 float * __restrict __ATTR_ALIGN__(32) lat,
							 const int32_t n) {

			      if(__builtin_expect(n<=0,0)) {return;}
			         const __m256 _0   = _mm256_setzero_ps();
			         const __m256 _1   = _mm256_set1_ps(1.0F);
			         const __m256 _2   = _mm256_set1_ps(1.0F);
			         const __m256 _54  = _mm256_set1_ps(54.0F);
			         const __m256 _3   = _mm256_set1_ps(3.0F);
			         const __m256 _0_5 = _mm256_set1_ps(0.5F);
			         const __m256 _0_3 = _mm256_set1_ps(0.3333333333333333333333333333333F);
				 const __m256 va   = _mm256_set1_ps(a);
				 const __m256 vb   = _mm256_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                                 __assume_aligned(pos_x,32);
				 __assume_aligned(pos_y,32);
				 __assume_aligned(pos_z,32);
				 __assume_aligned(alt,32);
				 __assume_aligned(lon,32);
				 __assume_aligned(lat,32);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                                 pos_x = (float*)__builtin_assume_aligned(pos_x,32);
				 pos_y = (float*)__builtin_assume_aligned(pos_y,32);
				 pos_z = (float*)__builtin_assume_aligned(pos_z,32);
				 alt   = (float*)__builtin_assume_aligned(alt,32);
				 lon   = (float*)__builtin_assume_aligned(lon,32);
				 lat   = (float*)__builtin_assume_aligned(lat,32);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			         for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
				 
                                      __m256  vf  = _mm256_div_ps(_mm256_sub_ps(va,vb),va);
			              __m256  vaa = _mm256_mul_ps(va,va);
			              __m256  ve2 = _mm256_sub_ps(_mm256_mul_ps(_2,vf),
			                                           _mm256_mul_ps(vf,vf));
			              __m256  vbb = _mm256_mul_ps(vb,vb);
				       _mm_prefetch((const char*)&pos_z[i+32],_MM_HINT_T0);
				      const __m256 vz = _mm256_load_ps(&pos_z[i]);
			              __m256  vzz = _mm256_mul_ps(vz,vz);
			              __m256  ve2 = _mm256_sub_ps(aa,bb);
				      _mm_prefetch((const char*)&pos_x[i+32],_MM_HINT_T0);
				      _mm_prefetch((const char*)&pos_y[i+32],_MM_HINT_T0);
				      const __m256 vx = _mm256_load_ps(&pos_x[i]);
				      const __m256 vy = _mm256_load_ps(&pos_y[i]);
			              __m256 vr  = _mm256_sqrt_ps(_mm256_mul_ps(vx,vx),
			                                            _mm256_mul_ps(vy,vy));
			              __m256  vrr = _mm256_mul_ps(vr,vr);
			              __m256  vff = _mm256_mul_ps(_54,
			                                   _mm256_mul_ps(vbb,vzz));
			              __m256  vep = _mm256_sqrt_ps(_mm256_sub_ps(
			                                   _mm256_div_ps(vaa,vbb),_1));
			              __m256  vt0 = _mm256_add_ps(vrr,_mm256_mul_ps(
			                                    _mm256_sub_ps(_1,ve2),vzz));
			              __m256   vg  = _mm256_sub_ps(vt0,_mm256_mul_ps(ve2,vee));
			              const __m256 vg2 = _mm256_mul_ps(vg,vg);
			              const __m256 vg3  = _mm256_mul_ps(vg2,vg);
			              const __m256 ve22 = _mm256_mul_ps(ve2,ve2);
			              __m256   vc = _mm256_div_ps(_mm256_mul_ps(ve22,
			                                     _mm256_mul_ps(vff,vrr)),vg3);
			              __m256   vt1 = _mm256_add_ps(_1,vc);
			              __m256   vt2 = _mm256_mul_ps(vc,vc);
			              __m256   vs  = _mm256_pow_ps(_mm256_add_ps(vt1,
			                                        _mm256_sqrt_ps(
					                           _mm256_fmadd_ps(_2,vc,vt2))),_0_3);
			              __m256   vt3 = _mm256_add_ps(vs,_mm256_add_ps(
			                                         _mm256_div_ps(_1,vs),_1));
                                      __m256   vt4 = _mm256_mul_ps(_mm256_mul_ps(_3,
			                                    _mm256_mul_ps(vt3,vt3)),vg2);
			              __m256   vpp = _mm256_div_ps(ff,vt4);
			              __m256   vt5 = _mm256_mul_ps(ve22,vpp);
			              __m256   vq  = _mm256_sqrt_ps(_mm256_fmadd_ps(_2,vt5,_1));
			              __m256   vt6 = _mm256_sub_ps(_0,vpp);
			              __m256   vt7 = _mm256_div_ps(_mm256_mul_ps(vt6,
			                                              _mm256_mul_ps(ve2,vr)),
						                      _mm256_add_ps(_1,vq));
			              __m256   vt8 = _mm256_mul_ps(_0_5,vaa);
			              __m256   vt9 = _mm256_mul_ps(vt8,_mm256_add_ps(_1,
			                                          _mm256_div_ps(_1,vq)));
			              __m256   vt10 = _mm256_mul_ps(vpp,_mm256_mul_ps(
			                                          _mm256_sub_ps(_1,ve2),vzz));
			              __m256   vt11 = _mm256_mul_ps(_1,_mm256_add_ps(_1,vq));
			              __m256   vt12 = _mm256_mul_ps(vpp,_mm256_mul_ps(
			                                         _mm256_sub_ps(_1,ve2),vzz));
			              __m256   vt13 = _mm256_div_ps(vt12,_mm256_mul_ps(vq,
			                                         _mm256_add_ps(_1,vq)));
			              __m256   vt14 = _mm256_mul_ps(_0_5,_mm256_mul_ps(vpp,vrr));
			              __m256   vt15 = _mm256_sub_ps(vt9,_mm256_sub_ps(vt13,vt14))
			              __m256   vr0  = _mm256_add_ps(vt7,_mm256_sqrt_ps(_mm256_max_ps(_0,vt15)));
			                               vt15 = _mm256_sub_ps(vr,_mm256_mul_ps(ve2,vr0));
			              const __m256 vt16 = _mm256_mul_ps(vt15,vt15);
			              __m256   vu = _mm256_sqrt_ps(_mm256_add_ps(vt16,vzz));
			              __m256   vv = _mm256_sqrt_ps(_mm256_add_ps(vt16,
			                                       _mm256_mul_ps(_mm256_sub_ps(_1,ve2),vzz)));
			                        vt14 = _mm256_mul_ps(_mm256_mul_ps(vbb,vbb),vz);
			                        vt13 = _mm256_mul_ps(va,vv);
			              __m256   vz0  = _mm256_div_ps(vt14,vt13);
			                        vt12 = _mm256_sub_ps(_1,_mm256_div_ps(vbb,vt13));
				       __mm256_store_ps(&alt[i], _mm256_mul_ps(vu,vt12));
			               __m256_store_ps(&lat[i],   _mm256_atan2_ps(_mm256_fmadd_ps(
			                          _mm256_mul_ps(vep,vep),vz0,vz),vr));
			               __m256_store_ps(&lon[i], _mm256_atan2_ps(vy,vx));
			         
			            
				 }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {

				          float aa = a*a;
					  float bb = b*b;
					  float f  = (a-b)/a;
					  float e2 = (2.0F*f-f*f);
					  float ep = cephes_sqrtf(aa/bb-1.0F);
					  const float z  = pos_z[i];
					  float zz = z*z;
					  const float x  = pos_x[i];
					  const float y  = pos_y[i];
					  float r  = cephes_sqrtf(x*x+y*y);
					  float rr = r*r;
					  float ee = aa-bb;
					  float ff = 54.0*bb*zz;
					  float g  = rr+(1.0-e2)*zz-e2*ee;
					  float c  = e2*e2*ff*rr/(g*g*g);
					  float t0 = 1.0+c;
					  float s  = cephes_powf(t0+cephes_sqrtf(c*c+2.0*c),0.3333333333333333333333333333333F);
					  float t1 = (s+1.0F/s+1.0F)*(s+1.0F/s+1.0F);
					  float pp = ff/(3.0F*t1*g*g);
					  float t2 = e2*e2*pp;
					  float q  = cephes_sqrtf(1.0F+2.0F*t2);
					  float t3 = -pp*e2*r/(1.0F+q);
					  float t4 = 0.5*aa*(1.0F+1.0F/q);
					  float t5 = pp*(1.0F-e2)*zz;
					  float t6 = q*(1.0F+q);
					  float t7 = 0.5*pp*rr;
					  float r0 = t3+cephes_sqrtf(std::max(0.0,t4-(t5/t6)-t7));
					  const float t8 = (r-e2*r0)*(r-e2*r0);
					  float u  = cephes_sqrtf(t0+zz);
					  float v  = cephes_sqrtf(t0+(1.0F-e2)*zz);
					  const float t9 = a*v;
					  float z0 = b*b*z/t9;
					  lat[i] = cephes_atan2f((z+ep*ep*z0),r);
					  alt[i] = u*(1.0F-b2/t9);
					  lon[i] = cephes_atan2f(y,x);
				    }
				 
			}


                    	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                static inline
			void geodetic_to_cart_ymm8r4(const __m256 a,
			                              const __m256 b,
						      const __m256 lat,
						      const __m256 lon,
						      const __m256 alt,
						      __m256 &pos_x,
						      __m256 &pos_y,
						      __m256 &pos_z) {
 
                              register const __m256 _0 = _mm256_setzero_ps();
			      register const __m256 _1 = _mm256_set1_ps(1.0F);
			      register __m256 t0   = _0;
			      register __m256 t1   = _0;
			      register __m256 t2   = _0;
			      register __m256 t3   = _0;
			      register __m256 t4   = _0;
			      register __m256 zmm0 = _0; //sin(lat)
			      register __m256 zmm1 = _0; //cos(lat)
			      register __m256 zmm2 = _0; //tan(lat)
			      register __m256 zmm3 = _0; //sin(lon)
			      register __m256 zmm4 = _0; //cos(lon)
			      register __m256 ve2  = _0;
			      register __m256 vom2 = _0;
			      register __m256 vd   = _0;
			      register __m256 vq   = _0;
			      register __m256 vad  = _0;
			      register __m256 vaa  = _0;
			      register __m256 vbb  = _0;
			      zmm0  = _mm256_sin_ps(lat);
			      vaa   = _mm256_mul_ps(a,a);
			      vbb   = _mm256_mul_ps(b,b);
			      zmm1  = _mm256_cos_ps(lat);
			      ve2   = _mm256_sub_ps(_1,_mm256_div_ps(vbb,vaa));
			      vom2  = _mm256_sub_ps(_1,ve2);
			      zmm2  = _mm256_tan_ps(lat);
			      t0    = _mm256_mul_ps(zmm2,zmm2);
			      vd    = _mm256_sqrt_ps(_mm256_fmadd_ps(vom2,t0,_1));
			      t1    = _mm256_mul_ps(ve2,_mm256_mul_ps(zmm0,zmm0));
			      vq    = _mm256_sqrt_ps(_mm256_sub_ps(_1,t1));
			      zmm3  = _mm256_sin_ps(lon);
			      vad   = _mm256_div_ps(a,vd);
			      zmm4  = _mm256_cos_ps(lon);
			      t2    = _mm256_mul_ps(zmm1,zmm4);
			      pos_x = _mm256_fmadd_ps(alt,t2,_mm256_mul_ps(vad,zmm4));
			      t3    = _mm256_mul_ps(zmm3,zmm1);
			      pos_y = _mm256_fmadd_ps(alt,t3,_mm256_mul_ps(vad,zmm3));
			      t4    = _mm256_div_ps(zmm0,vq);
			      t3    = _mm256_mul_ps(a,_mm256_mul_ps(vom2,t4));
			      pos_z = _mm256_fmadd_ps(alt,zmm0,t3);
			      
			}


                    	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			geodetic_to_cart_u_ymm8r4_looped(const float a,
			                                  const float b,
							  const float * __restrict lat,
							  const float * __restrict lon,
							  const float * __restrict alt,
							  float * __restrict pos_x,
							  float * __restrict pos_y,
							  float * __restrict pos_z,
							  const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			         register const __m256 _0 = _mm256_setzero_ps();
			         register const __m256 _1 = _mm256_set1_ps(1.0F);
				 register const __m256 va   = _mm256_set1_ps(a);
				 register const __m256 vb   = _mm256_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
			            _mm_prefetch((const char*)&lat[i+32],_MM_HINT_T0);
			            register const __m256 vlat = _mm256_loadu_ps(&lat[i]);
                                    register __m256 zmm0  = _mm256_sin_ps(vlat);
			            register __m256 vaa   = _mm256_mul_ps(va,va);
			            register __m256 vbb   = _mm256_mul_ps(vb,vb);
				    register __m256 zmm1  = _mm256_cos_ps(vlat);
			            register __m256 ve2   = _mm256_sub_ps(_1,_mm256_div_ps(vbb,vaa));
			            register __m256 vom2  = _mm256_sub_ps(_1,ve2);
			            register __m256 zmm2  = _mm256_tan_ps(vlat);
			            register __m256 t0    = _mm256_mul_ps(zmm2,zmm2);
			            register __m256 vd    = _mm256_sqrt_ps(_mm256_fmadd_ps(vom2,t0,_1));
			            register __m256 t1    = _mm256_mul_ps(ve2,_mm256_mul_ps(zmm0,zmm0));
			            register __m256 vq    = _mm256_sqrt_ps(_mm256_sub_ps(_1,t1));
			            _mm_prefetch((const char*)&lon[i+32],_MM_HINT_T0);
				    register const __m256 vlon = _mm256_loadu_ps(&lon[i]);
			            register __m256 zmm3  = _mm256_sin_ps(vlon);
			            register __m256 vad   = _mm256_div_ps(va,vd);
			            register __m256 zmm4  = _mm256_cos_ps(vlon);
			            register __m256 t2    = _mm256_mul_ps(zmm1,zmm4);
				    _mm_prefetch((const char*)&alt[i+32],_MM_HINT_T0);
				    register const __m256 valt = _mm256_loadu_ps(&alt[i]);
				    _mm256_storeu_ps(&pos_x[i],
				                     _mm256_fmadd_ps(valt,t2,_mm256_mul_ps(vad,zmm4)));
			            t3    = _mm256_mul_ps(zmm3,zmm1);
			            _mm256_storeu_ps(&pos_y[i],
				                     _mm256_fmadd_ps(valt,t3,_mm256_mul_ps(vad,zmm3)));
			            t4    = _mm256_div_ps(zmm0,vq);
			            t3    = _mm256_mul_ps(va,_mm256_mul_ps(vom2,t4));
			            _mm256_storeu_ps(&pos_z[i],
				                     _mm256_fmadd_ps(valt,zmm0,t3));
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
				       
                                        register const  float glat = lat[i];
					register const  float s0 = cephes_sinf(glat);
					register const  float glon = lon[i];
					register const  float s1 = cephes_cosf(glat);
					register float ee  = 1.0F-(b*b)/(a*a);
					register const float s3 = cephes_sinf(glon);
					register float om2 = 1.0F-ee;
					register float q  = cephes_sqrtf(1.0-ee*s0*s0);
					register const float salt = alt[i];
					register const float s4 = cephes_cosf(glon);
					register const float s2 = cephes_tanf(glat);
					register float d = cephes_sqrtf(1.0F-ee*s2*s2);
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
			geodetic_to_cart_a_ymm8r4_looped(const float a,
			                                  const float b,
							  float * __restrict __ATTR_ALIGN__(32) lat,
							  float * __restrict __ATTR_ALIGN__(32) lon,
							  float * __restrict __ATTR_ALIGN__(32) alt,
							  float * __restrict __ATTR_ALIGN__(32) pos_x,
							  float * __restrict __ATTR_ALIGN__(32) pos_y,
							  float * __restrict __ATTR_ALIGN__(32) pos_z,
							  const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			         register const __m256 _0 = _mm256_setzero_ps();
			         register const __m256 _1 = _mm256_set1_ps(1.0F);
				 register const __m256 va   = _mm256_set1_ps(a);
				 register const __m256 vb   = _mm256_set1_ps(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                                 __assume_aligned(pos_x,32);
				 __assume_aligned(pos_y,32);
				 __assume_aligned(pos_z,32);
				 __assume_aligned(alt,32);
				 __assume_aligned(lon,32);
				 __assume_aligned(lat,32);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                                 pos_x = (float*)__builtin_assume_aligned(pos_x,32);
				 pos_y = (float*)__builtin_assume_aligned(pos_y,32);
				 pos_z = (float*)__builtin_assume_aligned(pos_z,32);
				 alt   = (float*)__builtin_assume_aligned(alt,32);
				 lon   = (float*)__builtin_assume_aligned(lon,32);
				 lat   = (float*)__builtin_assume_aligned(lat,32);
#endif				 
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
			            _mm_prefetch((const char*)&lat[i+32],_MM_HINT_T0);
			            register const __m256 vlat = _mm256_load_ps(&lat[i]);
                                    register __m256 zmm0  = _mm256_sin_ps(vlat);
			            register __m256 vaa   = _mm256_mul_ps(va,va);
			            register __m256 vbb   = _mm256_mul_ps(vb,vb);
				    register __m256 zmm1  = _mm256_cos_ps(vlat);
			            register __m256 ve2   = _mm256_sub_ps(_1,_mm256_div_ps(vbb,vaa));
			            register __m256 vom2  = _mm256_sub_ps(_1,ve2);
			            register __m256 zmm2  = _mm256_tan_ps(vlat);
			            register __m256 t0    = _mm256_mul_ps(zmm2,zmm2);
			            register __m256 vd    = _mm256_sqrt_ps(_mm256_fmadd_ps(vom2,t0,_1));
			            register __m256 t1    = _mm256_mul_ps(ve2,_mm256_mul_ps(zmm0,zmm0));
			            register __m256 vq    = _mm256_sqrt_ps(_mm256_sub_ps(_1,t1));
			            _mm_prefetch((const char*)&lon[i+32],_MM_HINT_T0);
				    register const __m256 vlon = _mm256_load_ps(&lon[i]);
			            register __m256 zmm3  = _mm256_sin_ps(vlon);
			            register __m256 vad   = _mm256_div_ps(va,vd);
			            register __m256 zmm4  = _mm256_cos_ps(vlon);
			            register __m256 t2    = _mm256_mul_ps(zmm1,zmm4);
				    _mm_prefetch((const char*)&alt[i+16],_MM_HINT_T0);
				    register const __m256 valt = _mm256_load_ps(&alt[i]);
				    _mm256_store_ps(&pos_x[i],
				                     _mm256_fmadd_ps(valt,t2,_mm256_mul_ps(vad,zmm4)));
			            t3    = _mm256_mul_ps(zmm3,zmm1);
			            _mm256_store_ps(&pos_y[i],
				                     _mm256_fmadd_ps(valt,t3,_mm256_mul_ps(vad,zmm3)));
			            t4    = _mm256_div_ps(zmm0,vq);
			            t3    = _mm256_mul_ps(va,_mm256_mul_ps(vom2,t4));
			            _mm256_store_ps(&pos_z[i],
				                     _mm256_fmadd_ps(valt,zmm0,t3));
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
				        register const  float glat = lat[i];
					register const  float s0 = cephes_sinf(glat);
					register const  float glon = lon[i];
					register const  float s1 = cephes_cosf(glat);
					register float ee  = 1.0F-(b*b)/(a*a);
					register const float s3 = cephes_sinf(glon);
					register float om2 = 1.0F-ee;
					register float q  = cephes_sqrtf(1.0-ee*s0*s0);
					register const float salt = alt[i];
					register const float s4 = cephes_cosf(glon);
					register const float s2 = cephes_tanf(glat);
					register float d = cephes_sqrtf(1.0F-ee*s2*s2);
					register float ad = a/d;
					pos_x[i] = ad*s4+salt*s3*s1;
					pos_y[i] = ad*s3+salt*s4*1;
					pos_z[i] = a*om2*s0/q+salt*s0;
                                       
				    }
		        }


			


						     
     }

}














#endif /*__GMS_GEODESY_YMM8R4_HPP__*/
