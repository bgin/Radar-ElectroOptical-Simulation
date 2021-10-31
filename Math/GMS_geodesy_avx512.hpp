
#ifndef __GMS_GEODESY_AVX512_HPP__
#define __GMS_GEODESY_AVX512_HPP__ 161020210952


namespace file_version {

    const unsigned int gGMS_GEODESY_AVX512_MAJOR = 1U;
    const unsigned int gGMS_GEODESY_AVX512_MINOR = 0U;
    const unsigned int gGMS_GEODESY_AVX512_MICRO = 0U;
    const unsigned int gGMS_GEODESY_AVX512_FULLVER =
      1000U*gGMS_GEODESY_AVX512_MAJOR+
      100U*gGMS_GEODESY_AVX512_MINOR+
      10U*gGMS_GEODESY_AVX512_MICRO;
    const char * const pgGMS_GEODESY_AVX512_CREATION_DATE = "16-10-2021 09:52 AM +00200 (SAT 16 OCT 2021 GMT+2)";
    const char * const pgGMS_GEODESY_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_GEODESY_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_GEODESY_AVX512_DESCRIPTION   = "Vectorized (AVX512) geodesic computation implementation."

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
			cart_to_geodetic_zmm8r8(     const __m512d pos_x, //input position x [km]
			                             const __m512d pos_y, //input position y [km]
						     const __m512d pos_z, //input position z [km]
						     const __m512d a, // semi-minor axis [km]
						     const __m512d b, // semi-major axis [km]
						     __m512d &alt, //output altitude [km]
						     __m512d &lon, //output longtitude [rad]
						     __m512d &lat) { // output latitiude [km]
                            const __m512d _0   = _mm512_setzero_pd();
			    const __m512d _1   = _mm512_set1_pd(1.0);
			    const __m512d _2   = _mm512_set1_pd(1.0);
			    const __m512d _54  = _mm512_set1_pd(54.0);
			    const __m512d _3   = _mm512_set1_pd(3.0);
			    const __m512d _0_5 = _mm512_set1_pd(0.5);
			    const __m512d _0_3 = _mm512_set1_pd(0.3333333333333333333333333333333);
			    __m512d vaa = _0;
			    __m512d vbb = _0;
			    __m512d vf  = _0;
			    __m512d ve2 = _0;
			    __m512d vep = _0;
			    __m512d vzz = _0;
			    __m512d vr  = _0;
                            __m512d vrr = _0;
			    __m512d vee = _0;
			    __m512d vff = _0;
			    __m512d vg  = _0;
			    __m512d vc  = _0;
			    __m512d vs  = _0;
			    __m512d vpp = _0;
			   
			    __m512d vq  = _0;
			    __m512d vr0 = _0;
			    __m512d vt0 = _0;
			    __m512d vt1 = _0;
			    __m512d vt2 = _0;
			    __m512d vt3 = _0;
			    __m512d vt4 = _0;
			    __m512d vt5 = _0;
			    __m512d vt6 = _0;
			    __m512d vt7 = _0;
			    __m512d vt8 = _0;
			    __m512d vt9 = _0;
			    __m512d vt10 = _0;
			    __m512d vt11 = _0;
			    __m512d vt12 = _0;
			    __m512d vt13 = _0;
			    __m512d vt14 = _0;
			    __m512d vt15 = _0;
			    __m512d vu  = _0;
			    __m512d vv  = _0;
			    __m512d vz0 = _0;
			   			    
			    vf  = _mm512_div_pd(_mm512_sub_pd(a,b),a);
			    vaa = _mm512_mul_pd(a,a);
			    ve2 = _mm512_sub_pd(_mm512_mul_pd(_2,vf),
			                        _mm512_mul_pd(vf,vf));
			    vbb = _mm512_mul_pd(b,b);
			    vzz = _mm512_mul_pd(pos_z,pos_z);
			    ve2 = _mm512_sub_pd(aa,bb);
			    vr  = _mm512_sqrt_pd(_mm512_mul_pd(pos_x,pos_x),
			                         _mm512_mul_pd(pos_y,pos_y));
			    vrr = _mm512_mul_pd(vr,vr);
			    vff = _mm512_mul_pd(_54,
			                        _mm512_mul_pd(vbb,vzz));
			    vep = _mm512_sqrt_pd(_mm512_sub_pd(
			                         _mm512_div_pd(vaa,vbb),_1));
			    vt0 = _mm512_add_pd(vrr,_mm512_mul_pd(
			                        _mm512_sub_pd(_1,ve2),vzz));
			    vg  = _mm512_sub_pd(vt0,_mm512_mul_pd(ve2,vee));
			    const __m512d vg2 = _mm512_mul_pd(vg,vg);
			    const __m512d vg3  = _mm512_mul_pd(vg2,vg);
			    const __m512d ve22 = _mm512_mul_pd(ve2,ve2);
			    vc = _mm512_div_pd(_mm512_mul_pd(ve22,
			                       _mm512_mul_pd(vff,vrr)),vg3);
			    vt1 = _mm512_add_pd(_1,vc);
			    vt2 = _mm512_mul_pd(vc,vc);
			    vs  = _mm512_pow(_mm512_add_pd(vt1,
			                     _mm512_sqrt_pd(
					     _mm512_fmadd_pd(_2,vc,vt2))),_0_3);
			    vt3 = _mm512_add_pd(vs,_mm512_add_pd(
			                           _mm512_div_pd(_1,vs),_1));
                            vt4 = _mm512_mul_pd(_mm512_mul_pd(_3,
			                        _mm512_mul_pd(vt3,vt3)),vg2);
			    vpp = _mm512_div_pd(ff,vt4);
			    vt5 = _mm512_mul_pd(ve22,vpp);
			    vq  = _mm512_sqrt_pd(_mm512_fmadd_pd(_2,vt5,_1));
			    vt6 = _mm512_sub_pd(_0,vpp);
			    vt7 = _mm512_div_pd(_mm512_mul_pd(vt6,
			                        _mm512_mul_pd(ve2,vr)),
						_mm512_add_pd(_1,vq));
			    vt8 = _mm512_mul_pd(_0_5,vaa);
			    vt9 = _mm512_mul_pd(vt8,_mm512_add_pd(_1,
			                                          _mm512_div_pd(_1,vq)));
			    vt10 = _mm512_mul_pd(vpp,_mm512_mul_pd(
			                                      _mm512_sub_pd(_1,ve2),vzz));
			    vt11 = _mm512_mul_pd(_1,_mm512_add_pd(_1,vq));
			    vt12 = _mm512_mul_pd(vpp,_mm512_mul_pd(
			                             _mm512_sub_pd(_1,ve2),vzz));
			    vt13 = _mm512_div_pd(vt12,_mm512_mul_pd(vq,
			                              _mm512_add_pd(_1,vq)));
			    vt14 = _mm512_mul_pd(_0_5,_mm512_mul_pd(vpp,vrr));
			    vt15 = _mm512_sub_pd(vt9,_mm512_sub_pd(vt13,vt14))
			    vr0  = _mm512_add_pd(vt7,_mm512_sqrt_pd(_mm512_max_pd(_0,vt15)));
			    vt15 = _mm512_sub_pd(vr,_mm512_mul_pd(ve2,vr0));
			    const __m512d vt16 = _mm512_mul_pd(vt15,vt15);
			    vu = _mm512_sqrt_pd(_mm512_add_pd(vt16,vzz));
			    vv = _mm512_sqrt_pd(_mm512_add_pd(vt16,
			                        _mm512_mul_pd(_mm512_sub_pd(_1,ve2),vzz)));
			    vt14 = _mm512_mul_pd(_mm512_mul_pd(vbb,vbb),vz);
			    vt13 = _mm512_mul_pd(a,vv);
			    vz0  = _mm512_div_pd(vt14,vt13);
			    vt12 = _mm512_sub_pd(_1,_mm512_div_pd(vbb,vt13));
			    alt = _mm512_mul_pd(vu,vt12);
			    lat = _mm512_atan2_pd(_mm512_fmadd_pd(
			                          _mm512_mul_pd(vep,vep),vz0,pos_z),vr);
			    lon = _mm512_atan2_pd(pos_y,pos_x);
			}


			
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			cart_to_geodetic_u_zmm8r8_looped(const double * __restrict pos_x,
			                                 const double * __restrict pos_y,
							 const double * __restrict pos_z,
							 const double a,
							 const double b,
							 double * __restrict alt,
							 double * __restrict lon,
							 double * __restrict lat,
							 const int32_t n) {

			      if(__builtin_expect(n<=0,0)) {return;}
			         const __m512d _0   = _mm512_setzero_pd();
			         const __m512d _1   = _mm512_set1_pd(1.0);
			         const __m512d _2   = _mm512_set1_pd(1.0);
			         const __m512d _54  = _mm512_set1_pd(54.0);
			         const __m512d _3   = _mm512_set1_pd(3.0);
			         const __m512d _0_5 = _mm512_set1_pd(0.5);
			         const __m512d _0_3 = _mm512_set1_pd(0.3333333333333333333333333333333);
				 const __m512d va   = _mm512_set1_pd(a);
				 const __m512d vb   = _mm512_set1_pd(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			         for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
				 
                                      __m512d  vf  = _mm512_div_pd(_mm512_sub_pd(va,vb),va);
			              __m512d  vaa = _mm512_mul_pd(va,va);
			              __m512d  ve2 = _mm512_sub_pd(_mm512_mul_pd(_2,vf),
			                                           _mm512_mul_pd(vf,vf));
			              __m512d  vbb = _mm512_mul_pd(vb,vb);
				       _mm_prefetch((const char*)&pos_z[i+8],_MM_HINT_T0);
				      const __m512d vz = _mm512_loadu_pd(&pos_z[i]);
			              __m512d  vzz = _mm512_mul_pd(vz,vz);
			              __m512d  ve2 = _mm512_sub_pd(aa,bb);
				      _mm_prefetch((const char*)&pos_x[i+8],_MM_HINT_T0);
				      _mm_prefetch((const char*)&pos_y[i+8],_MM_HINT_T0);
				      const __m512d vx = _mm512_loadu_pd(&pos_x[i]);
				      const __m512d vy = _mm512_loadu_pd(&pos_y[i]);
			              __m512d vr  = _mm512_sqrt_pd(_mm512_mul_pd(vx,vx),
			                                            _mm512_mul_pd(vy,vy));
			              __m512d  vrr = _mm512_mul_pd(vr,vr);
			              __m512d  vff = _mm512_mul_pd(_54,
			                                   _mm512_mul_pd(vbb,vzz));
			              __m512d  vep = _mm512_sqrt_pd(_mm512_sub_pd(
			                                   _mm512_div_pd(vaa,vbb),_1));
			              __m512d  vt0 = _mm512_add_pd(vrr,_mm512_mul_pd(
			                                    _mm512_sub_pd(_1,ve2),vzz));
			              __m512d   vg  = _mm512_sub_pd(vt0,_mm512_mul_pd(ve2,vee));
			              const __m512d vg2 = _mm512_mul_pd(vg,vg);
			              const __m512d vg3  = _mm512_mul_pd(vg2,vg);
			              const __m512d ve22 = _mm512_mul_pd(ve2,ve2);
			              __m512d   vc = _mm512_div_pd(_mm512_mul_pd(ve22,
			                                     _mm512_mul_pd(vff,vrr)),vg3);
			              __m512d   vt1 = _mm512_add_pd(_1,vc);
			              __m512d   vt2 = _mm512_mul_pd(vc,vc);
			              __m512d   vs  = _mm512_pow(_mm512_add_pd(vt1,
			                                        _mm512_sqrt_pd(
					                           _mm512_fmadd_pd(_2,vc,vt2))),_0_3);
			              __m512d   vt3 = _mm512_add_pd(vs,_mm512_add_pd(
			                                         _mm512_div_pd(_1,vs),_1));
                                      __m512d   vt4 = _mm512_mul_pd(_mm512_mul_pd(_3,
			                                    _mm512_mul_pd(vt3,vt3)),vg2);
			              __m512d   vpp = _mm512_div_pd(ff,vt4);
			              __m512d   vt5 = _mm512_mul_pd(ve22,vpp);
			              __m512d   vq  = _mm512_sqrt_pd(_mm512_fmadd_pd(_2,vt5,_1));
			              __m512d   vt6 = _mm512_sub_pd(_0,vpp);
			              __m512d   vt7 = _mm512_div_pd(_mm512_mul_pd(vt6,
			                                              _mm512_mul_pd(ve2,vr)),
						                      _mm512_add_pd(_1,vq));
			              __m512d   vt8 = _mm512_mul_pd(_0_5,vaa);
			              __m512d   vt9 = _mm512_mul_pd(vt8,_mm512_add_pd(_1,
			                                          _mm512_div_pd(_1,vq)));
			              __m512d   vt10 = _mm512_mul_pd(vpp,_mm512_mul_pd(
			                                          _mm512_sub_pd(_1,ve2),vzz));
			              __m512d   vt11 = _mm512_mul_pd(_1,_mm512_add_pd(_1,vq));
			              __m512d   vt12 = _mm512_mul_pd(vpp,_mm512_mul_pd(
			                                         _mm512_sub_pd(_1,ve2),vzz));
			              __m512d   vt13 = _mm512_div_pd(vt12,_mm512_mul_pd(vq,
			                                         _mm512_add_pd(_1,vq)));
			              __m512d   vt14 = _mm512_mul_pd(_0_5,_mm512_mul_pd(vpp,vrr));
			              __m512d   vt15 = _mm512_sub_pd(vt9,_mm512_sub_pd(vt13,vt14))
			              __m512d   vr0  = _mm512_add_pd(vt7,_mm512_sqrt_pd(_mm512_max_pd(_0,vt15)));
			                               vt15 = _mm512_sub_pd(vr,_mm512_mul_pd(ve2,vr0));
			              const __m512d vt16 = _mm512_mul_pd(vt15,vt15);
			              __m512d   vu = _mm512_sqrt_pd(_mm512_add_pd(vt16,vzz));
			              __m512d   vv = _mm512_sqrt_pd(_mm512_add_pd(vt16,
			                                       _mm512_mul_pd(_mm512_sub_pd(_1,ve2),vzz)));
			                        vt14 = _mm512_mul_pd(_mm512_mul_pd(vbb,vbb),vz);
			                        vt13 = _mm512_mul_pd(va,vv);
			              __m512d   vz0  = _mm512_div_pd(vt14,vt13);
			                        vt12 = _mm512_sub_pd(_1,_mm512_div_pd(vbb,vt13));
				       __mm512_storeu_pd(&alt[i], _mm512_mul_pd(vu,vt12));
			               __m512_storeu_pd(&lat[i],   _mm512_atan2_pd(_mm512_fmadd_pd(
			                          _mm512_mul_pd(vep,vep),vz0,vz),vr));
			               __m512_storeu_pd(&lon[i], _mm512_atan2_pd(vy,vx));
				 }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {

				          double aa = a*a;
					  double bb = b*b;
					  double f  = (a-b)/a;
					  double e2 = (2.0*f-f*f);
					  double ep = std::sqrt(aa/bb-1.0);
					  const double z  = pos_z[i];
					  double zz = z*z;
					  const double x  = pos_x[i];
					  const double y  = pos_y[i];
					  double r  = std::sqrt(x*x+y*y);
					  double rr = r*r;
					  double ee = aa-bb;
					  double ff = 54.0*bb*zz;
					  double g  = rr+(1.0-e2)*zz-e2*ee;
					  double c  = e2*e2*ff*rr/(g*g*g);
					  double t0 = 1.0+c;
					  double s  = std::pow(t0+std::sqrt(c*c+2.0*c),0.3333333333333333333333333333333);
					  double t1 = (s+1.0/s+1.0)*(s+1.0/s+1.0);
					  double pp = ff/(3.0*t1*g*g);
					  double t2 = e2*e2*pp;
					  double q  = std::sqrt(1.0+2.0*t2);
					  double t3 = -pp*e2*r/(1.0+q);
					  double t4 = 0.5*aa*(1.0+1.0/q);
					  double t5 = pp*(1.0-e2)*zz;
					  double t6 = q*(1.0+q);
					  double t7 = 0.5*pp*rr;
					  double r0 = t3+std::sqrt(std::max(0.0,t4-(t5/t6)-t7));
					  const double t8 = (r-e2*r0)*(r-e2*r0);
					  double u  = std::sqrt(t0+zz);
					  double v  = std::sqrt(t0+(1.0-e2)*zz);
					  const double t9 = a*v;
					  double z0 = b*b*z/t9;
					  lat[i] = std::atan2((z+ep*ep*z0),r);
					  alt[i] = u*(1.0-b2/t9);
					  lon[i] = std::atan2(y,x);
				    }
				 
			}



			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			cart_to_geodetic_a_zmm8r8_looped(double * __restrict __ATTR_ALIGN__(64) pos_x,
			                                 double * __restrict __ATTR_ALIGN__(64) pos_y,
							 double * __restrict __ATTR_ALIGN__(64) pos_z,
							 const double a,
							 const double b,
							 double * __restrict __ATTR_ALIGN__(64) alt,
							 double * __restrict __ATTR_ALIGN__(64) lon,
							 double * __restrict __ATTR_ALIGN__(64) lat,
							 const int32_t n) {

			      if(__builtin_expect(n<=0,0)) {return;}
			         const __m512d _0   = _mm512_setzero_pd();
			         const __m512d _1   = _mm512_set1_pd(1.0);
			         const __m512d _2   = _mm512_set1_pd(1.0);
			         const __m512d _54  = _mm512_set1_pd(54.0);
			         const __m512d _3   = _mm512_set1_pd(3.0);
			         const __m512d _0_5 = _mm512_set1_pd(0.5);
			         const __m512d _0_3 = _mm512_set1_pd(0.3333333333333333333333333333333);
				 const __m512d va   = _mm512_set1_pd(a);
				 const __m512d vb   = _mm512_set1_pd(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                                 __assume_aligned(pos_x,64);
				 __assume_aligned(pos_y,64);
				 __assume_aligned(pos_z,64);
				 __assume_aligned(alt,64);
				 __assume_aligned(lon,64);
				 __assume_aligned(lat,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                                 pos_x = (double*)__builtin_assume_aligned(pos_x,64);
				 pos_y = (double*)__builtin_assume_aligned(pos_y,64);
				 pos_z = (double*)__builtin_assume_aligned(pos_z,64);
				 alt   = (double*)__builtin_assume_aligned(alt,64);
				 lon   = (double*)__builtin_assume_aligned(lon,64);
				 lat   = (double*)__builtin_assume_aligned(lat,64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			         for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
				 
                                      __m512d  vf  = _mm512_div_pd(_mm512_sub_pd(va,vb),va);
			              __m512d  vaa = _mm512_mul_pd(va,va);
			              __m512d  ve2 = _mm512_sub_pd(_mm512_mul_pd(_2,vf),
			                                           _mm512_mul_pd(vf,vf));
			              __m512d  vbb = _mm512_mul_pd(vb,vb);
				       _mm_prefetch((const char*)&pos_z[i+8],_MM_HINT_T0);
				      const __m512d vz = _mm512_loadu_pd(&pos_z[i]);
			              __m512d  vzz = _mm512_mul_pd(vz,vz);
			              __m512d  ve2 = _mm512_sub_pd(aa,bb);
				      _mm_prefetch((const char*)&pos_x[i+8],_MM_HINT_T0);
				      _mm_prefetch((const char*)&pos_y[i+8],_MM_HINT_T0);
				      const __m512d vx = _mm512_loadu_pd(&pos_x[i]);
				      const __m512d vy = _mm512_loadu_pd(&pos_y[i]);
			              __m512d vr  = _mm512_sqrt_pd(_mm512_mul_pd(vx,vx),
			                                            _mm512_mul_pd(vy,vy));
			              __m512d  vrr = _mm512_mul_pd(vr,vr);
			              __m512d  vff = _mm512_mul_pd(_54,
			                                   _mm512_mul_pd(vbb,vzz));
			              __m512d  vep = _mm512_sqrt_pd(_mm512_sub_pd(
			                                   _mm512_div_pd(vaa,vbb),_1));
			              __m512d  vt0 = _mm512_add_pd(vrr,_mm512_mul_pd(
			                                    _mm512_sub_pd(_1,ve2),vzz));
			              __m512d   vg  = _mm512_sub_pd(vt0,_mm512_mul_pd(ve2,vee));
			              const __m512d vg2 = _mm512_mul_pd(vg,vg);
			              const __m512d vg3  = _mm512_mul_pd(vg2,vg);
			              const __m512d ve22 = _mm512_mul_pd(ve2,ve2);
			              __m512d   vc = _mm512_div_pd(_mm512_mul_pd(ve22,
			                                     _mm512_mul_pd(vff,vrr)),vg3);
			              __m512d   vt1 = _mm512_add_pd(_1,vc);
			              __m512d   vt2 = _mm512_mul_pd(vc,vc);
			              __m512d   vs  = _mm512_pow(_mm512_add_pd(vt1,
			                                        _mm512_sqrt_pd(
					                           _mm512_fmadd_pd(_2,vc,vt2))),_0_3);
			              __m512d   vt3 = _mm512_add_pd(vs,_mm512_add_pd(
			                                         _mm512_div_pd(_1,vs),_1));
                                      __m512d   vt4 = _mm512_mul_pd(_mm512_mul_pd(_3,
			                                    _mm512_mul_pd(vt3,vt3)),vg2);
			              __m512d   vpp = _mm512_div_pd(ff,vt4);
			              __m512d   vt5 = _mm512_mul_pd(ve22,vpp);
			              __m512d   vq  = _mm512_sqrt_pd(_mm512_fmadd_pd(_2,vt5,_1));
			              __m512d   vt6 = _mm512_sub_pd(_0,vpp);
			              __m512d   vt7 = _mm512_div_pd(_mm512_mul_pd(vt6,
			                                              _mm512_mul_pd(ve2,vr)),
						                      _mm512_add_pd(_1,vq));
			              __m512d   vt8 = _mm512_mul_pd(_0_5,vaa);
			              __m512d   vt9 = _mm512_mul_pd(vt8,_mm512_add_pd(_1,
			                                          _mm512_div_pd(_1,vq)));
			              __m512d   vt10 = _mm512_mul_pd(vpp,_mm512_mul_pd(
			                                          _mm512_sub_pd(_1,ve2),vzz));
			              __m512d   vt11 = _mm512_mul_pd(_1,_mm512_add_pd(_1,vq));
			              __m512d   vt12 = _mm512_mul_pd(vpp,_mm512_mul_pd(
			                                         _mm512_sub_pd(_1,ve2),vzz));
			              __m512d   vt13 = _mm512_div_pd(vt12,_mm512_mul_pd(vq,
			                                         _mm512_add_pd(_1,vq)));
			              __m512d   vt14 = _mm512_mul_pd(_0_5,_mm512_mul_pd(vpp,vrr));
			              __m512d   vt15 = _mm512_sub_pd(vt9,_mm512_sub_pd(vt13,vt14))
			              __m512d   vr0  = _mm512_add_pd(vt7,_mm512_sqrt_pd(_mm512_max_pd(_0,vt15)));
			                               vt15 = _mm512_sub_pd(vr,_mm512_mul_pd(ve2,vr0));
			              const __m512d vt16 = _mm512_mul_pd(vt15,vt15);
			              __m512d   vu = _mm512_sqrt_pd(_mm512_add_pd(vt16,vzz));
			              __m512d   vv = _mm512_sqrt_pd(_mm512_add_pd(vt16,
			                                       _mm512_mul_pd(_mm512_sub_pd(_1,ve2),vzz)));
			                        vt14 = _mm512_mul_pd(_mm512_mul_pd(vbb,vbb),vz);
			                        vt13 = _mm512_mul_pd(va,vv);
			              __m512d   vz0  = _mm512_div_pd(vt14,vt13);
			                        vt12 = _mm512_sub_pd(_1,_mm512_div_pd(vbb,vt13));
				       __mm512_storeu_pd(&alt[i], _mm512_mul_pd(vu,vt12));
			               __m512_storeu_pd(&lat[i],   _mm512_atan2_pd(_mm512_fmadd_pd(
			                          _mm512_mul_pd(vep,vep),vz0,vz),vr));
			               __m512_storeu_pd(&lon[i], _mm512_atan2_pd(vy,vx));
				 }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {

				          double aa = a*a;
					  double bb = b*b;
					  double f  = (a-b)/a;
					  double e2 = (2.0*f-f*f);
					  double ep = std::sqrt(aa/bb-1.0);
					  const double z  = pos_z[i];
					  double zz = z*z;
					  const double x  = pos_x[i];
					  const double y  = pos_y[i];
					  double r  = std::sqrt(x*x+y*y);
					  double rr = r*r;
					  double ee = aa-bb;
					  double ff = 54.0*bb*zz;
					  double g  = rr+(1.0-e2)*zz-e2*ee;
					  double c  = e2*e2*ff*rr/(g*g*g);
					  double t0 = 1.0+c;
					  double s  = std::pow(t0+std::sqrt(c*c+2.0*c),0.3333333333333333333333333333333);
					  double t1 = (s+1.0/s+1.0)*(s+1.0/s+1.0);
					  double pp = ff/(3.0*t1*g*g);
					  double t2 = e2*e2*pp;
					  double q  = std::sqrt(1.0+2.0*t2);
					  double t3 = -pp*e2*r/(1.0+q);
					  double t4 = 0.5*aa*(1.0+1.0/q);
					  double t5 = pp*(1.0-e2)*zz;
					  double t6 = q*(1.0+q);
					  double t7 = 0.5*pp*rr;
					  double r0 = t3+std::sqrt(std::max(0.0,t4-(t5/t6)-t7));
					  const double t8 = (r-e2*r0)*(r-e2*r0);
					  double u  = std::sqrt(t0+zz);
					  double v  = std::sqrt(t0+(1.0-e2)*zz);
					  const double t9 = a*v;
					  double z0 = b*b*z/t9;
					  lat[i] = std::atan2((z+ep*ep*z0),r);
					  alt[i] = u*(1.0-b2/t9);
					  lon[i] = std::atan2(y,x);
				    }
				 
			}


                    	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                static inline
			void geodetic_to_cart_zmm8r8(const __m512d a,
			                             const __m512d b,
						     const __m512d lat,
						     const __m512d lon,
						     const __m512d alt,
						     __m512d &pos_x,
						     __m512d &pos_y,
						     __m512d &pos_z) {
 
                              register const __m512d _0 = _mm512_setzero_pd();
			      register const __m512d _1 = _mm512_set1_pd(1.0);
			      register __m512d t0   = _0;
			      register __m512d t1   = _0;
			      register __m512d t2   = _0;
			      register __m512d t3   = _0;
			      register __m512d t4   = _0;
			      register __m512d zmm0 = _0; //sin(lat)
			      register __m512d zmm1 = _0; //cos(lat)
			      register __m512d zmm2 = _0; //tan(lat)
			      register __m512d zmm3 = _0; //sin(lon)
			      register __m512d zmm4 = _0; //cos(lon)
			      register __m512d ve2  = _0;
			      register __m512d vom2 = _0;
			      register __m512d vd   = _0;
			      register __m512d vq   = _0;
			      register __m512d vad  = _0;
			      register __m512d vaa  = _0;
			      register __m512d vbb  = _0;
			      zmm0  = _mm512_sin_pd(lat);
			      vaa   = _mm512_mul_pd(a,a);
			      vbb   = _mm512_mul_pd(b,b);
			      zmm1  = _mm512_cos_pd(lat);
			      ve2   = _mm512_sub_pd(_1,_mm512_div_pd(vbb,vaa));
			      vom2  = _mm512_sub_pd(_1,ve2);
			      zmm2  = _mm512_tan_pd(lat);
			      t0    = _mm512_mul_pd(zmm2,zmm2);
			      vd    = _mm512_sqrt_pd(_mm512_fmadd_pd(vom2,t0,_1));
			      t1    = _mm512_mul_pd(ve2,_mm512_mul_pd(zmm0,zmm0));
			      vq    = _mm512_sqrt_pd(_mm512_sub_pd(_1,t1));
			      zmm3  = _mm512_sin_pd(lon);
			      vad   = _mm512_div_pd(a,vd);
			      zmm4  = _mm512_cos_pd(lon);
			      t2    = _mm512_mul_pd(zmm1,zmm4);
			      pos_x = _mm512_fmadd_pd(alt,t2,_mm512_mul_pd(vad,zmm4));
			      t3    = _mm512_mul_pd(zmm3,zmm1);
			      pos_y = _mm512_fmadd_pd(alt,t3,_mm512_mul_pd(vad,zmm3));
			      t4    = _mm512_div_pd(zmm0,vq);
			      t3    = _mm512_mul_pd(a,_mm512_mul_pd(vom2,t4));
			      pos_z = _mm512_fmadd_pd(alt,zmm0,t3);
			      
			}


                    	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
			void
			geodetic_to_cart_u_zmm8r8_looped(const double a,
			                                 const double b,
							 const double * __restrict lat,
							 const double * __restrict lon,
							 const double * __restrict alt,
							 double * __restrict pos_x,
							 double * __restrict pos_y,
							 double * __restrict pos_z,
							 const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			         register const __m512d _0 = _mm512_setzero_pd();
			         register const __m512d _1 = _mm512_set1_pd(1.0);
				 register const __m512d va   = _mm512_set1_pd(a);
				 register const __m512d vb   = _mm512_set1_pd(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
			            _mm_prefetch((const char*)&lat[i+8],_MM_HINT_T0);
			            register const __m512d vlat = _mm512_loadu_pd(&lat[i]);
                                    register __m512d zmm0  = _mm512_sin_pd(vlat);
			            register __m512d vaa   = _mm512_mul_pd(va,va);
			            register __m512d vbb   = _mm512_mul_pd(vb,vb);
				    register __m512d zmm1  = _mm512_cos_pd(vlat);
			            register __m512d ve2   = _mm512_sub_pd(_1,_mm512_div_pd(vbb,vaa));
			            register __m512d vom2  = _mm512_sub_pd(_1,ve2);
			            register __m512d zmm2  = _mm512_tan_pd(vlat);
			            register __m512d t0    = _mm512_mul_pd(zmm2,zmm2);
			            register __m512d vd    = _mm512_sqrt_pd(_mm512_fmadd_pd(vom2,t0,_1));
			            register __m512d t1    = _mm512_mul_pd(ve2,_mm512_mul_pd(zmm0,zmm0));
			            register __m512d vq    = _mm512_sqrt_pd(_mm512_sub_pd(_1,t1));
			            _mm_prefetch((const char*)&lon[i+8],_MM_HINT_T0);
				    register const __m512d vlon = _mm512_loadu_pd(&lon[i]);
			            register __m512d zmm3  = _mm512_sin_pd(vlon);
			            register __m512d vad   = _mm512_div_pd(va,vd);
			            register __m512d zmm4  = _mm512_cos_pd(vlon);
			            register __m512d t2    = _mm512_mul_pd(zmm1,zmm4);
				    _mm_prefetch((const char*)&alt[i+8],_MM_HINT_T0);
				    register const __m512d valt = _mm512_loadu_pd(&alt[i]);
				    _mm512_storeu_pd(&pos_x[i],
				                     _mm512_fmadd_pd(valt,t2,_mm512_mul_pd(vad,zmm4)));
			            t3    = _mm512_mul_pd(zmm3,zmm1);
			            _mm512_storeu_pd(&pos_y[i],
				                     _mm512_fmadd_pd(valt,t3,_mm512_mul_pd(vad,zmm3)));
			            t4    = _mm512_div_pd(zmm0,vq);
			            t3    = _mm512_mul_pd(va,_mm512_mul_pd(vom2,t4));
			            _mm512_storeu_pd(&pos_z[i],
				                     _mm512_fmadd_pd(valt,zmm0,t3));
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
				       
                                        register const  double glat = lat[i];
					register const double s0 = std::sin(glat);
					register const  double glon = lon[i];
					register const double s1 = std::cos(glat);
					register double ee  = 1.0-(b*b)/(a*a);
					register const double s3 = std::sin(glon);
					register double om2 = 1.0-ee;
					register double q  = std::sqrt(1.0-ee*s0*s0);
					register const double salt = alt[i];
					register const double s4 = std::cos(glon);
					register const double s2 = std::tan(glat);
					register double d = std::sqrt(1.0-ee*s2*s2);
					register double ad = a/d;
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
			geodetic_to_cart_a_zmm8r8_looped(const double a,
			                                 const double b,
							 double * __restrict __ATTR_ALIGN__(64) lat,
							 double * __restrict __ATTR_ALIGN__(64) lon,
							 double * __restrict __ATTR_ALIGN__(64) alt,
							 double * __restrict __ATTR_ALIGN__(64) pos_x,
							 double * __restrict __ATTR_ALIGN__(64) pos_y,
							 double * __restrict __ATTR_ALIGN__(64) pos_z,
							 const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			         register const __m512d _0 = _mm512_setzero_pd();
			         register const __m512d _1 = _mm512_set1_pd(1.0);
				 register const __m512d va   = _mm512_set1_pd(a);
				 register const __m512d vb   = _mm512_set1_pd(b);
				 int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                                 __assume_aligned(pos_x,64);
				 __assume_aligned(pos_y,64);
				 __assume_aligned(pos_z,64);
				 __assume_aligned(alt,64);
				 __assume_aligned(lon,64);
				 __assume_aligned(lat,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                                 pos_x = (double*)__builtin_assume_aligned(pos_x,64);
				 pos_y = (double*)__builtin_assume_aligned(pos_y,64);
				 pos_z = (double*)__builtin_assume_aligned(pos_z,64);
				 alt   = (double*)__builtin_assume_aligned(alt,64);
				 lon   = (double*)__builtin_assume_aligned(lon,64);
				 lat   = (double*)__builtin_assume_aligned(lat,64);
#endif				 
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
			            _mm_prefetch((const char*)&lat[i+8],_MM_HINT_T0);
			            register const __m512d vlat = _mm512_load_pd(&lat[i]);
                                    register __m512d zmm0  = _mm512_sin_pd(vlat);
			            register __m512d vaa   = _mm512_mul_pd(va,va);
			            register __m512d vbb   = _mm512_mul_pd(vb,vb);
				    register __m512d zmm1  = _mm512_cos_pd(vlat);
			            register __m512d ve2   = _mm512_sub_pd(_1,_mm512_div_pd(vbb,vaa));
			            register __m512d vom2  = _mm512_sub_pd(_1,ve2);
			            register __m512d zmm2  = _mm512_tan_pd(vlat);
			            register __m512d t0    = _mm512_mul_pd(zmm2,zmm2);
			            register __m512d vd    = _mm512_sqrt_pd(_mm512_fmadd_pd(vom2,t0,_1));
			            register __m512d t1    = _mm512_mul_pd(ve2,_mm512_mul_pd(zmm0,zmm0));
			            register __m512d vq    = _mm512_sqrt_pd(_mm512_sub_pd(_1,t1));
			            _mm_prefetch((const char*)&lon[i+8],_MM_HINT_T0);
				    register const __m512d vlon = _mm512_load_pd(&lon[i]);
			            register __m512d zmm3  = _mm512_sin_pd(vlon);
			            register __m512d vad   = _mm512_div_pd(va,vd);
			            register __m512d zmm4  = _mm512_cos_pd(vlon);
			            register __m512d t2    = _mm512_mul_pd(zmm1,zmm4);
				    _mm_prefetch((const char*)&alt[i+8],_MM_HINT_T0);
				    register const __m512d valt = _mm512_load_pd(&alt[i]);
				    _mm512_store_pd(&pos_x[i],
				                     _mm512_fmadd_pd(valt,t2,_mm512_mul_pd(vad,zmm4)));
			            t3    = _mm512_mul_pd(zmm3,zmm1);
			            _mm512_store_pd(&pos_y[i],
				                     _mm512_fmadd_pd(valt,t3,_mm512_mul_pd(vad,zmm3)));
			            t4    = _mm512_div_pd(zmm0,vq);
			            t3    = _mm512_mul_pd(va,_mm512_mul_pd(vom2,t4));
			            _mm512_store_pd(&pos_z[i],
				                     _mm512_fmadd_pd(valt,zmm0,t3));
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
				       
                                        register const  double glat = lat[i];
					register const double s0 = std::sin(glat);
					register const  double glon = lon[i];
					register const double s1 = std::cos(glat);
					register double ee  = 1.0-(b*b)/(a*a);
					register const double s3 = std::sin(glon);
					register double om2 = 1.0-ee;
					register double q  = std::sqrt(1.0-ee*s0*s0);
					register const double salt = alt[i];
					register const double s4 = std::cos(glon);
					register const double s2 = std::tan(glat);
					register double d = std::sqrt(1.0-ee*s2*s2);
					register double ad = a/d;
					pos_x[i] = ad*s4+salt*s3*s1;
					pos_y[i] = ad*s3+salt*s4*1;
					pos_z[i] = a*om2*s0/q+salt*s0;
				    }
		        }


	              

                        /*
                                Reference: https://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/
                           */
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                static inline
			void forward_metod_zmm8r8(const __m512d axis,      //ellipsoid semi-maxjor axis
			                         const __m512d flat,      //elipsoid flattening [dimensionless]
						 const __m512d vp1lat,    //vector of 8 starting-points latitude [rad]
						 const __m512d vp1lon,    //vector of 8 starting-points longtitude [rad]
						 const __m512d azvf,      //vector of 8 forward azimutes [rad]
						 const __m512d dstv,      //vector of 8 distances vp1-to-vp2 [m]
						 __m512d &vp2lat,         //vector of 8 endpoints latitude [rad]
						 __m512d &vp2lon,         //vector of 8 endpoints longtitude [rad]
						 __m512d &azvb) {        //backward facing vector of 8 azimutes vp2-to-vp1 [rad]

                              const __m512d _3_14   = _mm512_set1_pd(3.1415926535897932384626);
			      const __m512d _0      = _mm512_setzero_pd();
			      const __m512d veps    = _mm512_set1_pd(0.5e-13);
			      const __m512d _1      = _mm512_set1_pd(1.0);
			      const __m512d _2      = _mm512_set1_pd(2.0);
			      const __m512d _0_375  = _mm512_set1_pd(0.375);
			      const __m512d _n3     = _mm512_set1_pd(-3.0);
			      const __m512d _4      = _mm512_set1_pd(4.0);
			      const __m512d _0_25   = _mm512_set1_pd(0.25);
			      const __m512d _0_0625 = _mm512_set1_pd(0.0625);
			      const __m512d _0_16   = _mm512_set1_pd(0.1666666666666666666667);
			      register __m512d svy = _0;
			      register __m512d cvy = _0;
			      register __m512d cvz = _0;
			      register __m512d ve  = _0;
			      register __m512d vx  = _0;
			      register __m512d vc  = _0;
			      register __m512d vy  = _0;
			      __m512d vr   = _0;
			      __m512d vsf  = _0;
			      __m512d vcf  = _0;
			      __m512d vcu  = _0;
			      __m512d vtu  = _0;
			      __m512d vsu  = _0;
			      __m512d vsa  = _0;
			      __m512d vcsa = _0;
			      __m512d vc2a = _0;
			      __m512d vd   = _0;
			      __m512d t0   = _0;
			      __m512d t1   = _0;
			      __m512d t2   = _0;
			      __mmask8 neqz = 0x0;
			      __mmask8 lteps =  0x0;

			      vr = _mm512_sub_pd(_1,flat);
			      vtu = _mm512_mul_pd(vr,_mm512_div_pd(
			                             _mm512_sin_pd(vp1lat),
						     _mm512_cos_pd(vp1lat)));
			      vcf = _mm512_cos_pd(azvf);
			      neqz = _mm512_cmp_pd_mask(vcf,_0,_CMP_NEQ_OQ);
			      azvb = _mm512_maskz_mul_pd(neqz,
			                         _mm512_atan2_pd(vtu,vcf),_2);
			      vsf = _mm512_sin_pd(azvf);
			      t0  = _mm512_fmadd_pd(vtu,vtu,_1);
			      vcu = _mm512_div_pd(_1,
			                          _mm512_sqrt_pd(t0));
			      vsu = _mm512_mul_pd(vtu,vcu);
			      t2  = _mm512_div_pd(dstv,vr);
			      vsa = _mm512_mul_pd(vcu,vsf);
			      t0  = _mm512_sub_pd(_0,vsa);
			      vc2a = _mm512_fmadd_pd(t0,vsa,_1);
			      t0 = _mm512_div_pd(_1,_mm512_mul_pd(vr,vr));
			      t1 = _mm512_sub_pd(t0,_1);
			      vx = _mm512_add_pd(_mm512_sqrt_pd(
			                         _mm512_fmadd_pd(t1,vc2a,_1)),_1);
			      vx = _mm512_div_pd(_mm512_sub_pd(vx,_2),vx);
			      vc = _mm512_sub_pd(_1,vx);
			      t1 = _mm512_mul_pd(vx,vx);
			      vd = _mm512_mul_pd(_mm512_fmsub_pd(_0_375,t0,_1),vx);
			      vc = _mm512_div_pd(_mm512_fmadd_pd(t1,_0_25,_1),vc);
			      t0 = _mm512_div_pd(axis,vc);
			      vtu = _mm512_div_pd(t2,t0);
			      vy = vtu;
	    
			label_10:
                                   vsy = _mm512_sin_pd(vy);
				   vcy = _mm512_cos_pd(vy);
				   vcz = _mm512_cos_pd(_mm512_add_pd(azvb,vy));
				   ve  = _mm512_fmsub_pd(_mm512_mul_pd(vcz,vcz),_2,_1)
				   vc  = vy;
				   vx  = _mm512_mul_pd(ve,vcy);
				   vy  = _mm512_sub_pd(_mm512_add_pd(ve,vec),_1);
				   register const __m512d c0 = _mm512_fmsub_pd(_mm512_mul_pd(vsy,vsy),_4,_3);
				   register const __m512d c1 = _mm512_fmadd_pd(_mm512_mul_pd(vy,vcz),
				                                               _mm512_mul_pd(vd,_0_16),vx);
				   register const __m512d c2 = _mm512_fmsub_pd(vd,_0_25,vcz);
				   register const __m512d c3 = _mm512_fmadd_pd(vsy,vd,vtu);
				   vy = _mm512_mul_pd(_mm512_mul_pd(c0,c1),
				                      _mm512_mul_pd(c2,c3));
				   lteps = _mm512_cmp_pd_mask(_mm512_sub_pd(vy,vc,_CMP_LE_OQ));
				   if(1==lteps) goto label_20;
				   goto label_10;
				   
			label_20:
			      
			       azvb = _mm512_fmsub_pd(_mm512_mul_pd(vcu,vcy),vcf,
			                              _mm512_mul_pd(vsu,vsy));
			       vc   = _mm512_mul_pd(vr,_mm512_sqrt_pd(
			                               _mm512_fmadd_pd(vsa,vsu,_mm512_mul_pd(azvb,azvb))));
			       t0   = _mm512_mul_pd(vcu,_mm512_mul_pd(vsy,vcf));
			       vd   = _mm512_fmadd_pd(vsu,vcy,t0);
			       vp2lat = _mm512_atan2_pd(vd,vc);
			       t0   = _mm512_mul_pd(vsu,_mm512_mul_pd(vsy,vcf));
			       vc   = _mm512_fmsub_pd(vcu,vcy,t0);
			       vx   = _mm512_atan2_pd(_mm512_mul_pd(vsy,vsf),vc);
			       t0   = _mm512_fmadd_pd(_n3,vc2a,_4);
			       t1   = _mm512_fmadd_pd(t0,flat,_4);
			       vc   = _mm512_mul_pd(t1,_mm512_mul_pd(vc2a,
			                               _mm512_mul_pd(flat,_0_16)));
			       t0   = _mm512_fmadd_pd(_mm512_mul_pd(ve,vcy),vc,vy);
			       t1   = _mm512_fmadd_pd(_mm512_mul_pd(t0,vsy),vc,vy);
			       vd   = _mm512_mul_pd(t1,vsa);
			       t2   = _mm512_mul_pd(_mm512_sub_pd(_1,vc),
			                            _mm512_mul_pd(vd,flat));
			       vp2lon = _mm512_sub_pd(_mm512_add_pd(vp1lon,vx),t2);
			       azvb = _mm512_add_pd(_mm512_atan2_pd(vsa,azvb),_3_14);
		       }



		       	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		        static inline
			void
			forward_method_u_zmm8r8_looped(const double axis,
			                               const double flat,
						       const double * __restrict plat1,
						       const double * __restrict plon1,
						       const double * __restrict pfaz,
						       const double * __restrict pdst,
						       double * __restrict plat2,
						       double * __restrict plon2,
						       double * __restrict pbaz,
						       const int32_t n) {

                              if(__builtin_expect(n<=0,0)) { return;}
			      
			      const __m512d _3_14   = _mm512_set1_pd(3.1415926535897932384626);
			      const __m512d _0      = _mm512_setzero_pd();
			      const __m512d veps    = _mm512_set1_pd(0.5e-13);
			      const __m512d _1      = _mm512_set1_pd(1.0);
			      const __m512d _2      = _mm512_set1_pd(2.0);
			      const __m512d _0_375  = _mm512_set1_pd(0.375);
			      const __m512d _n3     = _mm512_set1_pd(-3.0);
			      const __m512d _4      = _mm512_set1_pd(4.0);
			      const __m512d _0_25   = _mm512_set1_pd(0.25);
			      const __m512d _0_0625 = _mm512_set1_pd(0.0625);
			      const __m512d _0_16   = _mm512_set1_pd(0.1666666666666666666667);
			      const __m512d vaxis   = _mm512_set1_pd(axis);
			      const __m512d vflat   = _mm512_set1_pd(flat);
			      register __m512d vsy  = _0;
			      register __m512d vcz  = _0;
			      register __m512d vcy  = _0;
			      register __m512d ve   = _0;
			      double baz;
			      double sy;
			      double cz;
			      double cy;
			      double e;
			      int32_t i;
			      __mmask8 neqz = 0x0;
			      __mmask8 lteps =  0x0;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                   __m512d vr = _mm512_sub_pd(_1,vflat);
				   _mm_prefetch((const char*)&plat1[i+8],_MM_HINT_T0);
				   const __m512d vlat1 = _mm512_loadu_pd(&plat1[i]);
			           __m512d vtu = _mm512_mul_pd(vr,_mm512_div_pd(
			                                          _mm512_sin_pd(vlat1),
						                  _mm512_cos_pd(vlat1)));
			           __m512d vcf = _mm512_cos_pd(pfaz);
			           neqz = _mm512_cmp_pd_mask(vcf,_0,_CMP_NEQ_OQ);
			            __m512d vbaz =  _mm512_maskz_mul_pd(neqz,
			                                         _mm512_atan2_pd(vtu,vcf),_2));
				   _mm_prefetch((const char*)&pfaz[i+8],_MM_HINT_T0);
				   const __m512d vfaz = _mm512_loadu_pd(&pfaz[i]);
			           __m512d vsf = _mm512_sin_pd(vfaz);
			           __m512d t0  = _mm512_fmadd_pd(vtu,vtu,_1);
			           __m512d vcu = _mm512_div_pd(_1,
			                          _mm512_sqrt_pd(t0));
			           __m512d vsu = _mm512_mul_pd(vtu,vcu);
				   _mm_prefetch((const char*)&pdst[i+8],_MM_HINT_T0);
				   const __m512d dstv = _mm512_loadu_pd(&pdst[i]);
			           __m512d t2  = _mm512_div_pd(dstv,vr);
			           __m512d vsa = _mm512_mul_pd(vcu,vsf);
			           t0  = _mm512_sub_pd(_0,vsa);
			           __m512d vc2a = _mm512_fmadd_pd(t0,vsa,_1);
			           t0 = _mm512_div_pd(_1,_mm512_mul_pd(vr,vr));
			           __m512d t1 = _mm512_sub_pd(t0,_1);
			           __m512d vx = _mm512_add_pd(_mm512_sqrt_pd(
			                         _mm512_fmadd_pd(t1,vc2a,_1)),_1);
			           vx = _mm512_div_pd(_mm512_sub_pd(vx,_2),vx);
			           __m512d vc = _mm512_sub_pd(_1,vx);
			           t1 = _mm512_mul_pd(vx,vx);
			           __m512d vd = _mm512_mul_pd(_mm512_fmsub_pd(_0_375,t0,_1),vx);
			           vc = _mm512_div_pd(_mm512_fmadd_pd(t1,_0_25,_1),vc);
			           t0 = _mm512_div_pd(axis,vc);
			           vtu = _mm512_div_pd(t2,t0);
			           __m512d vy = vtu;
	    
			label_10:
                                   register  vsy = _mm512_sin_pd(vy);
				   register  vcy = _mm512_cos_pd(vy);
				   register  vcz = _mm512_cos_pd(_mm512_add_pd(vbaz,vy));
				   register  ve  = _mm512_fmsub_pd(_mm512_mul_pd(vcz,vcz),_2,_1)
				   vc  = vy;
				   vx  = _mm512_mul_pd(ve,vcy);
				   vy  = _mm512_sub_pd(_mm512_add_pd(ve,vec),_1);
				   register const __m512d c0 = _mm512_fmsub_pd(_mm512_mul_pd(vsy,vsy),_4,_3);
				   register const __m512d c1 = _mm512_fmadd_pd(_mm512_mul_pd(vy,vcz),
				                                               _mm512_mul_pd(vd,_0_16),vx);
				   register const __m512d c2 = _mm512_fmsub_pd(vd,_0_25,vcz);
				   register const __m512d c3 = _mm512_fmadd_pd(vsy,vd,vtu);
				   vy = _mm512_mul_pd(_mm512_mul_pd(c0,c1),
				                      _mm512_mul_pd(c2,c3));
				   lteps = _mm512_cmp_pd_mask(_mm512_sub_pd(vy,vc,_CMP_LE_OQ));
				   if(1==lteps) goto label_20;
				   goto label_10;
				   
			label_20:
			      
			       vbaz = _mm512_fmsub_pd(_mm512_mul_pd(vcu,vcy),vcf,
			                              _mm512_mul_pd(vsu,vsy));
			       vc   = _mm512_mul_pd(vr,_mm512_sqrt_pd(
			                               _mm512_fmadd_pd(vsa,vsu,_mm512_mul_pd(vbaz,vbaz))));
			       t0   = _mm512_mul_pd(vcu,_mm512_mul_pd(vsy,vcf));
			       vd   = _mm512_fmadd_pd(vsu,vcy,t0);
			       _mm512_storeu_pd(&plat2[i], _mm512_atan2_pd(vd,vc));
			       t0   = _mm512_mul_pd(vsu,_mm512_mul_pd(vsy,vcf));
			       vc   = _mm512_fmsub_pd(vcu,vcy,t0);
			       vx   = _mm512_atan2_pd(_mm512_mul_pd(vsy,vsf),vc);
			       t0   = _mm512_fmadd_pd(_n3,vc2a,_4);
			       t1   = _mm512_fmadd_pd(t0,flat,_4);
			       vc   = _mm512_mul_pd(t1,_mm512_mul_pd(vc2a,
			                               _mm512_mul_pd(flat,_0_16)));
			       t0   = _mm512_fmadd_pd(_mm512_mul_pd(ve,vcy),vc,vy);
			       t1   = _mm512_fmadd_pd(_mm512_mul_pd(t0,vsy),vc,vy);
			       vd   = _mm512_mul_pd(t1,vsa);
			       t2   = _mm512_mul_pd(_mm512_sub_pd(_1,vc),
			                            _mm512_mul_pd(vd,flat));
			       _mm_prefetch((const char*)&plon1[i+8],_MM_HINT_T0);
			       const __m512d vlon1 = _mm512_loadu_pd(&plon1[i]);
			       _mm512_storeu_pd(&plon2[i], _mm512_sub_pd(_mm512_add_pd(vlon1,vx),t2));
			       vbaz = _mm512_add_pd(_mm512_atan2_pd(vsa,vbaz),_3_14);
			       _mm512_storeu_pd(&pbaz[i],vbaz);
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
                                        const double r = 1.0-flat;
					const double lat1 = plat1[i];
					double tu = std::sin(lat1)/std::cos(lat1);
					const double faz = pfaz[i];
					double sf = std::sin(faz);
					double cf = std::cos(faz);
					if(cf != 0.0){
                                           baz = std::atan2(tu,cf)*2.0;
					}else {
                                           baz = 0.0;
					}
					double cu  = 1.0/std::sqrt(tu*tu+1.0);
					double su  = tu*cu;
					double sa  = cu*sf;
                                        double c2a = -sa*sa + 1.0;
                                        double x   = std::sqrt((1.0/r/r-1.0)*c2a+1.0) + 1.0;
                                        x   = (x-2.0)/x;
                                        double c   = 1.0 - x;
                                        c   = (x*x*0.25+1.0)/c;
                                        double d   = (0.375*x*x-1.0)*x;
					const double s = pdst[i];
                                        tu  = s/r/a/c;
                                        double y   = tu;
				label_10:
				          sy = std::sin(y);
                                          cy = std::cos(y);
                                          cz = std::cos(baz+y);
                                          e  = cz*cz*2.0 - 1.0;
                                          c  = y;
                                          x  = e*cy;
                                          y  = e + e - 1.0;
                                          y  = (((sy*sy*4.0_wp-3.0_wp)*y*cz*d*0.1666666666666666666667+x)*d*0.25-cz)*sy*d + tu;
					  if(std::abs(y-c)<=0.5e-13) goto label_20;
					     goto label_10;
				label_20:
				         baz   = cu*cy*cf - su*sy;
                                         c     = r*std::sqrt(sa*sa+baz*baz);
                                         d     = su*cy + cu*sy*cf;
					 plat2[i] = std::atan2(d,c);
					 c     = cu*cy - su*sy*cf;
                                         x     = std::atan2(sy*sf,c);
                                         c     = ((-3.0*c2a+4.0)*f+4.0)*c2a*f*0.0625;
                                         d     = ((e*cy*c+cz)*sy*c+y)*sa;
					 plon2[i] = plon1[i]+x-(1.0-c)*d*flat
					 baz = std::atan2(sa,baz)+3.1415926535897932384626;
					 pbaz[i] = baz;
				    }
		      }



		       	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		        static inline
			void
			forward_method_a_zmm8r8_looped(const double axis,
			                               const double flat,
						       double * __restrict __ATTR_ALIGN__(64) plat1,
						       double * __restrict __ATTR_ALIGN__(64) plon1,
						       double * __restrict __ATTR_ALIGN__(64) pfaz,
						       double * __restrict __ATTR_ALIGN__(64) pdst,
						       double * __restrict __ATTR_ALIGN__(64) plat2,
						       double * __restrict __ATTR_ALIGN__(64) plon2,
						       double * __restrict __ATTR_ALIGN__(64) pbaz,
						       const int32_t n) {

                              if(__builtin_expect(n<=0,0)) { return;}
			      
			      const __m512d _3_14   = _mm512_set1_pd(3.1415926535897932384626);
			      const __m512d _0      = _mm512_setzero_pd();
			      const __m512d veps    = _mm512_set1_pd(0.5e-13);
			      const __m512d _1      = _mm512_set1_pd(1.0);
			      const __m512d _2      = _mm512_set1_pd(2.0);
			      const __m512d _0_375  = _mm512_set1_pd(0.375);
			      const __m512d _n3     = _mm512_set1_pd(-3.0);
			      const __m512d _4      = _mm512_set1_pd(4.0);
			      const __m512d _0_25   = _mm512_set1_pd(0.25);
			      const __m512d _0_0625 = _mm512_set1_pd(0.0625);
			      const __m512d _0_16   = _mm512_set1_pd(0.1666666666666666666667);
			      const __m512d vaxis   = _mm512_set1_pd(axis);
			      const __m512d vflat   = _mm512_set1_pd(flat);
			      register __m512d vsy  = _0;
			      register __m512d vcz  = _0;
			      register __m512d vcy  = _0;
			      register __m512d ve   = _0;
			      double baz;
			      double sy;
			      double cz;
			      double cy;
			      double e;
			      int32_t i;
			      __mmask8 neqz = 0x0;
			      __mmask8 lteps =  0x0;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                              __assume_aligned(plat1,64);
			      __assume_aligned(plon1,64);
			      __assume_aligned(pfaz,64);
			      __assume_aligned(pdst,64);
			      __assume_aligned(plat2,64);
			      __assume_aligned(plon2,64);
			      __assume_aligned(pbaz,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                              plat1 = (double*)__builtin_assume_aligned(plat1,64);
			      plon1 = (double*)__builtin_assume_aligned(plon1,64);
			      pfaz  = (double*)__builtin_assume_aligned(pfaz,64);
			      pdst  = (double*)__builtin_assume_aligned(pdst,64);
			      plat2 = (double*)__builtin_assume_aligned(plat2,64);
			      plon2 = (double*)__builtin_assume_aligned(plon2,64);
			      pbaz  = (double*)__builtin_assume_aligned(pbaz,64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                   __m512d vr = _mm512_sub_pd(_1,vflat);
				   _mm_prefetch((const char*)&plat1[i+8],_MM_HINT_T0);
				   const __m512d vlat1 = _mm512_load_pd(&plat1[i]);
			           __m512d vtu = _mm512_mul_pd(vr,_mm512_div_pd(
			                                          _mm512_sin_pd(vlat1),
						                  _mm512_cos_pd(vlat1)));
			           __m512d vcf = _mm512_cos_pd(pfaz);
			           neqz = _mm512_cmp_pd_mask(vcf,_0,_CMP_NEQ_OQ);
			            __m512d vbaz =  _mm512_maskz_mul_pd(neqz,
			                                         _mm512_atan2_pd(vtu,vcf),_2));
				   _mm_prefetch((const char*)&pfaz[i+8],_MM_HINT_T0);
				   const __m512d vfaz = _mm512_load_pd(&pfaz[i]);
			           __m512d vsf = _mm512_sin_pd(vfaz);
			           __m512d t0  = _mm512_fmadd_pd(vtu,vtu,_1);
			           __m512d vcu = _mm512_div_pd(_1,
			                          _mm512_sqrt_pd(t0));
			           __m512d vsu = _mm512_mul_pd(vtu,vcu);
				   _mm_prefetch((const char*)&pdst[i+8],_MM_HINT_T0);
				   const __m512d dstv = _mm512_load_pd(&pdst[i]);
			           __m512d t2  = _mm512_div_pd(dstv,vr);
			           __m512d vsa = _mm512_mul_pd(vcu,vsf);
			           t0  = _mm512_sub_pd(_0,vsa);
			           __m512d vc2a = _mm512_fmadd_pd(t0,vsa,_1);
			           t0 = _mm512_div_pd(_1,_mm512_mul_pd(vr,vr));
			           __m512d t1 = _mm512_sub_pd(t0,_1);
			           __m512d vx = _mm512_add_pd(_mm512_sqrt_pd(
			                         _mm512_fmadd_pd(t1,vc2a,_1)),_1);
			           vx = _mm512_div_pd(_mm512_sub_pd(vx,_2),vx);
			           __m512d vc = _mm512_sub_pd(_1,vx);
			           t1 = _mm512_mul_pd(vx,vx);
			           __m512d vd = _mm512_mul_pd(_mm512_fmsub_pd(_0_375,t0,_1),vx);
			           vc = _mm512_div_pd(_mm512_fmadd_pd(t1,_0_25,_1),vc);
			           t0 = _mm512_div_pd(axis,vc);
			           vtu = _mm512_div_pd(t2,t0);
			           __m512d vy = vtu;
	    
			label_10:
                                   register  vsy = _mm512_sin_pd(vy);
				   register  vcy = _mm512_cos_pd(vy);
				   register  vcz = _mm512_cos_pd(_mm512_add_pd(vbaz,vy));
				   register  ve  = _mm512_fmsub_pd(_mm512_mul_pd(vcz,vcz),_2,_1)
				   vc  = vy;
				   vx  = _mm512_mul_pd(ve,vcy);
				   vy  = _mm512_sub_pd(_mm512_add_pd(ve,vec),_1);
				   register const __m512d c0 = _mm512_fmsub_pd(_mm512_mul_pd(vsy,vsy),_4,_3);
				   register const __m512d c1 = _mm512_fmadd_pd(_mm512_mul_pd(vy,vcz),
				                                               _mm512_mul_pd(vd,_0_16),vx);
				   register const __m512d c2 = _mm512_fmsub_pd(vd,_0_25,vcz);
				   register const __m512d c3 = _mm512_fmadd_pd(vsy,vd,vtu);
				   vy = _mm512_mul_pd(_mm512_mul_pd(c0,c1),
				                      _mm512_mul_pd(c2,c3));
				   lteps = _mm512_cmp_pd_mask(_mm512_sub_pd(vy,vc,_CMP_LE_OQ));
				   if(1==lteps) goto label_20;
				   goto label_10;
				   
			label_20:
			      
			       vbaz = _mm512_fmsub_pd(_mm512_mul_pd(vcu,vcy),vcf,
			                              _mm512_mul_pd(vsu,vsy));
			       vc   = _mm512_mul_pd(vr,_mm512_sqrt_pd(
			                               _mm512_fmadd_pd(vsa,vsu,_mm512_mul_pd(vbaz,vbaz))));
			       t0   = _mm512_mul_pd(vcu,_mm512_mul_pd(vsy,vcf));
			       vd   = _mm512_fmadd_pd(vsu,vcy,t0);
			       _mm512_store_pd(&plat2[i], _mm512_atan2_pd(vd,vc));
			       t0   = _mm512_mul_pd(vsu,_mm512_mul_pd(vsy,vcf));
			       vc   = _mm512_fmsub_pd(vcu,vcy,t0);
			       vx   = _mm512_atan2_pd(_mm512_mul_pd(vsy,vsf),vc);
			       t0   = _mm512_fmadd_pd(_n3,vc2a,_4);
			       t1   = _mm512_fmadd_pd(t0,flat,_4);
			       vc   = _mm512_mul_pd(t1,_mm512_mul_pd(vc2a,
			                               _mm512_mul_pd(flat,_0_16)));
			       t0   = _mm512_fmadd_pd(_mm512_mul_pd(ve,vcy),vc,vy);
			       t1   = _mm512_fmadd_pd(_mm512_mul_pd(t0,vsy),vc,vy);
			       vd   = _mm512_mul_pd(t1,vsa);
			       t2   = _mm512_mul_pd(_mm512_sub_pd(_1,vc),
			                            _mm512_mul_pd(vd,flat));
			       _mm_prefetch((const char*)&plon1[i+8],_MM_HINT_T0);
			       const __m512d vlon1 = _mm512_load_pd(&plon1[i]);
			       _mm512_store_pd(&plon2[i], _mm512_sub_pd(_mm512_add_pd(vlon1,vx),t2));
			       vbaz = _mm512_add_pd(_mm512_atan2_pd(vsa,vbaz),_3_14);
			       _mm512_store_pd(&pbaz[i],vbaz);
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
                                        const double r = 1.0-flat;
					const double lat1 = plat1[i];
					double tu = std::sin(lat1)/std::cos(lat1);
					const double faz = pfaz[i];
					double sf = std::sin(faz);
					double cf = std::cos(faz);
					if(cf != 0.0){
                                           baz = std::atan2(tu,cf)*2.0;
					}else {
                                           baz = 0.0;
					}
					double cu  = 1.0/std::sqrt(tu*tu+1.0);
					double su  = tu*cu;
					double sa  = cu*sf;
                                        double c2a = -sa*sa + 1.0;
                                        double x   = std::sqrt((1.0/r/r-1.0)*c2a+1.0) + 1.0;
                                        x   = (x-2.0)/x;
                                        double c   = 1.0 - x;
                                        c   = (x*x*0.25+1.0)/c;
                                        double d   = (0.375*x*x-1.0)*x;
					const double s = pdst[i];
                                        tu  = s/r/a/c;
                                        double y   = tu;
				label_10:
				          sy = std::sin(y);
                                          cy = std::cos(y);
                                          cz = std::cos(baz+y);
                                          e  = cz*cz*2.0 - 1.0;
                                          c  = y;
                                          x  = e*cy;
                                          y  = e + e - 1.0;
                                          y  = (((sy*sy*4.0_wp-3.0_wp)*y*cz*d*0.1666666666666666666667+x)*d*0.25-cz)*sy*d + tu;
					  if(std::abs(y-c)<=0.5e-13) goto label_20;
					     goto label_10;
				label_20:
				         baz   = cu*cy*cf - su*sy;
                                         c     = r*std::sqrt(sa*sa+baz*baz);
                                         d     = su*cy + cu*sy*cf;
					 plat2[i] = std::atan2(d,c);
					 c     = cu*cy - su*sy*cf;
                                         x     = std::atan2(sy*sf,c);
                                         c     = ((-3.0*c2a+4.0)*f+4.0)*c2a*f*0.0625;
                                         d     = ((e*cy*c+cz)*sy*c+y)*sa;
					 plon2[i] = plon1[i]+x-(1.0-c)*d*flat
					 baz = std::atan2(sa,baz)+3.1415926535897932384626;
					 pbaz[i] = baz;
				    }
		      }


		      /*
                           @Reference:
                                          http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
                        */
		      	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                static inline
			__m512d
			spheroid_distance_zmm8r8(const __m512d vr,
			                         const __m512d vlon1,
						 const __m512d vlat1,
						 const __m512d vlon2,
						 const __m512d vlat2) {

                               const __m512d _0 = _mm512_setzero_pd();
			       register __m512d d    = _0;
			       register __m512d zmm0 = _0; 
			       register __m512d zmm1 = _0;
			       register __m512d zmm2 = _0;
			       register __m512d zmm3 = _0;
			       register __m512d zmm4 = _0;
			       register __m512d zmm5 = _0;
			       register __m512d zmm6 = _0;
			       register __m512d t0   = _0;
			       register __m512d t1   = _0;
			       register __m512d t2   = _0;
                               register __m512d t3   = _0;
			       register __m512d t4   = _0;
			       zmm0 = _mm512_cos_pd(vlat1);
			       zmm4 = _mm512_sub_pd(vlon1,vlon2);
			       zmm1 = _mm512_sin_pd(vlat1);
			       zmm2 = _mm512_cos_pd(vlat2);
			       zmm3 = _mm512_sin_pd(vlat2);
			       zmm5 = _mm512_cos_pd(zmm4);
			       t0   = _mm512_mul_pd(zmm0,
			                        _mm512_mul_pd(zmm2,zmm5));
			       zmm6 = _mm512_sin_pd(zmm4);
			       t1   = _mm512_fmadd_pd(zmm1,zmm3,t0); 
			       t0   = _mm512_mul_pd(zmm2,zmm6);
			       t2   = _mm512_mul_pd(t0,t0); 
			       t3   = _mm512_mul_pd(zmm1,
			                        _mm512_mul_pd(zmm2,zmm5));
			       t4   = _mm512_fmsub_pd(zmm0,zmm3,t3);
			       t3   = _mm512_mul_pd(t4,t4); 
			       t0   = _mm512_sqrt_pd(_mm512_add_pd(t2,t3));
			       d    = _mm512_mul_pd(vr,_mm512_atan2_pd(t0,t1));
			       return (d);
			}


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			spheroid_distance_u_zmm8r8_looped(const double r,
			                                  double * __restrict plon1,
							  double * __restrict plat1,
							  double * __restrict plon2,
							  double * __restrict plat2,
							  double * __restrict pd,
							  const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			      register __m512d vr   = _mm512_set1_pd(r);
                              int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                  _mm_prefetch((const char*)&plat1[i+8],_MM_HINT_T0);
				  register const __m512d vlat1 = _mm512_loadu_pd(&plat1[i]);
				  register __m512d zmm0 = _mm512_cos_pd(vlat1);
				  _mm_prefetch((const char*)&plon1[i+8],_MM_HINT_T0);
				  register __m512d zmm1 = _mm512_sin_pd(vlat1);
				  _mm_prefetch((const char*)&plat2[i+8],_MM_HINT_T0);
				  register const __m512d vlat2 = _mm512_loadu_pd(&plat2[i]);
				  register __m512d zmm3 = _mm512_sin_pd(vlat2);
				  register const __m512d vlon1 = _mm512_loadu_pd(&plon1[i]);
				  _mm_prefetch((const char*)&plon2[i+8],_MM_HINT_T0);
				  register const __m512d vlon2 = _mm512_loadu_pd(&plon2[i]);
			          register __m512d zmm4 = _mm512_sub_pd(vlon1,vlon2);
			          register __m512d zmm2 = _mm512_cos_pd(vlat2);
			      	  register __m512d zmm5 = _mm512_cos_pd(zmm4);
			          register __m512d t0   = _mm512_mul_pd(zmm0,
			                                            _mm512_mul_pd(zmm2,zmm5));
			          register __m512d zmm6 = _mm512_sin_pd(zmm4);
			          register __m512d t1   = _mm512_fmadd_pd(zmm1,zmm3,t0); 
			          t0   = _mm512_mul_pd(zmm2,zmm6);
			          register __m512d t2   = _mm512_mul_pd(t0,t0); 
			          register __m512d t3   = _mm512_mul_pd(zmm1,
			                                            _mm512_mul_pd(zmm2,zmm5));
			          register __m512d t4   = _mm512_fmsub_pd(zmm0,zmm3,t3);
			          t3   = _mm512_mul_pd(t4,t4); 
			          t0   = _mm512_sqrt_pd(_mm512_add_pd(t2,t3));
			          register __m512d vd    = _mm512_mul_pd(vr,_mm512_atan2_pd(t0,t1));
				  _mm512_storeu_pd(&pd[i],vd);
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
				         const double lat1 = plat1[i];
					 double c1   = std::cos(lat1);
					 const double lon1 = plon1[i];
				         double s1   = std::sin(lat1);
					 const double lat2 = plat2[i];
				         double c2   = std::cos(lat2);
					 const double lon1 = plon1[i];
			                 double s2    = std::sin(lat2);
					 double lon2  = plon2[i];
					 double delon = lon1-lon2;
				         double clon  = std::cos(delon);
				         double t0    = s1*s2+c1*c2*clon;
				         double slon  = std::sin(delon);
					 double t1    = (c2*slon)*(c2*slon);
					 double t2    = (c1*s2-s1*c2*clon)*(c1*s2-s1*c2*clon);
					 double t3    = std::sqrt(t1+t2);
					 const double d = r*std::atan2(t3,t0);
					 pd[i] = d;
			            }
		        }


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			spheroid_distance_a_zmm8r8_looped(const double r,
			                                  double * __restrict __ATTR_ALIGN__(64) plon1,
							  double * __restrict __ATTR_ALIGN__(64) plat1,
							  double * __restrict __ATTR_ALIGN__(64) plon2,
							  double * __restrict __ATTR_ALIGN__(64) plat2,
							  double * __restrict __ATTR_ALIGN__(64) pd,
							  const int32_t n) {

                              if(__builtin_expect(n<=0,0)) {return;}
			      register __m512d vr   = _mm512_set1_pd(r);
			      int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                                __assume_aligned(plon1,64);
				__assume_aligned(plat1,64);
				__assume_aligned(plon2,64);
				__assume_aligned(plat2,64);
				__assume_aligned(pd,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                                plon1 = (double*)__builtin_assume_aligned(plon1,64);
				plat1 = (double*)__builtin_assume_aligned(plat1,64);
				plon2 = (double*)__builtin_assume_aligned(plon2,64);
				plat2 = (double*)__builtin_assume_aligned(plat2,64);
				pd    = (double*)__builtin_assume_aligned(pd,64);
#endif			      
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                  _mm_prefetch((const char*)&plat1[i+8],_MM_HINT_T0);
				  register const __m512d vlat1 = _mm512_load_pd(&plat1[i]);
				  register __m512d zmm0 = _mm512_cos_pd(vlat1);
				  _mm_prefetch((const char*)&plon1[i+8],_MM_HINT_T0);
				  register __m512d zmm1 = _mm512_sin_pd(vlat1);
				  _mm_prefetch((const char*)&plat2[i+8],_MM_HINT_T0);
				  register const __m512d vlat2 = _mm512_load_pd(&plat2[i]);
				  register __m512d zmm3 = _mm512_sin_pd(vlat2);
				  register const __m512d vlon1 = _mm512_load_pd(&plon1[i]);
				  _mm_prefetch((const char*)&plon2[i+8],_MM_HINT_T0);
				  register const __m512d vlon2 = _mm512_load_pd(&plon2[i]);
			          register __m512d zmm4 = _mm512_sub_pd(vlon1,vlon2);
			          register __m512d zmm2 = _mm512_cos_pd(vlat2);
			      	  register __m512d zmm5 = _mm512_cos_pd(zmm4);
			          register __m512d t0   = _mm512_mul_pd(zmm0,
			                                            _mm512_mul_pd(zmm2,zmm5));
			          register __m512d zmm6 = _mm512_sin_pd(zmm4);
			          register __m512d t1   = _mm512_fmadd_pd(zmm1,zmm3,t0); 
			          t0   = _mm512_mul_pd(zmm2,zmm6);
			          register __m512d t2   = _mm512_mul_pd(t0,t0); 
			          register __m512d t3   = _mm512_mul_pd(zmm1,
			                                            _mm512_mul_pd(zmm2,zmm5));
			          register __m512d t4   = _mm512_fmsub_pd(zmm0,zmm3,t3);
			          t3   = _mm512_mul_pd(t4,t4); 
			          t0   = _mm512_sqrt_pd(_mm512_add_pd(t2,t3));
			          register __m512d vd    = _mm512_mul_pd(vr,_mm512_atan2_pd(t0,t1));
				  _mm512_store_pd(&pd[i],vd);
			      }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
				         const double lat1 = plat1[i];
					 double c1   = std::cos(lat1);
					 const double lon1 = plon1[i];
				         double s1   = std::sin(lat1);
					 const double lat2 = plat2[i];
				         double c2   = std::cos(lat2);
					 const double lon1 = plon1[i];
			                 double s2    = std::sin(lat2);
					 double lon2  = plon2[i];
					 double delon = lon1-lon2;
				         double clon  = std::cos(delon);
				         double t0    = s1*s2+c1*c2*clon;
				         double slon  = std::sin(delon);
					 double t1    = (c2*slon)*(c2*slon);
					 double t2    = (c1*s2-s1*c2*clon)*(c1*s2-s1*c2*clon);
					 double t3    = std::sqrt(t1+t2);
					 const double d = r*std::atan2(t3,t0);
					 pd[i] = d;
			            }
		        }


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                static inline
			__m512d
			geocentric_radius_zmm8r8(const __m512d va,
			                         const __m512d vb,
						 const __m512d vlat) {

			    const __m512d _0 = _mm512_setzero_pd();
			    register __m512d vr = _0;
			    register __m512d vclat = _0;
			    register __m512d vslat = _0;
			    register __m512d vaa   = _0;
			    register __m512d vbb   = _0;
			    register __m512d vden  = _0;
			    register __m512d vnum  = _0;
			    register __m512d t0    = _0;
			    register __m512d t1    = _0;
			    register __m512d t2    = _0;
                            __mmask8 aeq0 = 0x0;
			    __mmask8 beq0 = 0x0;
			    aeq0 = _mm512_cmp_pd_mask(va,_0,_CMP_EQ_OQ);
			    beq0 = _mm512_cmp_pd_mask(vb,_0,_CMP_EQ_OQ);
			    if(aeq0 && beq0) { return (vr);}
			    t0   = _mm512_cos_pd(vlat);
			    vaa  = _mm512_mul_pd(va,va);
			    vclat = _mm512_mul_pd(t0,t0);
			    t0   = _mm512_sin_pd(vlat);
			    vbb  = _mm512_mul_pd(vb,vb);
			    vslat = _mm512_mul_pd(t0,t0);
			    t1  = __mm512_mul_pd(vslat,_mm512_mul_pd(vbb,vbb));
			    vnum = _mm512_fmadd_pd(vclat,_mm512_mul_pd(vaa,vaa),t1);
			    t2 = _mm512_mul_pd(vslat,vbb);
			    vden = _mm512_fmadd_pd(vclat,vaa,t2);
			    vr = _mm512_sqrt_pd(vnum,vden);
			    return (vr);
		       }


		       
		       	__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			geocentric_radius_u_zmm8r8_looped(const double a,
			                                  const double * __restrict pb,
							  const double * __restrict plat,
							  double * __restrict pr,
							  const int32_t n) {
                               if(__builtin_expect(n<=0),0) {return;}
                               const __m512d _0 = _mm512_setzero_pd();
			       const register __m512d va = _mm512_set1_pd(a);
			       int32_t i;
			       // Error checking code removed!!
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                  _mm_prefetch((const char*)&plat[i+8],_MM_HINT_T0);
				  register const __m512d vlat = _mm512_loadu_pd(plat[i]);
				  register __m512d t0   = _mm512_cos_pd(vlat);
			          register __m512d vaa  = _mm512_mul_pd(va,va);
			          register __m512d vclat = _mm512_mul_pd(t0,t0);
			          t0   = _mm512_sin_pd(vlat);
			          _mm_prefetch((const char*)&pb[i+8],_MM_HINT_T0);
				  register const __m512d vb = _mm512_loadu_pd(&pb[i]);
			          register __m512d vbb  = _mm512_mul_pd(vb,vb);
			          register __m512d vslat = _mm512_mul_pd(t0,t0);
			          register __m512d t1  = __mm512_mul_pd(vslat,_mm512_mul_pd(vbb,vbb));
			          register __m512d vnum = _mm512_fmadd_pd(vclat,_mm512_mul_pd(vaa,vaa),t1);
			          register __m512d t2 = _mm512_mul_pd(vslat,vbb);
			          register __m512d vden = _mm512_fmadd_pd(vclat,vaa,t2);
			          const register __m512d vr = _mm512_sqrt_pd(vnum,vden);
				  _mm512_storeu_pd(&pr[i],vr);
			     }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
                                        const double lat = plat[i];
					double t0  = std::cos(lat);
					double clat = t0*t0;
					const double b = pb[i];
					double t1 = std::sin(lat);
					double slat = t1*t1;
					double aa = a*a;
					double bb = b*b;
					double num = clat*aa*aa+slat*bb*bb;
					double den = clat*aa+slat*bb;
					const double r = std::sqrt(num/den);
					pr[i] = r;
				    }
			      
			}

                        __ATTR_ALWAYS_INLINE__
			__ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			geocentric_radius_a_zmm8r8_looped(const double a,
			                                  const double * __restrict __ATTR_ALIGN__(64) pb,
							  const double * __restrict __ATTR_ALIGN__(64) plat,
							  double * __restrict __ATTR_ALIGN__(64) pr,
							  const int32_t n) {
                               if(__builtin_expect(n<=0),0) {return;}
                               const __m512d _0 = _mm512_setzero_pd();
			       const register __m512d va = _mm512_set1_pd(a);
			       int32_t i;
			       // Error checking code removed!!
#if defined(__INTEL_COMPILER) || defined(__ICC)
                               __assume_aligned(pb,64);
			       __assume_aligned(plat,64);
			       __assume_aligned(pr,64);
#elif  defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                               pb   = (double*)__builtin_assume_aligned(pb,64);
			       plat = (double*)__builtin_assume_aligned(plat,64);
			       pr   = (double*)__builtin_assume_aligned(pr,64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
                              for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                  _mm_prefetch((const char*)&plat[i+8],_MM_HINT_T0);
				  register const __m512d vlat = _mm512_load_pd(plat[i]);
				  register __m512d t0   = _mm512_cos_pd(vlat);
			          register __m512d vaa  = _mm512_mul_pd(va,va);
			          register __m512d vclat = _mm512_mul_pd(t0,t0);
			          t0   = _mm512_sin_pd(vlat);
			          _mm_prefetch((const char*)&pb[i+8],_MM_HINT_T0);
				  register const __m512d vb = _mm512_load_pd(&pb[i]);
			          register __m512d vbb  = _mm512_mul_pd(vb,vb);
			          register __m512d vslat = _mm512_mul_pd(t0,t0);
			          register __m512d t1  = __mm512_mul_pd(vslat,_mm512_mul_pd(vbb,vbb));
			          register __m512d vnum = _mm512_fmadd_pd(vclat,_mm512_mul_pd(vaa,vaa),t1);
			          register __m512d t2 = _mm512_mul_pd(vslat,vbb);
			          register __m512d vden = _mm512_fmadd_pd(vclat,vaa,t2);
			          const register __m512d vr = _mm512_sqrt_pd(vnum,vden);
				  _mm512_store_pd(&pr[i],vr);
			     }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
                                    for(; i != n; ++i) {
                                        const double lat = plat[i];
					double t0  = std::cos(lat);
					double clat = t0*t0;
					const double b = pb[i];
					double t1 = std::sin(lat);
					double slat = t1*t1;
					double aa = a*a;
					double bb = b*b;
					double num = clat*aa*aa+slat*bb*bb;
					double den = clat*aa+slat*bb;
					const double r = std::sqrt(num/den);
					pr[i] = r;
				    }
			      
			}
						     
     }

}














#endif /*__GMS_GEODESY_AVX512_HPP__*/
