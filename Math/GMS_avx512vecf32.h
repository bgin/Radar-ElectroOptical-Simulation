

#ifndef __GMS_AVX512VECF32_H__
#define __GMS_AVX512VECF32_H__


namespace file_info {

    const unsigned int gGMS_AVX512VECF32_MAJOR = 1U;
    const unsigned int gGMS_AVX512VECF32_MINOR = 0U;
    const unsigned int gGMS_AVX512VECF32_MICRO = 0U;
    const unsigned int gGMS_AVX512VECF32_FULLVER =
       1000U*gGMS_AVX512VECF32_MAJOR+
       100U*gGMS_AVX512VECF32_MINOR +
       10U*gGMS_AVX512VECF32_MICRO;
    const char * const pgGMS_AVX512VECF32_CREATE_DATE = "15-01-2020 10:44 +00200 (WED 15 JAN 2020 GMT+2)";
    const char * const pgGMS_AVX512VECF32_BUILD_DATE  = __DATE__ " " __TIME__;
    const char * const pgGMS_AVX512VECF32_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_AVX512VECF32_SYNOPSYS    = "Wrapper class around __m512 data type.";
}

#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

          namespace math {

                struct AVX512Vec16 {

		     __m512  __ATTR_ALIGN__(64) m_v16;

		     __ATTR_COLD__
		     __ATTR_ALIGN__(16)
		  inline  AVX512Vec16() {
                          m_v16 = _mm512_setzero_ps();
		     }

		     __ATTR_COLD__
		     __ATTR_ALIGN__(16)
                  inline AVX512Vec16(const float * __restrict __ATTR_ALIGN__(64) v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        v = (const float*)__builtin_assume_aligned(v,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(v,64);
#endif
                        m_v16 = _mm512_load_ps(&v[0]);
		  }

		    __ATTR_COLD__
		    __ATTR_ALIGN__(16)
		    __ATTR_VECTORCALL__
		  inline AVX512Vec16( const float s0,
		                      const float s1,
				      const float s2,
				      const float s3,
				      const float s4,
				      const float s5,
				      const float s6,
				      const float s7,
				      const float s8,
				      const float s9,
				      const float s10;
				      const float s11,
				      const float s12,
				      const float s13,
				      const float s14,
				      const float s15) {

                      m_v16 = _mm512_setr_ps(s0,s1,s2,s3,
		                             s4,s5,s6,s7,
					     s8,s9,s10,s11,
					     s12,s13,s14,s15);
		 }

		   __ATTR_COLD__
		   __ATTR_ALIGN__(16)
		 inline explicit AVX512Vec16(const float s0) {

		      m_v16 = _mm512_set1_ps(s0);
		 }

		 __ATTR_COLD__
		 __ATTR_ALIGN__(16)
		 __ATTR_VECTORCALL__
	       inline AVX512Vec16(const __m512 v) {

	              m_v16 = v;
	        }

		 __ATTR_COLD__
		 __ATTR_ALIGN__(16)
		 __ATTR_VECTORCALL__
	       inline AVX512Vec16(const __m512i v) {

	              m_v16 = _mm512_castsi512_ps(v);
	        }

		 __ATTR_COLD__
		 __ATTR_ALIGN__(16)
	       inline AVX512Vec16(const AVX512Vec16 &x) {

	              m_v16 = x.m_v16;
	        }

		__ATTR_COLD__
		__ATTR_ALIGN__(16)
		__ATTR_VECTORCALL__
	      inline AVX512Vec16(const __m256 lo,
	                         const __m256 hi) {

		      m_v16 = _mm512_insert32x8_ps(m_v16,lo,0);
		      m_v16 = _mm512_insert32x8_ps(m_v16,hi,1);
	        }

              ~AVX512Vec16() {}

	      __ATTR_COLD__
	      __ATTR_ALIGN__(16)
	    inline __m512 lo_part() const {

	             return (_mm512_extract32x8_ps(m_v16,0));
	       }

	      __ATTR_COLD__
	      __ATTR_ALIGN__(16)
	    inline __m512 hi_part() const {
  
                     return (_mm512_extract32x8_ps(m_v16,1));
	       }

	            //
		    // Load-store functions
		    //
		    // Address aligned on 64-byte boundary

	     __ATTR_COLD__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
	   inline AVX512Vec16 &
	   load_a(const float * __restrict __ATTR_ALIGN__(64) addr) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                    addr = (const float*)__builtin_assume_aligned(addr,64);
#elif defined __ICC || defined __INTEL_COMPILER
                    __assume_aligned(addr,64);
#endif
                    m_v16 = _mm512_load_ps(&addr[0]);
		    return (*this);
	       }

	     __ATTR_COLD__
	     __ATTR_ALIGN__(16)
	     __ATTR_VECTORCALL__
	  inline AVX512Vec16 &
	  load_u(const float * __restrict addr) {

	            m_v16 = _mm512_loadu_ps(&addr[0]);
		    return (*this);
	       }

              // Address argument should be aligned on 64-byte boundary
	    __ATTR_COLD__
	    __ATTR_ALIGN__(16)
	    inline void
	    store_a(float * __restrict __ATTR_ALIGN__(64) addr) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                    addr = (const float*)__builtin_assume_aligned(addr,64);
#elif defined __ICC || defined __INTEL_COMPILER
                    __assume_aligned(addr,64);
#endif
                    _mm512_store_ps(&addr[0],m_v16);
	       }

	     __ATTR_COLD__
	     __ATTR_ALIGN__(16)
	     inline void
	     store_u(float * __restrict addr) {

                    _mm512_storeu_ps(&addr[0],m_v16);
	       }

	     __ATTR_COLD__
	     __ATTR_ALIGN__(16)
	    inline void
	    stream(float * __restrict addr) {

	            _mm512_stream_ps(&addr[0],m_v16);
		    _mm_sfence();
	       }

             __ATTR_COLD__
             __ATTR_ALIGN__(16)
             inline float
	     extract_f32x1(const uint32_t pos) {
                  float t[16] __ATTR_ALIGN__(64) = {};
		  store_a(t);
		  return (t[pos&0xF]);
	       }

              // Inserts a single float in location 0...15
	      __ATTR_COLD__
              __ATTR_ALIGN__(32)
	      inline AVX512Vec16 &
	      insert_f32x1(const uint32_t pos,
	                   const float s0) {
                    float t[16] __ATTR_ALIGN__(64) = {};
		    store_a(t);
                    t[pos&0xF] = s0;
		    load_a(t);
		    return (*this);
	       }

	       __ATTR_HOT__
	       __ATTR_ALIGN__(16)
	       inline AVX512Vec16 &
	       operator=(const AVX512Vec16 &x) {
                    if(__builtin_expect(this==&x),0)
		       return (*this);
                    m_v16 = x.m_v16
		    return (*this);
	       }

	       __ATTR_HOT__
	       __ATTR_ALIGN__(16)
	       __ATTR_VECTORCALL__
	       inline AVX512Vec16 &
	       operator=(const AVX512vec16 x) {
                    if(__builtin_expect(this==&x),0)
		       return (*this);
                    m_v16 = x.m_v16
		    return (*this);
	       }

	       operator __m512 () const { return (m_v16);}
	       
	   } __ATTR_ALIGN__(64);



	   __ATTR_COLD__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   __m256 extract_f32x8(const AVX512Vec16 x,
	                        const int32_t imm8) {
                    return (_mm512_extractf32x8_ps(x,imm8));
		}

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 select(const AVX512Vec16 x,
	                      const AVX512Vec16 y,
			      const __mmask16 pred) {
                     return (_mm512_mask_blend_ps(pred,x,y));
		}

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 max(const AVX512Vec16 x,
	                   const AVX512Vec16 y) {
                     return (_mm512_max_ps(x,y));
		}

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 min(const AVX512Vec16 x,
	                   const AVX512Vec16 y) {
                     return (_mm512_min_ps(x,y));
		}

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 abs(const AVX512Vec16 x) {
                     return (_mm512_abs_ps(x));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 sqrt(const AVX512Vec16 x) {
                     return (_mm512_sqrt_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 rsqrt(const AVX512Vec16 x) {
                     return (_mm512_rsqrt_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 cbrt(const AVX512Vec16 x) {
                     return (_mm512_cbrt_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   float reduce_add(const AVX512Vec16 x) {
                     return (_mm512_reduce_add_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   float reduce_mul(const AVX512Vec16 x) {
                     return (_mm512_reduce_mul_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   float reduce_max(const AVX512Vec16 x) {
                      return (_mm512_reduce_max_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   float reduce_min(const AVX512Vec16 x) {
                      return (_mm512_reduce_min_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 ceil(const AVX512Vec16 x) {
                      return (_mm512_ceil_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 floor(const AVX512Vec16 x) {
                      return (_mm512_floor_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 round(const AVX512Vec16 x,
	                     const int32_t imm) {
                      return (_mm512_roundscale_ps(x,imm));
	      }

	   // Masked trigo vector functions.

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_acos(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_acos_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_acosh(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_acosh_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_asin(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_asin_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_asinh(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_asinh_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_atan(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_atan_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_atan2(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y,
				 const AVX512Vec16 z) {
                      return (_mm512_mask_atan2_ps(x,k,y,z));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_atanh(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_atanh_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_cosh(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_cosh_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_sin(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_sin_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_sind(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_sind_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_sincos( AVX512Vec16* cos_res,
	                            const AVX512Vec16  sin_src,
				    const AVX512Vec16  cos_src,
	                            const __mmask16 k,
				    const AVX512Vec16 a) {
                      return (_mm512_mask_sincos_ps(&cos_res,sin_src,cos_src,k,a));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_sinh(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_sinh_ps(x,k,y));
	      }


	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_tan(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_tan_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_tanh(const AVX512Vec16 x,
	                         const __mmask16 k,
				 const AVX512Vec16 y) {
                      return (_mm512_mask_tanh_ps(x,k,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_permutevar(const AVX512Vec16 src,
	                               const __mmask16 k,
				       const AVX512Vec16 a,
				       const __m512i b) {
                       return (_mm512_mask_permutevar_ps(src,k,a,b));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	    AVX512Vec16 nearbyint(const AVX512Vec16 x) {
                       return (_mm512_nearnyint_ps(x));
	     }

	   

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 sin(const AVX512Vec16 x) {
                      return (_mm512_sin_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 sind(const AVX512Vec16 x) {
                      return (_mm512_sind_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 cos(const AVX512Vec16 x) {
                      return (_mm512_cos_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 cosd(const AVX512Vec16 x) {
                      return (_mm512_cosd_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 sinh(const AVX512Vec16 x) {
                      return (_mm512_sinh_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 cosh(const AVX512Vec16 x) {
                      return (_mm512_cosh_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 tan(const AVX512Vec16 x) {
                      return (_mm512_tan_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 tanh(const AVX512Vec16 x) {
                      return (_mm512_tanh_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 asin(const AVX512Vec16 x) {
                      return (_mm512_asin_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 asinh(const AVX512Vec16 x) {
                      return (_mm512_asinh_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 acos(const AVX512Vec16 x) {
                      return (_mm512_acos_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 acosh(const AVX512Vec16 x) {
                      return (_mm512_acosh_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 atan(const AVX512Vec16 x) {
                      return (_mm512_atan_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 atanh(const AVX512Vec16 x) {
                      return (_mm512_atanh_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 log(const AVX512Vec16 x) {
                      return (_mm512_log_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 exp(const AVX512Vec16 x) {
                      return (_mm512_exp_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 atan2(const AVX512Vec16 x,
	                     const AVX512Vec16 y) {
                      return (_mm512_atan2_ps(x,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 hypot(const AVX512Vec16 x,
	                     const AVX512Vec16 y) {
                      return (_mm512_hypot_ps(x,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 fmadd(const AVX512Vec16 x,
	                     const AVX512Vec16 y,
			     const AVX512Vec16 z) {
                      return (_mm512_fmadd_ps(x,y,z));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 fmadsubb(const AVX512Vec16 x,
	                        const AVX512Vec16 y,
			        const AVX512Vec16 z) {
                      return (_mm512_fmadsubb_ps(x,y,z));
	       }

           __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 fmsub(const AVX512Vec16 x,
	                        const AVX512Vec16 y,
			        const AVX512Vec16 z) {
                      return (_mm512_fmsub_ps(x,y,z));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 fmsub(const AVX512Vec16 x,
	                        const AVX512Vec16 y,
			        const AVX512Vec16 z) {
                      return (_mm512_fmsub_ps(x,y,z));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 fmsubadd(const AVX512Vec16 x,
	                        const AVX512Vec16 y,
			        const AVX512Vec16 z) {
                      return (_mm512_fmsubadd_ps(x,y,z));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 fnmadd(const AVX512Vec16 x,
	                        const AVX512Vec16 y,
			        const AVX512Vec16 z) {
                      return (_mm512_fnmadd_ps(x,y,z));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 fnmsub(const AVX512Vec16 x,
	                        const AVX512Vec16 y,
			        const AVX512Vec16 z) {
                      return (_mm512_fnmsub_ps(x,y,z));
	       }

	       // Statistical/Probability Functions

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 cdfnorm(const AVX512Vec16 x) {
                       return (_mm512_cdfnorm_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_cdfnorm(const AVX512Vec16 x,
	                            const __mmask16 k,
				    const AVX512Vec16 a) {
                       return (_mm512_mask_cdfnorm_ps(x,k,a));
	      }

	   
           __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 cdfnorminv(const AVX512Vec16 x) {
                       return (_mm512_cdfnorminv_ps(x));
	      }

	   
	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_cdfnorminv(const AVX512Vec16 x,
	                            const __mmask16 k,
				    const AVX512Vec16 a) {
                       return (_mm512_mask_cdfnorminv_ps(x,k,a));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 erf(const AVX512Vec16 x) {
                       return (_mm512_erf_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_erf(const AVX512Vec16 x,
	                            const __mmask16 k,
				    const AVX512Vec16 a) {
                       return (_mm512_mask_erf_ps(x,k,a));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 erfc(const AVX512Vec16 x) {
                       return (_mm512_erfc_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_erfc(const AVX512Vec16 x,
	                            const __mmask16 k,
				    const AVX512Vec16 a) {
                       return (_mm512_mask_erfc_ps(x,k,a));
	      }

	    __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 erfinv(const AVX512Vec16 x) {
                       return (_mm512_erfinv_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_erfinv(const AVX512Vec16 x,
	                            const __mmask16 k,
				    const AVX512Vec16 a) {
                       return (_mm512_mask_erfinv_ps(x,k,a));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 erfcinv(const AVX512Vec16 x) {
                       return (_mm512_erfcinv_ps(x));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 mask_erfcinv(const AVX512Vec16 x,
	                            const __mmask16 k,
				    const AVX512Vec16 a) {
                       return (_mm512_mask_erfcinv_ps(x,k,a));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+(const AVX512Vec16 x,
	                         const AVX512Vec16 y) {
                       return (_mm512_add_ps(x,y));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+(const AVX512Vec16 x,
	                         const __m512 y) {
                       return (_mm512_add_ps(x,y));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+(const __m512 x,
	                         const AVX512Vec16 y) {
                       return (_mm512_add_ps(x,y));
	       }

	   
	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16  operator+(const AVX512Vec16 x,
	                          const float s) {
                        return (_mm512_add_ps(x,
			        _mm512_set1_ps(s)));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+(const float s,
	                         const AVX512Vec16 x) {
                      return (_mm512_add_ps(
		                       _mm512_set1_ps(s),
				       x));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+=(AVX512Vec16 x,
	                          const AVX512Vec16 y) {
                      x = x+y;
		      return (x);
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+=(AVX512Vec16 x,
	                          const __m512 y) {
                      x = x+y;
		      return (x);
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+=(const __m512 x,
	                          AVX512Vec16 y) {
                     y = x+y;
		     return (y);
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+=(AVX512Vec16 x,
	                          const float s) {
                     x = x+s;
		     return (x);
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator+=(const float s,
	                          AVX512Vec16 x) {
                     x = s+x;
		     return (x);
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-(const AVX512Vec16 x,
	                         const AVX512Vec16 y) {
                       return (_mm512_sub_ps(x,y));
	       }
 
	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-(const AVX512Vec16 x,
	                         const __m512 y) {
                       return (_mm512_sub_ps(x,y));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-(const __m512 x,
	                         const AVX512Vec16 y) {
                       return (_mm512_sub_ps(x,y));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-(const AVX512Vec16 x,
	                         const float s) {
                       return (_mm512_sub_ps(x,
		                     _mm512_set1_ps(s)));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-(const float s,
	                         AVX512Vec16 x) {
                       return (_mm512_sub_ps(
		                     _mm512_set1_ps(s),
				     x));
	       }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-=(AVX512Vec16 x,
	                          const AVX512Vec16 y) {
                       x = x-y;
		       return (x);
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-=(AVX512Vec16 x,
	                          const __m512 y) {

                       x = x-y;
		       return (x);
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-=(const __m512 x,
	                          AVX512Vec16 y) {
                       y = x-y;
		       return (y);
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-=(AVX512Vec16 x,
	                          const float s) {
                       x = x-s;
		       return (x);
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator-=(const float s,
	                          AVX512Vec16 x) {
                       x = s-x;
		       return (x);
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*(const AVX512Vec16 x,
	                         const AVX512Vec16 y) {
                       return (_mm512_mul_ps(x,y));
	      }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*(const AVX512Vec16 x,
	                         const __m512 y) {
                       return (_mm512_mul_ps(x,y));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*(const __m512 y,
	                         const AVX512Vec16 x) {
                       return (_mm512_mul_ps(y,x));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*(const AVX512Vec16 x,
	                         const float s) {
                       return (_mm512_mul_ps(x,
		                         _mm512_set1_ps(s)));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*(const float s,
	                         const AVX512Vec16 x) {
                       return (_mm512_mul_ps(
		                           _mm512_set1_ps(s),
					     x));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*=( AVX512Vec16 x,
	                          const AVX512Vec16 y) {

		        x = x*y;
			return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*=(AVX512Vec16 x,
	                          const __m512 y) {

			x = x*y;
			return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*=(const __m512 x,
	                          AVX512Vec16 y) {

			y = x*y;
			return (y);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*=(AVX512Vec16 x,
	                          const float s) {

			x = x*s;
			return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator*=(const float s,
	                          AVX512Vec16 x) {

			 x = s*x;
			 return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/(const AVX512Vec16 x,
	                         const AVX512Vec16 y) {

		       return (_mm512_div_ps(x,y));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/(const AVX512Vec16 x,
	                         const __m512 y) {

			return (_mm512_div_ps(x,y));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/(const __m512 x,
	                         const AVX512Vec16 y) {

		         return (_mm512_div_ps(x,y));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/(const AVX512Vec16 x,
	                         const float s) {

			 return (_mm512_div_ps(x,
			               _mm512_set1_ps(s)));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/(const float s,
	                         const AVX512Vec16 x) {

			 return (_mm512_div_ps(
			                    _mm512_set1_ps(s),x));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/=(AVX512Vec16 x,
	                          const AVX512Vec16 y) {

			 x = x/y;
			 return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/=(AVX512Vec16 x,
	                          const __m512 y) {

			  x = x/y;
			  return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/=(const __m512 x,
	                          AVX512Vec16 y) {

			  y = x/y;
			  return (y);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/=(AVX512Vec16 x,
	                          const float s) {

			  x = x/s;
			  return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	   AVX512Vec16 operator/=(const float s,
	                          AVX512Vec16 x) {

			  x = s/x;
			  return (x);
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator==(const AVX512Vec16 x,
	                       const AVX512Vec16 y) {

			  return (_mm512_cmp_ps_mask(x,y,_CMP_EQ_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator==(const AVX512Vec16 x,
	                       const float s) {

			   return (_mm512_cmp_ps_mask(x,
			                        _mm512_set1_ps(s),
						          _CMP_EQ_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator==(const float s,
	                       const AVX512Vec16 x) {

                           return (_mm512_cmp_ps_mask(
			                      _mm512_set1_ps(s),x,
					                  _CMP_EQ_OQ));
					      
	     }

	  
             __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator!=(const AVX512Vec16 x,
	                       const AVX512Vec16 y) {

			  return (_mm512_cmp_ps_mask(x,y,_CMP_NEQ_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator!=(const AVX512Vec16 x,
	                       const float s) {

			   return (_mm512_cmp_ps_mask(x,
			                        _mm512_set1_ps(s),
						          _CMP_NEQ_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator!=(const float s,
	                       const AVX512Vec16 x) {

                           return (_mm512_cmp_ps_mask(
			                      _mm512_set1_ps(s),x,
					                  _CMP_NEQ_OQ));
					      
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator>=(const AVX512Vec16 x,
	                       const AVX512Vec16 y) {

			  return (_mm512_cmp_ps_mask(x,y,_CMP_GE_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator>=(const AVX512Vec16 x,
	                       const float s) {

			   return (_mm512_cmp_ps_mask(x,
			                        _mm512_set1_ps(s),
						          _CMP_GE_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator>=(const float s,
	                       const AVX512Vec16 x) {

                           return (_mm512_cmp_ps_mask(
			                      _mm512_set1_ps(s),x,
					                  _CMP_GE_OQ));
					      
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator>(const AVX512Vec16 x,
	                       const AVX512Vec16 y) {

			  return (_mm512_cmp_ps_mask(x,y,_CMP_GT_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator>(const AVX512Vec16 x,
	                       const float s) {

			   return (_mm512_cmp_ps_mask(x,
			                        _mm512_set1_ps(s),
						          _CMP_GT_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator>(const float s,
	                       const AVX512Vec16 x) {

                           return (_mm512_cmp_ps_mask(
			                      _mm512_set1_ps(s),x,
					                  _CMP_GT_OQ));
					      
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator<=(const AVX512Vec16 x,
	                       const AVX512Vec16 y) {

			  return (_mm512_cmp_ps_mask(x,y,_CMP_LE_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator<=(const AVX512Vec16 x,
	                       const float s) {

			   return (_mm512_cmp_ps_mask(x,
			                        _mm512_set1_ps(s),
						          _CMP_LE_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator<=(const float s,
	                       const AVX512Vec16 x) {

                           return (_mm512_cmp_ps_mask(
			                      _mm512_set1_ps(s),x,
					                  _CMP_LE_OQ));
					      
	     }

	    __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator<(const AVX512Vec16 x,
	                       const AVX512Vec16 y) {

			  return (_mm512_cmp_ps_mask(x,y,_CMP_LT_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator<(const AVX512Vec16 x,
	                       const float s) {

			   return (_mm512_cmp_ps_mask(x,
			                        _mm512_set1_ps(s),
						          _CMP_LT_OQ));
	     }

	   __ATTR_HOT__
	   __ATTR_ALIGN__(16)
	   __ATTR_VECTORCALL__
	   static inline
	  __mmask16 operator<(const float s,
	                       const AVX512Vec16 x) {

                           return (_mm512_cmp_ps_mask(
			                      _mm512_set1_ps(s),x,
					                  _CMP_LT_OQ));
					      
	     }
	     

	  
    } // math

} // gms





#endif /*__GMS_AVX512VECF32_H__*/
