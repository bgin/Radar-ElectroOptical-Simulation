
#ifndef __GMS_ROTATIONS_AVX2_HELPERS_H__
#define __GMS_ROTATIONS_AVX2_HELPERS_H__



#include <immintrin.h>


                               const __m256 v8_0      = _mm256_set1_ps(0.0F);
		               const __m256 v8_2      = _mm256_set1_ps(2.0f);
                               const __m256 v8_n1     = _mm256_set1_ps(-1.0f);
			       const __m256 v8_n2     = _mm256_set1_ps(-2.0F);
			       const __m256 v8_1o2    = _mm256_set1_ps(0.5F);
			       const __m256 v8_spi    = _mm256_set1_ps(1.7724538509055160272982F);
			       const __m256 v8_s6pi   = _mm256_set1_ps(1.381976597885341917061F);
			       const __m256 v8_a      = _mm256_set1_ps(1.9257490199582527754939F);
			       const __m256 v8_ap     = _mm256_set1_ps(2.1450293971110256000775F);
			       const __m256 v8_sc     = _mm256_set1_ps(0.8977727869612861128953F);
			       const __m256 v8_beta   = _mm256_set1_ps(0.9628745099791263877469F);
			       const __m256 v8_r1     = _mm256_set1_ps(1.3306700394914687909256F);
			       const __m256 v8_r2     = _mm256_set1_ps(1.4142135623730950488017F);
			       const __m256 v8_pi12   = _mm256_set1_ps(0.2617993877991494365386F);
			       const __m256 v8_prek   = _mm256_set1_ps(1.6434564029725030125017F);
			       const __m256 v8_pi     = _mm256_set1_ps(3.1415926535897932384626F);
			       const __m256 v8_2pi    = _mm256_set1_ps(6.2831853071795864769253F);
			      
			       const __m256d v4_1o2    = _mm256_set1_pd(0.5);
                               const __m256d v4_0      = _mm256_set1_pd(0.0);
		               const __m256d v4_2      = _mm256_set1_pd(2.0);
			       const __m256d v4_n1     = _mm256_set1_pd(-1.0);
			       const __m256d  v4_pi     = _mm256_set1_pd(3.1415926535897932384626);
			       const __m256d  v4_2pi    = _mm256_set1_pd(6.2831853071795864769253);
			       const __m256d v4_spi    = _mm256_set1_pd(1.7724538509055160272982);
			       const __m256d v4_s6pi   = _mm256_set1_pd(1.381976597885341917061);
			       const __m256d v4_a      = _mm256_set1_pd(1.9257490199582527754939);
			       const __m256d v4_ap     = _mm256_set1_pd(2.1450293971110256000775);
			       const __m256d v4_sc     = _mm256_set1_pd(0.8977727869612861128953);
			       const __m256d v4_beta   = _mm256_set1_pd(0.9628745099791263877469);
			       const __m256d v4_r1     = _mm256_set1_pd(1.3306700394914687909256);
			       const __m256d v4_r2     = _mm256_set1_pd(1.4142135623730950488017);
			       const __m256d v8_pi12   = _mm256_set1_pd(0.2617993877991494365386);
			       const __m256d v4_prek   = _mm256_set1_pd(1.6434564029725030125017);

                               __attribute__((always_inline))
			       __attribute__((aligned(32)))
			       static inline
			       __m256
			       ymm8r4_abs_ymm8r4(const __m256 v) {
                        	     const __m256 mask = _mm256_set1_ps(0x7FFFFFFF);
	                             return (_mm256_and_ps(x,mask));
                               }

			       __attribute__((always_inline))
			       __attribute__((aligned(32)))
			       static inline
			       __m256d
			       ymm4r8_abs_ymm4r8(const __m256d v) {
                        	     const __m256d mask = _mm256_set1_pd(0x7FFFFFFFFFFFFFFF);
	                             return (_mm256_and_pd(x,mask));
                               }

			       __attribute__((always_inline))
			       __attribute__((aligned(32)))
			       static inline
			       __m256
			       ymm8r4_sign_ymm8r4(const __m256 va,
			                          const __m256 vb) {
                                   register __m256 vret = v8_0;
				   register __m256 t0   = ymm8r4_abs_ymm8r4(va);
				   register __m256 mask = _mm256_cmp_ps(vb,v8_0,_CMP_GE_OQ);
				   vret = _mm256_blendv_ps(t0,_mm256_sub_ps(v8_0,t0),mask);
                                   return (vret);
			       }


			       __attribute__((always_inline))
			       __attribute__((aligned(32)))
			       static inline
			       __m256d
			       ymm4r8_sign_ymm4r8(const __m256d va,
			                          const __m256d vb) {
                                   register __m256d vret = v4_0;
				   const register __m256d t0   = ymm4r8_abs_ymm4r8(va);
				   const register __m256d mask = _mm256_cmp_pd(vb,v4_0,_CMP_GE_OQ);
				   vret = _mm256_blendv_pd(t0,_mm256_sub_pd(v4_0,t0),mask);
                                   return (vret);
			       }

			    
                               
			       __attribute__((always_inline))
			       __attribute__((aligned(32)))
			       static inline
			       __m256
			       ymm8r4_truncate(const __m256 v) {
                                  
                                    return (_mm256_cvtepi32_ps(_mm256_cvttps_epi32(v)));
			       }


			       __attribute__((always_inline))
			       __attribute__((aligned(32)))
			       static inline
			       __m256 fmod_ymm8r4(const __m256 a,
			                          const __m256 b) {
                                   const register __m256 v = ymm8r4_truncate(_mm256_div_ps(a,b));
                                   return (_mm256_sub_ps(a,_mm256_mul_ps(v,b)));  
			      }
                             
			       
                             /*  __attribute__((always_inline))
			       __attribute__((aligned(32)))
			      static inline
		              __m512 fmod_zmm16r4(const __m512 a,
		                                  const __m512 b) {

                                     __m512 v = _mm512_sub_ps(a,_mm512_mul_ps(
			             _mm512_div_round_ps(a,b,_MM_FROUND_TO_ZERO|_MM_FROUND_NO_EXEC),b));
			             return (v);
			  
		               }*/

		      
                            /*  __attribute__((always_inline))
			       __attribute__((aligned(32)))
		              static inline
		              __m512d fmod_zmm8r8(const __m512d a,
		                                  const __m512d b) {

                                    __m512d v = _mm512_sub_pd(a,_mm512_mul_pd(
			             _mm512_div_round_pd(a,b,_MM_FROUND_TO_ZERO|_MM_FROUND_NO_EXEC),b));
			       return (v);
			  
		            }*/


			    


			    __attribute__((always_inline))
			    __attribute__((aligned(32)))
		           static inline
		           __m256d norm2_ymm4r8(const __m256d y,
		                                const __m256d z,
					        const __m256d w) {

                                 const __m256d t0 = _mm256_mul_pd(y,y);
			         const __m256d t1 = _mm256_mul_pd(z,z);
			         const __m256d t2 = _mm256_mul_pd(w,w);
			         const __m256d v  = _mm256_add_pd(t0,_mm256_add_pd(t1,t2));
			         return (_mm256_sqrt_pd(v));
			    
		           }


		          __attribute__((always_inline))
			  __attribute__((aligned(32)))
		          static inline
		          __m256 norm2_ymm8r4(  const __m256 y,
		                                const __m256 z,
					        const __m256 w) {

                                const __m256 t0 = _mm256_mul_ps(y,y);
			        const __m256 t1 = _mm256_mul_ps(z,z);
			        const __m256 t2 = _mm256_mul_ps(w,w);
			        const __m256 v  = _mm256_add_ps(t0,_mm256_add_ps(t1,t2));
			        return (_mm256_sqrt_ps(v));
			    
		          }


			 __attribute__((always_inline))
			 __attribute__((aligned(32)))
		         static inline
		         __m256 clip_ymm8r4( const __m256 x,
		                             const __m256 lo,
					     const __m256 hi) {

                               return (_mm256_max_ps(lo,_mm256_min_ps(x,hi)));
		          }


		         __attribute__((always_inline))
			 __attribute__((aligned(32)))
		         static inline
		          __m256d clip_ymm4r8(const __m256d x,
		                              const __m256d lo,
					      const __m256d hi) {

                               return (_mm256_max_pd(lo,_mm256_min_pd(x,hi)));
		          }













#endif /*__GMS_ROTATIONS_AVX2_HELPERS_H__*/
