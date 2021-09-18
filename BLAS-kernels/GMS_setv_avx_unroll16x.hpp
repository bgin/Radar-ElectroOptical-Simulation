
#ifndef __GMS_SETV_AVX_UNROLL16X_HPP__
#define __GMS_SETV_AVX_UNROLL16X_HPP__ 180920211123


/*
     This version is based on BLIS framework original implementation.
     Few changes were introduced in order to more easily adapt it 
     to 'Guided Missile Simulation' project.
     
   The original authors license
   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.
   Copyright (C) 2016 - 2019, Advanced Micro Devices, Inc.
   Copyright (C) 2018 - 2020, The University of Texas at Austin. All rights reserved.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     
*/

namespace file_version {

    const unsigned int gGMS_SETV_AVX_UNROLL16X_MAJOR = 1U;
    const unsigned int gGMS_SETV_AVX_UNROLL16X_MINOR = 0U;
    const unsigned int gGMS_SETV_AVX_UNROLL16X_MICRO = 0U;
    const unsigned int gGMS_SETV_AVX_UNROLL16X_FULLVER =
      1000U*gGMS_SETV_AVX_UNROLL16X_MAJOR+
      100U*gGMS_SETV_AVX_UNROLL16X_MINOR+
      10U*gGMS_SETVV_AVX_UNROLL16X_MICRO;
    const char * const pgGMS_SETV_AVX_UNROLL16X_CREATION_DATE = "18-09-2021 11:23 AM +00200 (SAT 18 SEP 2021 GMT+2)";
    const char * const pgGMS_SETV_AVX_UNROLL16X_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_SETV_AVX_UNROLL16X_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_SETV_AVX_UNROLL16X_DESCRIPTION   = "AVX optimized SETV kernels."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace  gms {

          namespace math {


	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void ssetv_u_ymm8r4_unroll16x(const int32_t n,
		                                 const float alpha,
						 float * __restrict x,
						 const int32_t incx) {

		         if(__builtin_expect(0==n,0)) { return;}
                         __m256 alphav;
			 int32_t i;
			 
			 if(__builtin_expect(1==incx,1)) {
			 
                            alphav = _mm256_broadcast_ss(alpha);
			    for(i = 0; (i+127) < n; i += 128) {
			    
                                 _mm256_storeu_ps(&x[i+0],alphav);
				 _mm256_storeu_ps(&x[i+8],alphav);
				 _mm256_storeu_ps(&x[i+16],alphav);
				 _mm256_storeu_ps(&x[i+24],alphav);
				 _mm256_storeu_ps(&x[i+32],alphav);
				 _mm256_storeu_ps(&x[i+40],alphav);
				 _mm256_storeu_ps(&x[i+48],alphav);
				 _mm256_storeu_ps(&x[i+56],alphav);
				 _mm256_storeu_ps(&x[i+64],alphav);
				 _mm256_storeu_ps(&x[i+72],alphav);
				 _mm256_storeu_ps(&x[i+80],alphav);
				 _mm256_storeu_ps(&x[i+88],alphav);
				 _mm256_storeu_ps(&x[i+96],alphav);
				 _mm256_storeu_ps(&x[i+104],alphav);
				 _mm256_storeu_ps(&x[i+112],alphav);
				 _mm256_storeu_ps(&x[i+120],alphav);
			    }

			    for(; (i+79) < n; i += 80) {

			         _mm256_storeu_ps(&x[i+0],alphav);
				 _mm256_storeu_ps(&x[i+8],alphav);
				 _mm256_storeu_ps(&x[i+16],alphav);
				 _mm256_storeu_ps(&x[i+24],alphav);
				 _mm256_storeu_ps(&x[i+32],alphav);
				 _mm256_storeu_ps(&x[i+40],alphav);
				 _mm256_storeu_ps(&x[i+48],alphav);
				 _mm256_storeu_ps(&x[i+56],alphav);
				 _mm256_storeu_ps(&x[i+64],alphav);
				 _mm256_storeu_ps(&x[i+72],alphav);
			    }

			    for(; (i+63) < n; i += 64) {
 
                                 _mm256_storeu_ps(&x[i+0],alphav);
				 _mm256_storeu_ps(&x[i+8],alphav);
				 _mm256_storeu_ps(&x[i+16],alphav);
				 _mm256_storeu_ps(&x[i+24],alphav);
				 _mm256_storeu_ps(&x[i+32],alphav);
				 _mm256_storeu_ps(&x[i+40],alphav);
				 _mm256_storeu_ps(&x[i+48],alphav);
				 _mm256_storeu_ps(&x[i+56],alphav);
			    }

			    for(; (i+31) < n; i += 32) {

                                 _mm256_storeu_ps(&x[i+0],alphav);
				 _mm256_storeu_ps(&x[i+8],alphav);
				 _mm256_storeu_ps(&x[i+16],alphav);
				 _mm256_storeu_ps(&x[i+24],alphav);
			    }

			    for(; (i+15) < n; i += 16) {

			         _mm256_storeu_ps(&x[i+0],alphav);
				 _mm256_storeu_ps(&x[i+8],alphav);
			    }

			    for(; i (i+7) < n; i += 8) {

                                   _mm256_storeu_ps(&x[i+0],alphav);
			    }

			    for(; (i+0) < n; i += 1) {

			          x[i] = alpha;
			    }
		      }
		      else {

                              for(i = 0; i != n; ++i) {

			           *x = alpha;
				    x += incx;
			      }
		      }
                      
	       }



	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void ssetv_a_ymm8r4_unroll16x(const int32_t n,
		                                 const float alpha,
						 float * __restrict __ATTR_ALIGN__(32) x,
						 const int32_t incx) {

		         if(__builtin_expect(0==n,0)) { return;}
                         __m256 alphav;
			 int32_t i;
			 
			 if(__builtin_expect(1==incx,1)) {
			 
                            alphav = _mm256_broadcast_ss(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (float*)__builtin_assume_aligned(x,32);
#endif			    
			    for(i = 0; (i+127) < n; i += 128) {
			    
                                 _mm256_store_ps(&x[i+0],alphav);
				 _mm256_store_ps(&x[i+8],alphav);
				 _mm256_store_ps(&x[i+16],alphav);
				 _mm256_store_ps(&x[i+24],alphav);
				 _mm256_store_ps(&x[i+32],alphav);
				 _mm256_store_ps(&x[i+40],alphav);
				 _mm256_store_ps(&x[i+48],alphav);
				 _mm256_store_ps(&x[i+56],alphav);
				 _mm256_store_ps(&x[i+64],alphav);
				 _mm256_store_ps(&x[i+72],alphav);
				 _mm256_store_ps(&x[i+80],alphav);
				 _mm256_store_ps(&x[i+88],alphav);
				 _mm256_store_ps(&x[i+96],alphav);
				 _mm256_store_ps(&x[i+104],alphav);
				 _mm256_store_ps(&x[i+112],alphav);
				 _mm256_store_ps(&x[i+120],alphav);
			    }

			    for(; (i+79) < n; i += 80) {

			         _mm256_store_ps(&x[i+0],alphav);
				 _mm256_store_ps(&x[i+8],alphav);
				 _mm256_store_ps(&x[i+16],alphav);
				 _mm256_store_ps(&x[i+24],alphav);
				 _mm256_store_ps(&x[i+32],alphav);
				 _mm256_store_ps(&x[i+40],alphav);
				 _mm256_store_ps(&x[i+48],alphav);
				 _mm256_store_ps(&x[i+56],alphav);
				 _mm256_store_ps(&x[i+64],alphav);
				 _mm256_store_ps(&x[i+72],alphav);
			    }

			    for(; (i+63) < n; i += 64) {
 
                                 _mm256_store_ps(&x[i+0],alphav);
				 _mm256_store_ps(&x[i+8],alphav);
				 _mm256_store_ps(&x[i+16],alphav);
				 _mm256_store_ps(&x[i+24],alphav);
				 _mm256_store_ps(&x[i+32],alphav);
				 _mm256_store_ps(&x[i+40],alphav);
				 _mm256_store_ps(&x[i+48],alphav);
				 _mm256_store_ps(&x[i+56],alphav);
			    }

			    for(; (i+31) < n; i += 32) {

                                 _mm256_store_ps(&x[i+0],alphav);
				 _mm256_store_ps(&x[i+8],alphav);
				 _mm256_store_ps(&x[i+16],alphav);
				 _mm256_store_ps(&x[i+24],alphav);
			    }

			    for(; (i+15) < n; i += 16) {

			         _mm256_store_ps(&x[i+0],alphav);
				 _mm256_store_ps(&x[i+8],alphav);
			    }

			    for(; i (i+7) < n; i += 8) {

                                   _mm256_store_ps(&x[i+0],alphav);
			    }

			    for(; (i+0) < n; i += 1) {

			          x[i] = alpha;
			    }
		      }
		      else {

                              for(i = 0; i != n; ++i) {

			           *x = alpha;
				    x += incx;
			      }
		      }
                      
	       }



	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void dsetv_u_ymm4r8_unroll16x(const int32_t n,
		                                 const double alpha,
						 double * __restrict x,
						 const int32_t incx) {

                         if(__builtin_expect(0==n,0)) { return;}
                         __m256d alphav;
			 int32_t i;
			 
			 if(__builtin_expect(1==incx,1)) {

                             alphav = _mm256_broadcast_sd(alpha);

			     for(i = 0; (i+63) < n; i += 64) {

                                 _mm256_storeu_pd(&x[i+0],alphav);
				 _mm256_storeu_pd(&x[i+4],alphav);
				 _mm256_storeu_pd(&x[i+8],alphav);
				 _mm256_storeu_pd(&x[i+12],alphav);
				 _mm256_storeu_pd(&x[i+16],alphav);
				 _mm256_storeu_pd(&x[i+20],alphav);
				 _mm256_storeu_pd(&x[i+24],alphav);
				 _mm256_storeu_pd(&x[i+28],alphav);
				 _mm256_storeu_pd(&x[i+32],alphav);
				 _mm256_storeu_pd(&x[i+36],alphav);
				 _mm256_storeu_pd(&x[i+40],alphav);
				 _mm256_storeu_pd(&x[i+44],alphav);
				 _mm256_storeu_pd(&x[i+48],alphav);
				 _mm256_storeu_pd(&x[i+52],alphav);
				 _mm256_storeu_pd(&x[i+56],alphav);
				 _mm256_storeu_pd(&x[i+60],alphav);
			     }

			     for(; (i+39) < n; i += 40) {

			          _mm256_storeu_pd(&x[i+0],alphav);
				  _mm256_storeu_pd(&x[i+4],alphav);
				  _mm256_storeu_pd(&x[i+8],alphav);
				  _mm256_storeu_pd(&x[i+12],alphav);
				  _mm256_storeu_pd(&x[i+16],alphav);
				  _mm256_storeu_pd(&x[i+20],alphav);
				  _mm256_storeu_pd(&x[i+24],alphav);
				  _mm256_storeu_pd(&x[i+28],alphav);
				  _mm256_storeu_pd(&x[i+32],alphav);
				  _mm256_storeu_pd(&x[i+36],alphav);
			     }

			     for(; (i+31) < n; i += 32) {

                                    _mm256_storeu_pd(&x[i+0],alphav);
				    _mm256_storeu_pd(&x[i+4],alphav);
				    _mm256_storeu_pd(&x[i+8],alphav);
				    _mm256_storeu_pd(&x[i+12],alphav);
				    _mm256_storeu_pd(&x[i+16],alphav);
				    _mm256_storeu_pd(&x[i+20],alphav);
				    _mm256_storeu_pd(&x[i+24],alphav);
				    _mm256_storeu_pd(&x[i+28],alphav);
			     }

			     for(; (i+15) < n; i += 16) {

                                    _mm256_storeu_pd(&x[i+0],alphav);
				    _mm256_storeu_pd(&x[i+4],alphav);
				    _mm256_storeu_pd(&x[i+8],alphav);
				    _mm256_storeu_pd(&x[i+12],alphav);
			     }

			     for(; (i+7) < n; i += 8) {

                                   _mm256_storeu_pd(&x[i+0],alphav);
				   _mm256_storeu_pd(&x[i+4],alphav);
			     }

			     for(; (i+3) < n; i += 4) {

			           _mm256_storeu_pd(&x[i+0],alphav);
			     }

			     for(; (i+0) < n; i += 1) {

			           x[i] = alpha;
			     }
			     
			 }
			 else {

			       for(i = 0; i != n; ++i) {
                                    *x = alpha;
				     x += incx;
			       }
			 }
		  }



		   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void dsetv_a_ymm4r8_unroll16x(const int32_t n,
		                                 const double alpha,
						 double * __restrict __ATTR_ALIGN__(32) x,
						 const int32_t incx) {

                         if(__builtin_expect(0==n,0)) { return;}
                         __m256d alphav;
			 int32_t i;
			 
			 if(__builtin_expect(1==incx,1)) {

                             alphav = _mm256_broadcast_sd(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (double*)__builtin_assume_aligned(x,32);
#endif
			     for(i = 0; (i+63) < n; i += 64) {

                                 _mm256_store_pd(&x[i+0],alphav);
				 _mm256_store_pd(&x[i+4],alphav);
				 _mm256_store_pd(&x[i+8],alphav);
				 _mm256_store_pd(&x[i+12],alphav);
				 _mm256_store_pd(&x[i+16],alphav);
				 _mm256_store_pd(&x[i+20],alphav);
				 _mm256_store_pd(&x[i+24],alphav);
				 _mm256_store_pd(&x[i+28],alphav);
				 _mm256_store_pd(&x[i+32],alphav);
				 _mm256_store_pd(&x[i+36],alphav);
				 _mm256_store_pd(&x[i+40],alphav);
				 _mm256_store_pd(&x[i+44],alphav);
				 _mm256_store_pd(&x[i+48],alphav);
				 _mm256_store_pd(&x[i+52],alphav);
				 _mm256_store_pd(&x[i+56],alphav);
				 _mm256_store_pd(&x[i+60],alphav);
			     }

			     for(; (i+39) < n; i += 40) {

			          _mm256_store_pd(&x[i+0],alphav);
				  _mm256_store_pd(&x[i+4],alphav);
				  _mm256_store_pd(&x[i+8],alphav);
				  _mm256_store_pd(&x[i+12],alphav);
				  _mm256_store_pd(&x[i+16],alphav);
				  _mm256_store_pd(&x[i+20],alphav);
				  _mm256_store_pd(&x[i+24],alphav);
				  _mm256_store_pd(&x[i+28],alphav);
				  _mm256_store_pd(&x[i+32],alphav);
				  _mm256_store_pd(&x[i+36],alphav);
			     }

			     for(; (i+31) < n; i += 32) {

                                    _mm256_store_pd(&x[i+0],alphav);
				    _mm256_store_pd(&x[i+4],alphav);
				    _mm256_store_pd(&x[i+8],alphav);
				    _mm256_store_pd(&x[i+12],alphav);
				    _mm256_store_pd(&x[i+16],alphav);
				    _mm256_store_pd(&x[i+20],alphav);
				    _mm256_store_pd(&x[i+24],alphav);
				    _mm256_store_pd(&x[i+28],alphav);
			     }

			     for(; (i+15) < n; i += 16) {

                                    _mm256_store_pd(&x[i+0],alphav);
				    _mm256_store_pd(&x[i+4],alphav);
				    _mm256_store_pd(&x[i+8],alphav);
				    _mm256_store_pd(&x[i+12],alphav);
			     }

			     for(; (i+7) < n; i += 8) {

                                   _mm256_store_pd(&x[i+0],alphav);
				   _mm256_store_pd(&x[i+4],alphav);
			     }

			     for(; (i+3) < n; i += 4) {

			           _mm256_store_pd(&x[i+0],alphav);
			     }

			     for(; (i+0) < n; i += 1) {

			           x[i] = alpha;
			     }
			     
			 }
			 else {

			       for(i = 0; i != n; ++i) {
                                    *x = alpha;
				     x += incx;
			       }
			 }
		  }


                  
    }

}















#endif /*__GMS_SETV_AVX_UNROLL16X_HPP__*/
