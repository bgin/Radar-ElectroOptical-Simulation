



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


#include <immintrin.h>
#include "GMS_setv_avx512_unroll16x.h"



void gms::math::ssetv_u_zmm16r4_unroll16x(const int32_t n,
		                                  const float alpha,
						                  float * __restrict x,
						                  const int32_t incx) 
{

            if(__builtin_expect(0==n,0)) { return;}
			    __m512 alphav;
			    int32_t i;

			 if(__builtin_expect(1==incx,1)) 
			 {
                            
                    alphav = _mm512_set1_ps(alpha);
			        for(i = 0; (i+255) < n; i += 256) 
					{
			    
                        _mm512_storeu_ps(&x[i+0],alphav);
				        _mm512_storeu_ps(&x[i+16],alphav);
				        _mm512_storeu_ps(&x[i+32],alphav);
				        _mm512_storeu_ps(&x[i+48],alphav);
				        _mm512_storeu_ps(&x[i+64],alphav);
				        _mm512_storeu_ps(&x[i+80],alphav);
				        _mm512_storeu_ps(&x[i+96],alphav);
				        _mm512_storeu_ps(&x[i+112],alphav);
				        _mm512_storeu_ps(&x[i+128],alphav);
				        _mm512_storeu_ps(&x[i+144],alphav);
				        _mm512_storeu_ps(&x[i+160],alphav);
				        _mm512_storeu_ps(&x[i+176],alphav);
				        _mm512_storeu_ps(&x[i+192],alphav);
				        _mm512_storeu_ps(&x[i+208],alphav);
				        _mm512_storeu_ps(&x[i+224],alphav);
				        _mm512_storeu_ps(&x[i+240],alphav);
			    }

			    for(; (i+175) < n; i += 176) 
				{

			        _mm512_storeu_ps(&x[i+0],alphav);
				    _mm512_storeu_ps(&x[i+16],alphav);
				    _mm512_storeu_ps(&x[i+32],alphav);
				    _mm512_storeu_ps(&x[i+48],alphav);
				    _mm512_storeu_ps(&x[i+64],alphav);
				    _mm512_storeu_ps(&x[i+80],alphav);
				    _mm512_storeu_ps(&x[i+96],alphav);
				    _mm512_storeu_ps(&x[i+112],alphav);
				    _mm512_storeu_ps(&x[i+128],alphav);
				    _mm512_storeu_ps(&x[i+144],alphav);
				    _mm512_storeu_ps(&x[i+160],alphav);
			    }

			    for(; (i+127) < n; i += 128) 
				{

			        _mm512_storeu_ps(&x[i+0],alphav);
				    _mm512_storeu_ps(&x[i+16],alphav);
				    _mm512_storeu_ps(&x[i+32],alphav);
				    _mm512_storeu_ps(&x[i+48],alphav);
				    _mm512_storeu_ps(&x[i+64],alphav);
				    _mm512_storeu_ps(&x[i+80],alphav);
				    _mm512_storeu_ps(&x[i+96],alphav);
				    _mm512_storeu_ps(&x[i+112],alphav);
			    }

			    for(; (i+63) < n; i += 64) 
				{

                    _mm512_storeu_ps(&x[i+0],alphav);
				    _mm512_storeu_ps(&x[i+16],alphav);
				    _mm512_storeu_ps(&x[i+32],alphav);
				    _mm512_storeu_ps(&x[i+48],alphav);
			    }

			    for(; (i+31) < n; i += 32) 
				{

			        _mm512_storeu_ps(&x[i+0],alphav);
				    _mm512_storeu_ps(&x[i+16],alphav);
			    }

			    for(; (i+15) < n; i += 16) 
				{

                    _mm512_storeu_ps(&x[i+0],alphav);
			    }

			    for(; (i+0) < n; i += 1) 
				{

                      x[i] = alpha;
			    }
		      }
		      else 
			  {

                    for(i = 0; i != n; i += 1) 
					{

			           *x = alpha;
					    x  += incx;
				       
			        }
		      }
		}



		 
void gms::math::ssetv_a_zmm16r4_unroll16x(const int32_t n,
		                                  const float alpha,
						                  float * __restrict __ATTR_ALIGN__(64) x,
						                   const int32_t incx) 
{

            if(__builtin_expect(0==n,0)) { return;}
			   __m512 alphav;
			   int32_t i;

			if(__builtin_expect(1==incx,1)) 
		    {
                            
                 alphav = _mm512_set1_ps(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
                __assume_aligned(x,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                x = (float*)__builtin_assume_aligned(x,64);
#endif				    
			    for(i = 0; (i+255) < n; i += 256) 
				{
			    
                    _mm512_store_ps(&x[i+0],alphav);
				    _mm512_store_ps(&x[i+16],alphav);
				    _mm512_store_ps(&x[i+32],alphav);
				    _mm512_store_ps(&x[i+48],alphav);
				    _mm512_store_ps(&x[i+64],alphav);
				    _mm512_store_ps(&x[i+80],alphav);
				    _mm512_store_ps(&x[i+96],alphav);
				    _mm512_store_ps(&x[i+112],alphav);
				    _mm512_store_ps(&x[i+128],alphav);
				    _mm512_store_ps(&x[i+144],alphav);
				    _mm512_store_ps(&x[i+160],alphav);
				    _mm512_store_ps(&x[i+176],alphav);
				    _mm512_store_ps(&x[i+192],alphav);
				    _mm512_store_ps(&x[i+208],alphav);
				    _mm512_store_ps(&x[i+224],alphav);
				    _mm512_store_ps(&x[i+240],alphav);
			    }

			    for(; (i+175) < n; i += 176) 
				{

			        _mm512_store_ps(&x[i+0],alphav);
				    _mm512_store_ps(&x[i+16],alphav);
				    _mm512_store_ps(&x[i+32],alphav);
				    _mm512_store_ps(&x[i+48],alphav);
				    _mm512_store_ps(&x[i+64],alphav);
				    _mm512_store_ps(&x[i+80],alphav);
				    _mm512_store_ps(&x[i+96],alphav);
				    _mm512_store_ps(&x[i+112],alphav);
				    _mm512_store_ps(&x[i+128],alphav);
				    _mm512_store_ps(&x[i+144],alphav);
				    _mm512_store_ps(&x[i+160],alphav);
			    }

			    for(; (i+127) < n; i += 128) 
				{

			        _mm512_store_ps(&x[i+0],alphav);
				    _mm512_store_ps(&x[i+16],alphav);
				    _mm512_store_ps(&x[i+32],alphav);
				    _mm512_store_ps(&x[i+48],alphav);
				    _mm512_store_ps(&x[i+64],alphav);
				    _mm512_store_ps(&x[i+80],alphav);
				     _mm512_store_ps(&x[i+96],alphav);
				    _mm512_store_ps(&x[i+112],alphav);
			    }

			    for(; (i+63) < n; i += 64) 
				{

                    _mm512_store_ps(&x[i+0],alphav);
				    _mm512_store_ps(&x[i+16],alphav);
				    _mm512_store_ps(&x[i+32],alphav);
				    _mm512_store_ps(&x[i+48],alphav);
			    }

			    for(; (i+31) < n; i += 32) 
				{

			        _mm512_store_ps(&x[i+0],alphav);
				    _mm512_store_ps(&x[i+16],alphav);
			    }

			    for(; (i+15) < n; i += 16) 
				{

                    _mm512_store_ps(&x[i+0],alphav);
			    }

			    for(; (i+0) < n; i += 1) 
				{

                    x[i] = alpha;
			    }
		      }
		      else {

                    for(i = 0; i != n; i += 1) 
					{

			            *x = alpha;
						 x += incx;
				       
			        }
		      }
}



		 
void gms::math::dsetv_u_zmm8r8_unroll16x(const int32_t n,
		                                 const double alpha,
						                 double * __restrict x,
						                 const int32_t incx) {

        if(__builtin_expect(0==n,0)) { return;}
			    __m512d alphav;
			    int32_t i;

		if(__builtin_expect(1==incx,1)) 
		{
                            
                alphav = _mm512_set1_pd(alpha);

            for(i = 0; (i+127) < n; i += 128) 
			{
			    
                 _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
				 _mm512_storeu_pd(&x[i+32],alphav);
				 _mm512_storeu_pd(&x[i+40],alphav);
				 _mm512_storeu_pd(&x[i+48],alphav);
				 _mm512_storeu_pd(&x[i+56],alphav);
				 _mm512_storeu_pd(&x[i+64],alphav);
				 _mm512_storeu_pd(&x[i+72],alphav);
				 _mm512_storeu_pd(&x[i+80],alphav);
				 _mm512_storeu_pd(&x[i+88],alphav);
				 _mm512_storeu_pd(&x[i+96],alphav);
				 _mm512_storeu_pd(&x[i+104],alphav);
				 _mm512_storeu_pd(&x[i+112],alphav);
				 _mm512_storeu_pd(&x[i+120],alphav);
			}

			for(; (i+79) < n; i += 80) 
			{

			     _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
				 _mm512_storeu_pd(&x[i+32],alphav);
				 _mm512_storeu_pd(&x[i+40],alphav);
				 _mm512_storeu_pd(&x[i+48],alphav);
				 _mm512_storeu_pd(&x[i+56],alphav);
				 _mm512_storeu_pd(&x[i+64],alphav);
				 _mm512_storeu_pd(&x[i+72],alphav);
			}

			for(; (i+63) < n; i += 64) 
			{
 
                 _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
				 _mm512_storeu_pd(&x[i+32],alphav);
				 _mm512_storeu_pd(&x[i+40],alphav);
				 _mm512_storeu_pd(&x[i+48],alphav);
				 _mm512_storeu_pd(&x[i+56],alphav);
			}

			for(; (i+31) < n; i += 32) 
			{

                 _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
			}

			for(; (i+15) < n; i += 16) 
			{

			     _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
			}

			for(; (i+7) < n; i += 8) 
			{

                 _mm512_storeu_pd(&x[i+0],alphav);
			}

			for(; (i+0) < n; i += 1) 
			{

			          x[i] = alpha;
			}
		}
		else 
		{

			for(i = 0; i != n; i += 1) 
			{

			    *x  = alpha;
				 x += incx;
			   
			}
		}
}
   


void gms::math::dsetv_a_zmm8r8_unroll16x(const int32_t n,
		                                 const double alpha,
						                 double * __restrict __ATTR_ALIGN__(64) x,
						                 const int32_t incx) 
{

        if(__builtin_expect(0==n,0)) { return;}
			__m512d alphav;
			int32_t i;

		if(__builtin_expect(1==incx,1)) 
	    {
                            
            alphav = _mm512_set1_pd(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
            __assume_aligned(x,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
             x = (double*)__builtin_assume_aligned(x,64);
#endif
        for(i = 0; (i+127) < n; i += 128) 
		{
			    
                 _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
				 _mm512_storeu_pd(&x[i+32],alphav);
				 _mm512_storeu_pd(&x[i+40],alphav);
				 _mm512_storeu_pd(&x[i+48],alphav);
				 _mm512_storeu_pd(&x[i+56],alphav);
				 _mm512_storeu_pd(&x[i+64],alphav);
				 _mm512_storeu_pd(&x[i+72],alphav);
				 _mm512_storeu_pd(&x[i+80],alphav);
				 _mm512_storeu_pd(&x[i+88],alphav);
				 _mm512_storeu_pd(&x[i+96],alphav);
				 _mm512_storeu_pd(&x[i+104],alphav);
				 _mm512_storeu_pd(&x[i+112],alphav);
				 _mm512_storeu_pd(&x[i+120],alphav);
		}

		for(; (i+79) < n; i += 80) 
		{

			     _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
				 _mm512_storeu_pd(&x[i+32],alphav);
				 _mm512_storeu_pd(&x[i+40],alphav);
				 _mm512_storeu_pd(&x[i+48],alphav);
				 _mm512_storeu_pd(&x[i+56],alphav);
				 _mm512_storeu_pd(&x[i+64],alphav);
				 _mm512_storeu_pd(&x[i+72],alphav);
		}

		for(; (i+63) < n; i += 64) 
		{
 
                 _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
				 _mm512_storeu_pd(&x[i+32],alphav);
				 _mm512_storeu_pd(&x[i+40],alphav);
				 _mm512_storeu_pd(&x[i+48],alphav);
				 _mm512_storeu_pd(&x[i+56],alphav);
		}

		for(; (i+31) < n; i += 32) 
		{

                 _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
				 _mm512_storeu_pd(&x[i+16],alphav);
				 _mm512_storeu_pd(&x[i+24],alphav);
		}

		for(; (i+15) < n; i += 16) 
		{

			     _mm512_storeu_pd(&x[i+0],alphav);
				 _mm512_storeu_pd(&x[i+8],alphav);
		}

		for(;  (i+7) < n; i += 8) 
		{

                 _mm512_storeu_pd(&x[i+0],alphav);
		}

		for(; (i+0) < n; i += 1) 
		{

			    x[i] = alpha;
		}
	}
	else 
	{

		for(i = 0; i != n; i += 1)
		 {

			 *x  = alpha;
			  x += incx;
			  
		}
	}
}
   
 
