
#ifndef __GMS_AVX512_TRANSPOSITION_16X16_H__
#define __GMS_AVX512_TRANSPOSITION_16X16_H__

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
        Various transposition (load/store)  avx512 utiliy functions.
        Emulation of Matrix 16x16
*/

#include <immintrin.h>
#include "GMS_config"

namespace gms {

           namespace math {

                      // In-place
	              __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_zmm16r4_16x16_ip(__m512 &x00,
		                                      __m512 &x01,
						      __m512 &x02,
						      __m512 &x03,
						      __m512 &x04,
						      __m512 &x05,
						      __m512 &x06,
						      __m512 &x07,
						      __m512 &x08,
						      __m512 &x09,
						      __m512 &x10,
						      __m512 &x11,
						      __m512 &x12,
						      __m512 &x13,
						      __m512 &x14,
						      __m512 &x15); 


		      __ATTR_REGCALL__
                     __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     void transpose_u_zmm16r4_16x16_ip( float * __restrict x00,
		                                        float * __restrict x01,
							float * __restrict x02,
							float * __restrict x03,
							float * __restrict x04,
							float * __restrict x05,
							float * __restrict x06,
							float * __restrict x07,
							float * __restrict x08,
							float * __restrict x09,
							float * __restrict x10,
							float * __restrict x11,
							float * __restrict x12,
							float * __restrict x13,
							float * __restrict x14,
							float * __restrict x15); 

		      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_a_zmm16r4_16x16_ip(float * __restrict __ATTR_ALIGN__(64) x00,
		                                        float * __restrict __ATTR_ALIGN__(64) x01,
							float * __restrict __ATTR_ALIGN__(64) x02,
							float * __restrict __ATTR_ALIGN__(64) x03,
							float * __restrict __ATTR_ALIGN__(64) x04,
							float * __restrict __ATTR_ALIGN__(64) x05,
							float * __restrict __ATTR_ALIGN__(64) x06,
							float * __restrict __ATTR_ALIGN__(64) x07,
							float * __restrict __ATTR_ALIGN__(64) x08,
							float * __restrict __ATTR_ALIGN__(64) x09,
							float * __restrict __ATTR_ALIGN__(64) x10,
							float * __restrict __ATTR_ALIGN__(64) x11,
							float * __restrict __ATTR_ALIGN__(64) x12,
							float * __restrict __ATTR_ALIGN__(64) x13,
							float * __restrict __ATTR_ALIGN__(64) x14,
							float * __restrict __ATTR_ALIGN__(64) x15); 
			


		    
                     __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_u_zmm16r4_16x16_ip(float * __restrict x); 

		   
                    
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_a_zmm16r4_16x16_ip(float * __restrict __ATTR_ALIGN__(64) x); 
		      
		      
		    
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_u_zmm16r4_16x16(float * __restrict x,
		                                     float * __restrict y);

		     
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_a_zmm16r4_16x16(float * __restrict __ATTR_ALIGN__(64)  x,
		                                     float * __restrict __ATTR_ALIGN__(64)  y); 

#include <cstdint>


		    
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     void transpose_u_zmm16r4_16x16_v2(float * __restrict x,
		                                        float * __restrict y,
							const int32_t n); 


                    
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_a_zmm16r4_16x16_v2(float * __restrict __ATTR_ALIGN__(64))) x,
		                                        float * __restrict __ATTR_ALIGN__(64))) y,
							const int32_t n); 

		 


		    
               
     }

}











#endif /*__GMS_AVX512_TRANSPOSITION_16X16_H__*/
