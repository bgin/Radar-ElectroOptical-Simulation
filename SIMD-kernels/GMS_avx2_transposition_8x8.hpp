

#ifndef __GMS_AVX2_TRANSPOSITION_8X8_HPP__
#define __GMS_AVX2_TRANSPOSITION_8X8_HPP__


#include <immintrin.h>
#include "GMS_config.h"



namespace gms {


          namespace math {

	             // In-place
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_ymm8r4_8x8_ip(__m256 &x0,
		                                   __m256 &x1,
						   __m256 &x2,
						   __m256 &x3,
						   __m256 &x4,
						   __m256 &x5,
						   __m256 &x6,
						   __m256 &x7) {

                          register __m256 y0;
			  register __m256 y1;
			  register __m256 y2;
			  register __m256 y3;
			  register __m256 y4;
			  register __m256 y5;
			  register __m256 y6;
			  register __m256 y7;
			  register __m256 z0;
			  register __m256 z1;
			  register __m256 z2;
			  register __m256 z3;
			  register __m256 z4;
			  register __m256 z5;
			  register __m256 z6;
			  register __m256 z7;

			  y0 = _mm256_unpacklo_ps(x0,x1);
			  y1 = _mm256_unpackhi_ps(x0,x1);
			  y2 = _mm256_unpacklo_ps(x2,x3);
			  z0 = _mm256_shuffle_ps(y0,y2,_MM_SHUFFLE(1,0,1,0));
			  z1 = _mm256_shuffle_ps(y0,y2,_MM_SHUFFLE(3,2,3,2));
			  y3 = _mm256_unpackhi_ps(x2,x3);
			  z2 = _mm256_shuffle_ps(y1,y3,_MM_SHUFFLE(1,0,1,0));
			  z3 = _mm256_shuffle_ps(y1,y3,_MM_SHUFFLE(3,2,3,3));
			  y4 = _mm256_unpacklo_ps(x4,x5);
			  y5 = _mm256_unpackhi_ps(x4,x5);
			  y6 = _mm256_unpacklo_ps(x6,x7);
			  z4 = _mm256_shuffle_ps(y4,y6,_MM_SHUFFLE(1,0,1,0));
			  z5 = _mm256_shuffle_ps(y4,y6,_MM_SHUFFLE(3,2,3,2));
			  y7 = _mm256_unpackhi_ps(x6,x7);
			  z6 = _mm256_shuffle_ps(y5,y7,_MM_SHUFFLE(1,0,1,0));
			  z7 = _mm256_shuffle_ps(y5,y7,_MM_SHUFFLE(3,2,3,2));

			  x0 = _mm256_permute2f128_ps(z0,z4,0x20);
			  x1 = _mm256_permute2f128_ps(z1,z5,0x20);
			  x2 = _mm256_permute2f128_ps(z2,z6,0x20);
			  x3 = _mm256_permute2f128_ps(z3,z7,0x20);
			  x4 = _mm256_permute2f128_ps(z0,z4,0x31);
			  x5 = _mm256_permute2f128_ps(z1,z5,0x31);
			  x6 = _mm256_permute2f128_ps(z2,z6,0x31);
			  x7 = _mm256_permute2f128_ps(z3,z7,0x31);

		 }

     }


}









#endif /*__GMS_AVX2_TRANSPOSITION_8X8_HPP__*/
