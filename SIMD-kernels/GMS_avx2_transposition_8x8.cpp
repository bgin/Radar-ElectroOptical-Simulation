






#include "GMS_avx2_transposition_8x8.h"




		      void gms::math::transpose_ymm8r4_8x8_ip(__m256 &x0,
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


	               // In-place
	            
		      void gms::math::transpose_u_ymm8r4_8x8_ip(float * __restrict x0,
		                                     float * __restrict x1,
						     float * __restrict x2,
						     float * __restrict x3,
						     float * __restrict x4,
						     float * __restrict x5,
						     float * __restrict x6,
						     float * __restrict x7) {

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

			  y0 = _mm256_loadu_ps(x0);
			  y1 = _mm256_loadu_ps(x1);
			  z0 = _mm256_unpacklo_ps(y0,y1);
			  z1 = _mm256_unpackhi_ps(y0,y1);
			  y2 = _mm256_loadu_ps(x2);
			  y3 = _mm256_loadu_ps(x3);
			  z2 = _mm256_unpacklo_ps(y2,y3);
			  z3 = _mm256_unpackhi_ps(y2,y3);
			  y4 = _mm256_loadu_ps(x4);
			  y5 = _mm256_loadu_ps(x5);
			  z4 = _mm256_unpacklo_ps(x4,x5);
			  z5 = _mm256_unpackhi_ps(x4,x5);
			  y6 = _mm256_loadu_ps(x6);
			  y7 = _mm256_loadu_ps(x7);
			  z6 = _mm256_unpacklo_ps(x6,x7);
			  z7 = _mm256_unpackhi_ps(x6,x7);
                          y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			  y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			  y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			  y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			  y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			  _mm256_storeu_ps(x0,_mm256_permute2f128_ps(y0,y4,0x20));
			  _mm256_storeu_ps(x4,_mm256_permute2f128_ps(y0,y4,0x31));
			  y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			  _mm256_storeu_ps(x1,_mm256_permute2f128_ps(y1,y5,0x20));
			  _mm256_storeu_ps(x5,_mm256_permute2f128_ps(y1,y5,0x31));
			  y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			  _mm256_storeu_ps(x2,_mm256_permute2f128_ps(y2,y6,0x20));
			  _mm256_storeu_ps(x6,_mm256_permute2f128_ps(y2,y6,0x31));
			  y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			  _mm256_storeu_ps(x3,_mm256_permute2f128_ps(y3,y7,0x20));
			  _mm256_storeu_ps(x7,_mm256_permute2f128_ps(y3,y7,0x31));
		 }


		    // In-place
	            
		      void gms::math::transpose_a_ymm8r4_8x8_ip(float * __restrict __ATTR_ALIGN__(32) x0,
		                                     float * __restrict __ATTR_ALIGN__(32) x1,
						     float * __restrict __ATTR_ALIGN__(32) x2,
						     float * __restrict __ATTR_ALIGN__(32) x3,
						     float * __restrict __ATTR_ALIGN__(32) x4,
						     float * __restrict __ATTR_ALIGN__(32) x5,
						     float * __restrict __ATTR_ALIGN__(32) x6,
						     float * __restrict __ATTR_ALIGN__(32) x7) {

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

			  y0 = _mm256_load_ps(x0);
			  y1 = _mm256_load_ps(x1);
			  z0 = _mm256_unpacklo_ps(y0,y1);
			  z1 = _mm256_unpackhi_ps(y0,y1);
			  y2 = _mm256_load_ps(x2);
			  y3 = _mm256_load_ps(x3);
			  z2 = _mm256_unpacklo_ps(y2,y3);
			  z3 = _mm256_unpackhi_ps(y2,y3);
			  y4 = _mm256_load_ps(x4);
			  y5 = _mm256_load_ps(x5);
			  z4 = _mm256_unpacklo_ps(x4,x5);
			  z5 = _mm256_unpackhi_ps(x4,x5);
			  y6 = _mm256_load_ps(x6);
			  y7 = _mm256_load_ps(x7);
			  z6 = _mm256_unpacklo_ps(x6,x7);
			  z7 = _mm256_unpackhi_ps(x6,x7);
                          y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			  y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			  y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			  y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			  y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			  _mm256_store_ps(x0,_mm256_permute2f128_ps(y0,y4,0x20));
			  _mm256_store_ps(x4,_mm256_permute2f128_ps(y0,y4,0x31));
			  y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			  _mm256_store_ps(x1,_mm256_permute2f128_ps(y1,y5,0x20));
			  _mm256_store_ps(x5,_mm256_permute2f128_ps(y1,y5,0x31));
			  y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			  _mm256_store_ps(x2,_mm256_permute2f128_ps(y2,y6,0x20));
			  _mm256_store_ps(x6,_mm256_permute2f128_ps(y2,y6,0x31));
			  y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			  _mm256_store_ps(x3,_mm256_permute2f128_ps(y3,y7,0x20));
			  _mm256_store_ps(x7,_mm256_permute2f128_ps(y3,y7,0x31));
		 }



		        // In-place
	            
		      void gms::math::transpose_u_ymm8r4_8x8_ip(float * __restrict x) {

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

			  y0 = _mm256_loadu_ps(&x[0*8]);
			  y1 = _mm256_loadu_ps(&x[1*8]);
			  z0 = _mm256_unpacklo_ps(y0,y1);
			  z1 = _mm256_unpackhi_ps(y0,y1);
			  y2 = _mm256_loadu_ps(&x[2*8]);
			  y3 = _mm256_loadu_ps(&x[3*8]);
			  z2 = _mm256_unpacklo_ps(y2,y3);
			  z3 = _mm256_unpackhi_ps(y2,y3);
			  y4 = _mm256_loadu_ps(&x[4*8]);
			  y5 = _mm256_loadu_ps(&x[5*8]);
			  z4 = _mm256_unpacklo_ps(x4,x5);
			  z5 = _mm256_unpackhi_ps(x4,x5);
			  y6 = _mm256_loadu_ps(&x[6*8]);
			  y7 = _mm256_loadu_ps(&x[7*8]);
			  z6 = _mm256_unpacklo_ps(x6,x7);
			  z7 = _mm256_unpackhi_ps(x6,x7);
                          y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			  y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			  y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			  y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			  y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			  _mm256_storeu_ps(&x[0*8],_mm256_permute2f128_ps(y0,y4,0x20));
			  _mm256_storeu_ps(&x[4*8],_mm256_permute2f128_ps(y0,y4,0x31));
			  y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			  _mm256_storeu_ps(&x[1*8],_mm256_permute2f128_ps(y1,y5,0x20));
			  _mm256_storeu_ps(&x[5*8],_mm256_permute2f128_ps(y1,y5,0x31));
			  y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			  _mm256_storeu_ps(&x[2*8],_mm256_permute2f128_ps(y2,y6,0x20));
			  _mm256_storeu_ps(&x[6*8],_mm256_permute2f128_ps(y2,y6,0x31));
			  y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			  _mm256_storeu_ps(&x[3*8],_mm256_permute2f128_ps(y3,y7,0x20));
			  _mm256_storeu_ps(&x[7*8],_mm256_permute2f128_ps(y3,y7,0x31));
		  }



		          // In-place
	             
		      void gms::math::transpose_u_ymm8r4_8x8_ip(float * __restrict __ATTR_ALIGN__(32) x) {

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

			  y0 = _mm256_load_ps(&x[0*8]);
			  y1 = _mm256_load_ps(&x[1*8]);
			  z0 = _mm256_unpacklo_ps(y0,y1);
			  z1 = _mm256_unpackhi_ps(y0,y1);
			  y2 = _mm256_load_ps(&x[2*8]);
			  y3 = _mm256_load_ps(&x[3*8]);
			  z2 = _mm256_unpacklo_ps(y2,y3);
			  z3 = _mm256_unpackhi_ps(y2,y3);
			  y4 = _mm256_load_ps(&x[4*8]);
			  y5 = _mm256_load_ps(&x[5*8]);
			  z4 = _mm256_unpacklo_ps(x4,x5);
			  z5 = _mm256_unpackhi_ps(x4,x5);
			  y6 = _mm256_load_ps(&x[6*8]);
			  y7 = _mm256_load_ps(&x[7*8]);
			  z6 = _mm256_unpacklo_ps(x6,x7);
			  z7 = _mm256_unpackhi_ps(x6,x7);
                          y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			  y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			  y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			  y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			  y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			  _mm256_store_ps(&x[0*8],_mm256_permute2f128_ps(y0,y4,0x20));
			  _mm256_store_ps(&x[4*8],_mm256_permute2f128_ps(y0,y4,0x31));
			  y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			  _mm256_store_ps(&x[1*8],_mm256_permute2f128_ps(y1,y5,0x20));
			  _mm256_store_ps(&x[5*8],_mm256_permute2f128_ps(y1,y5,0x31));
			  y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			  _mm256_store_ps(&x[2*8],_mm256_permute2f128_ps(y2,y6,0x20));
			  _mm256_store_ps(&x[6*8],_mm256_permute2f128_ps(y2,y6,0x31));
			  y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			  _mm256_store_ps(&x[3*8],_mm256_permute2f128_ps(y3,y7,0x20));
			  _mm256_store_ps(&x[7*8],_mm256_permute2f128_ps(y3,y7,0x31));
		  }



		    
		      void gms::math::transpose_u_ymm8r4_8x8(float * __restrict x,
		                                  float * __restrict y) {

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

			  y0 = _mm256_loadu_ps(&x[0*8]);
			  y1 = _mm256_loadu_ps(&x[1*8]);
			  z0 = _mm256_unpacklo_ps(y0,y1);
			  z1 = _mm256_unpackhi_ps(y0,y1);
			  y2 = _mm256_loadu_ps(&x[2*8]);
			  y3 = _mm256_loadu_ps(&x[3*8]);
			  z2 = _mm256_unpacklo_ps(y2,y3);
			  z3 = _mm256_unpackhi_ps(y2,y3);
			  y4 = _mm256_loadu_ps(&x[4*8]);
			  y5 = _mm256_loadu_ps(&x[5*8]);
			  z4 = _mm256_unpacklo_ps(x4,x5);
			  z5 = _mm256_unpackhi_ps(x4,x5);
			  y6 = _mm256_loadu_ps(&x[6*8]);
			  y7 = _mm256_loadu_ps(&x[7*8]);
			  z6 = _mm256_unpacklo_ps(x6,x7);
			  z7 = _mm256_unpackhi_ps(x6,x7);
                          y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			  y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			  y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			  y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			  y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			  _mm256_storeu_ps(&y[0*8],_mm256_permute2f128_ps(y0,y4,0x20));
			  _mm256_storeu_ps(&y[4*8],_mm256_permute2f128_ps(y0,y4,0x31));
			  y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			  _mm256_storeu_ps(&y[1*8],_mm256_permute2f128_ps(y1,y5,0x20));
			  _mm256_storeu_ps(&y[5*8],_mm256_permute2f128_ps(y1,y5,0x31));
			  y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			  _mm256_storeu_ps(&y[2*8],_mm256_permute2f128_ps(y2,y6,0x20));
			  _mm256_storeu_ps(&y[6*8],_mm256_permute2f128_ps(y2,y6,0x31));
			  y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			  _mm256_storeu_ps(&y[3*8],_mm256_permute2f128_ps(y3,y7,0x20));
			  _mm256_storeu_ps(&y[7*8],_mm256_permute2f128_ps(y3,y7,0x31));
		  }



		    
		      void gms::math::transpose_a_ymm8r4_8x8(float * __restrict __ATTR_ALIGN__(32) x,
		                                  float * __restrict __ATTR_ALIGN__(32) y) {

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

			  y0 = _mm256_load_ps(&x[0*8]);
			  y1 = _mm256_load_ps(&x[1*8]);
			  z0 = _mm256_unpacklo_ps(y0,y1);
			  z1 = _mm256_unpackhi_ps(y0,y1);
			  y2 = _mm256_load_ps(&x[2*8]);
			  y3 = _mm256_load_ps(&x[3*8]);
			  z2 = _mm256_unpacklo_ps(y2,y3);
			  z3 = _mm256_unpackhi_ps(y2,y3);
			  y4 = _mm256_load_ps(&x[4*8]);
			  y5 = _mm256_load_ps(&x[5*8]);
			  z4 = _mm256_unpacklo_ps(x4,x5);
			  z5 = _mm256_unpackhi_ps(x4,x5);
			  y6 = _mm256_load_ps(&x[6*8]);
			  y7 = _mm256_load_ps(&x[7*8]);
			  z6 = _mm256_unpacklo_ps(x6,x7);
			  z7 = _mm256_unpackhi_ps(x6,x7);
                          y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			  y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			  y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			  y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			  y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			  _mm256_store_ps(&y[0*8],_mm256_permute2f128_ps(y0,y4,0x20));
			  _mm256_store_ps(&y[4*8],_mm256_permute2f128_ps(y0,y4,0x31));
			  y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			  _mm256_store_ps(&y[1*8],_mm256_permute2f128_ps(y1,y5,0x20));
			  _mm256_store_ps(&y[5*8],_mm256_permute2f128_ps(y1,y5,0x31));
			  y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			  _mm256_store_ps(&y[2*8],_mm256_permute2f128_ps(y2,y6,0x20));
			  _mm256_store_ps(&y[6*8],_mm256_permute2f128_ps(y2,y6,0x31));
			  y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			  _mm256_store_ps(&y[3*8],_mm256_permute2f128_ps(y3,y7,0x20));
			  _mm256_store_ps(&y[7*8],_mm256_permute2f128_ps(y3,y7,0x31));
		  }

#include <cstdint>


                     
		      void gms::math::transpose_u_ymm8r4_8x8_v2(float * __restrict x,
		                                     float * __restrict y,
						     const int32_t n) {
                          if(__builtin_expect((n%8)!=0,1)) {return;}
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
                          constexpr int32_t stride = 64;
			  int32_t i,j;
			  j = 0;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif   
                         for(i = 0; i != n; ++i) {   			  
			     y0 = _mm256_loadu_ps(&x[0*8+j]);
			     y1 = _mm256_loadu_ps(&x[1*8+j]);
			     z0 = _mm256_unpacklo_ps(y0,y1);
			     z1 = _mm256_unpackhi_ps(y0,y1);
			     y2 = _mm256_loadu_ps(&x[2*8+j]);
			     y3 = _mm256_loadu_ps(&x[3*8+j]);
			     z2 = _mm256_unpacklo_ps(y2,y3);
			     z3 = _mm256_unpackhi_ps(y2,y3);
			     y4 = _mm256_loadu_ps(&x[4*8+j]);
			     y5 = _mm256_loadu_ps(&x[5*8+j]);
			     z4 = _mm256_unpacklo_ps(x4,x5);
			     z5 = _mm256_unpackhi_ps(x4,x5);
			     y6 = _mm256_loadu_ps(&x[6*8+j]);
			     y7 = _mm256_loadu_ps(&x[7*8+j]);
			     z6 = _mm256_unpacklo_ps(x6,x7);
			     z7 = _mm256_unpackhi_ps(x6,x7);
                             y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			     y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			     y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			     y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			     y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			     _mm256_storeu_ps(&y[0*8+j],_mm256_permute2f128_ps(y0,y4,0x20));
			     _mm256_storeu_ps(&y[4*8+j],_mm256_permute2f128_ps(y0,y4,0x31));
			     y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			     _mm256_storeu_ps(&y[1*8+j],_mm256_permute2f128_ps(y1,y5,0x20));
			     _mm256_storeu_ps(&y[5*8+j],_mm256_permute2f128_ps(y1,y5,0x31));
			     y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			     _mm256_storeu_ps(&y[2*8+j],_mm256_permute2f128_ps(y2,y6,0x20));
			     _mm256_storeu_ps(&y[6*8+j],_mm256_permute2f128_ps(y2,y6,0x31));
			     y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			     _mm256_storeu_ps(&y[3*8+j],_mm256_permute2f128_ps(y3,y7,0x20));
			     _mm256_storeu_ps(&y[7*8+j],_mm256_permute2f128_ps(y3,y7,0x31));
                             j += stride;
			 }
		  }




		     
		      void gms::math::transpose_a_ymm8r4_8x8_v2(float * __restrict __ATTR_ALIGN__(32) x,
		                                     float * __restrict __ATTR_ALIGN__(32) y,
						     const int32_t n) {
                          if(__builtin_expect((n%8)!=0,1)) {return;}
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
                          constexpr int32_t stride = 64;
			  int32_t i,j;
			  j = 0;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                           __assume_aligned(x,32);
			   __assume_aligned(y,32);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                           x = (float*)__builtin_assume_aligned(x,32);
			   y = (float*)__nuiltin_assume_aligned(y,32);
#endif			  
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif   
                         for(i = 0; i != n; ++i) {   			  
			     y0 = _mm256_load_ps(&x[0*8+j]);
			     y1 = _mm256_load_ps(&x[1*8+j]);
			     z0 = _mm256_unpacklo_ps(y0,y1);
			     z1 = _mm256_unpackhi_ps(y0,y1);
			     y2 = _mm256_load_ps(&x[2*8+j]);
			     y3 = _mm256_load_ps(&x[3*8+j]);
			     z2 = _mm256_unpacklo_ps(y2,y3);
			     z3 = _mm256_unpackhi_ps(y2,y3);
			     y4 = _mm256_load_ps(&x[4*8+j]);
			     y5 = _mm256_load_ps(&x[5*8+j]);
			     z4 = _mm256_unpacklo_ps(x4,x5);
			     z5 = _mm256_unpackhi_ps(x4,x5);
			     y6 = _mm256_load_ps(&x[6*8+j]);
			     y7 = _mm256_load_ps(&x[7*8+j]);
			     z6 = _mm256_unpacklo_ps(x6,x7);
			     z7 = _mm256_unpackhi_ps(x6,x7);
                             y0 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(1,0,1,0));
			     y1 = _mm256_shuffle_ps(z0,z2,_MM_SHUFFLE(3,2,3,2));
			     y2 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(1,0,1,0));
			     y3 = _mm256_shuffle_ps(z1,z3,_MM_SHUFFLE(3,2,3,2));
			     y4 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(1,0,1,0));
			     _mm256_store_ps(&y[0*8+j],_mm256_permute2f128_ps(y0,y4,0x20));
			     _mm256_store_ps(&y[4*8+j],_mm256_permute2f128_ps(y0,y4,0x31));
			     y5 = _mm256_shuffle_ps(z4,z6,_MM_SHUFFLE(3,2,3,2));
			     _mm256_store_ps(&y[1*8+j],_mm256_permute2f128_ps(y1,y5,0x20));
			     _mm256_store_ps(&y[5*8+j],_mm256_permute2f128_ps(y1,y5,0x31));
			     y6 = _mm256_shuffle_ps(z5,z7,_MM_SHUFFLE(1,0,1,0));
			     _mm256_store_ps(&y[2*8+j],_mm256_permute2f128_ps(y2,y6,0x20));
			     _mm256_store_ps(&y[6*8+j],_mm256_permute2f128_ps(y2,y6,0x31));
			     y7 = _mm256_shuffle_ps(z6,z7,_MM_SHUFFLE(3,2,3,2));
			     _mm256_store_ps(&y[3*8+j],_mm256_permute2f128_ps(y3,y7,0x20));
			     _mm256_store_ps(&y[7*8+j],_mm256_permute2f128_ps(y3,y7,0x31));
                             j += stride;
			 }
		  }




		 

