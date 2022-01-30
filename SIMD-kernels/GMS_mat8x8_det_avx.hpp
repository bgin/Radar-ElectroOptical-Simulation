

#ifndef __GMS_MAT8X8_DET_AVX_HPP__
#define __GMS_MAT8X8_DET_AVX_HPP__

#include <immintrin.h>
#include "GMS_config.h"

namespace gms {


         namespace  math {


#define NUM_MATS 8
#define NUM_MATS_HALF 4
#define MAT_WIDTH 4
#define MAT_HEIGHT 4
#define MAT_SIZE (MAT_WIDTH * MAT_HEIGHT)
	              /*
                           Based on Intel example.
                       */
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      mat8x8_det_a_ymm8r4(const float * __restrict input,
		                          float * __restrict output) {

	              __m256	a11, a12, a13, a14;
	              __m256	a21, a22, a23, a24;
	              __m256	a31, a32, a33, a34;
	              __m256	a41, a42, a43, a44;
	              __m256    ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	              __m256	ymma, ymmb, ymmc, ymmd, ymme, ymmf, ymmg, ymmh;

	


		      //Read the first eight rows of the 8x8 block of floats into ymm registers
		      	//Begin the transpose
		      ymm0 = _mm256_load_ps((float *) (input + (0*MAT_SIZE)));	//8 (first 2 rows) from matrix 1 
		      ymm1 = _mm256_load_ps((float *) (input + (1*MAT_SIZE)));	//8 (first 2 rows) from matrix 2
		      ymma = _mm256_unpacklo_ps(ymm0, ymm1);	
		      ymmb = _mm256_unpackhi_ps(ymm0, ymm1);	
		      ymm2 = _mm256_load_ps((float *) (input + (2*MAT_SIZE)));	//8 (first 2 rows) from matrix 3 
		      ymm3 = _mm256_load_ps((float *) (input + (3*MAT_SIZE)));	//8 (first 2 rows) from matrix 4
		      ymmc = _mm256_unpacklo_ps(ymm2, ymm3);	
		      ymmd = _mm256_unpackhi_ps(ymm2, ymm3);	
		      ymm4 = _mm256_load_ps((float *) (input + (4*MAT_SIZE)));	//8 (first 2 rows) from matrix 5 
		      ymm5 = _mm256_load_ps((float *) (input + (5*MAT_SIZE)));	//8 (first 2 rows) from matrix 6
		      ymme = _mm256_unpacklo_ps(ymm4, ymm5);	
		      ymmf = _mm256_unpackhi_ps(ymm4, ymm5);	
		      ymm6 = _mm256_load_ps((float *) (input + (6*MAT_SIZE)));	//8 (first 2 rows) from matrix 7 
		      ymm7 = _mm256_load_ps((float *) (input + (7*MAT_SIZE)));	//8 (first 2 rows) from matrix 8 
                      ymmg = _mm256_unpacklo_ps(ymm6, ymm7);	
		      ymmh = _mm256_unpackhi_ps(ymm6, ymm7);	
			
		//Create output rows 0, 1, 4 and 5
		      ymm7 = _mm256_shuffle_ps(ymma, ymmc, 0x4e);	
		      ymm5 = _mm256_blend_ps(ymm7, ymma, 0x33);	
		      ymm3 = _mm256_blend_ps(ymm7, ymmc, 0xcc);	
		
		      ymm6 = _mm256_shuffle_ps(ymme, ymmg, 0x4e);	
		      ymm4 = _mm256_blend_ps(ymm6, ymme, 0x33);	
		      ymm2 = _mm256_blend_ps(ymm6, ymmg, 0xcc);	

		//Last step to create rows 0, 1, 4, and 5
		     a11 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	// 11, 11, 11, 11, 11, 11, 11, 11
		     a12 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	// 12, 12, 12, 12, 12, 12, 12, 12
		     a21 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	// 13, 13, 13, 13, 13, 13, 13, 13
		     a22 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);	// 14, 14, 14, 14, 14, 14, 14, 14

		//Create output rows 2, 3, 6 and 7
		     ymm7 = _mm256_shuffle_ps(ymmb, ymmd, 0x4e);	
		     ymm5 = _mm256_blend_ps(ymm7, ymmb, 0x33);	
		     ymm3 = _mm256_blend_ps(ymm7, ymmd, 0xcc);

		     ymm6 = _mm256_shuffle_ps(ymmf, ymmh, 0x4e);	
		     ymm4 = _mm256_blend_ps(ymm6, ymmf, 0x33);	
		     ymm2 = _mm256_blend_ps(ymm6, ymmh, 0xcc);	

		//Last step to create rows 2, 3, 6, and 7
		     a13 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	
		     a14 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	
		     a23 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	
		     a24 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);

		//Read the first eight rows of the 8x8 block of floats into ymm registers
		     ymm0 = _mm256_load_ps((float *) (input + (0*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 1 
		     ymm1 = _mm256_load_ps((float *) (input + (1*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 2
		     ymma = _mm256_unpacklo_ps(ymm0, ymm1);	
		     ymmb = _mm256_unpackhi_ps(ymm0, ymm1);	
		     ymm2 = _mm256_load_ps((float *) (input + (2*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 3 
		     ymm3 = _mm256_load_ps((float *) (input + (3*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 4
		     ymmc = _mm256_unpacklo_ps(ymm2, ymm3);	
		     ymmd = _mm256_unpackhi_ps(ymm2, ymm3);	
		     ymm4 = _mm256_load_ps((float *) (input + (4*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 5 
		     ymm5 = _mm256_load_ps((float *) (input + (5*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 6
		     ymme = _mm256_unpacklo_ps(ymm4, ymm5);	
		     ymmf = _mm256_unpackhi_ps(ymm4, ymm5);	
		     ymm6 = _mm256_load_ps((float *) (input + (6*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 7 
		     ymm7 = _mm256_load_ps((float *) (input + (7*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 8 
                     ymmg = _mm256_unpacklo_ps(ymm6, ymm7);	
		     ymmh = _mm256_unpackhi_ps(ymm6, ymm7);	
	
				//Create output rows 0, 1, 4 and 5
		     ymm7 = _mm256_shuffle_ps(ymma, ymmc, 0x4e);	 
		     ymm5 = _mm256_blend_ps(ymm7, ymma, 0x33);	
		     ymm3 = _mm256_blend_ps(ymm7, ymmc, 0xcc);	 
		
		     ymm6 = _mm256_shuffle_ps(ymme, ymmg, 0x4e);	 
		     ymm4 = _mm256_blend_ps(ymm6, ymme, 0x33);	 
		     ymm2 = _mm256_blend_ps(ymm6, ymmg, 0xcc);	 

		//Last step to create rows 0, 1, 4, and 5
		     a31 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	
		     a32 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	
		     a41 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	
		     a42 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);	

		//Create output rows 2, 3, 6 and 7
		     ymm7 = _mm256_shuffle_ps(ymmb, ymmd, 0x4e);	
		     ymm5 = _mm256_blend_ps(ymm7, ymmb, 0x33);	
		     ymm3 = _mm256_blend_ps(ymm7, ymmd, 0xcc);	

		     ymm6 = _mm256_shuffle_ps(ymmf, ymmh, 0x4e);	
		     ymm4 = _mm256_blend_ps(ymm6, ymmf, 0x33);	
		     ymm2 = _mm256_blend_ps(ymm6, ymmh, 0xcc);	

		//Last step to create rows 2, 3, 6, and 7
		     a33 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	//33, 33, 33, 33, 33, 33, 33, 33
		     a34 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	
		     a43 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	
		     a44 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);	

		/// ******Determinants for a11*****************
		     __m256      a11a22, a33a44, a43a34, a11a23, a32a44, a34a42, a11a24, a32a43, a42a33;
		     __m256	a11a22_D, a11a23_D, a11a24_D, a11_D;
		     __m256	a33a44_a43a34, a32a44_a34a42, a32a43_a42a33;

		// calculate a11a22(a33a44 - a43a34) for 8 matrices
		     a11a22 = _mm256_mul_ps(a11, a22); 
		     a33a44 = _mm256_mul_ps(a33, a44);
		     a43a34 = _mm256_mul_ps(a43, a34);
		     a33a44_a43a34 = _mm256_sub_ps(a33a44, a43a34);
		     a11a22_D = _mm256_mul_ps(a11a22, a33a44_a43a34); 
		
		// calculate a11a23(a32a44 - a34a42) for 8 matrices
		     a11a23 = _mm256_mul_ps(a11, a23);
		     a32a44 = _mm256_mul_ps(a32, a44);
		     a34a42 = _mm256_mul_ps(a34, a42);
		     a32a44_a34a42 = _mm256_sub_ps(a32a44, a34a42);
		     a11a23_D = _mm256_mul_ps(a11a23, a32a44_a34a42); 
		
		// calculate a11a24(a32a43 - a42a33) for 8 matrices
		     a11a24 = _mm256_mul_ps(a11, a24);
		     a32a43 = _mm256_mul_ps(a32, a43);
		     a42a33 = _mm256_mul_ps(a42, a33);
		     a32a43_a42a33 = _mm256_sub_ps(a32a43, a42a33);
		     a11a24_D = _mm256_mul_ps(a11a24, a32a43_a42a33); 

		// calculate partial determinant for 8 matrices: (a11a22(a33a44 - a43a34) + 
		//a11a23(a32a44 - a34a42) + a11a24(a32a43 - a42a33)) 
		     a11_D = _mm256_add_ps(_mm256_sub_ps(a11a22_D, a11a23_D), a11a24_D);

		/// ******Determinants for a12*****************//
		     __m256	a12a21, a12a23, a31a44, a41a34, a12a24, a31a43, a41a33;
		     __m256	a12a21_D, a12a23_D, a12a24_D, a12_D;
		     __m256	a31a44_a41a34, a31a43_a41a33;

		// calculate a12a21(a33a44 - a43a34) for 8 matrices
		     a12a21 = _mm256_mul_ps(a12, a21);
		     a12a21_D = _mm256_mul_ps(a12a21, a33a44_a43a34); 
		
		// calculate a12a23(a31a44 - a41a34) for 8 matrices
		     a12a23 = _mm256_mul_ps(a12, a23);
		     a31a44 = _mm256_mul_ps(a31, a44);
		     a41a34 = _mm256_mul_ps(a41, a34);
		     a31a44_a41a34 = _mm256_sub_ps(a31a44, a41a34);
		     a12a23_D = _mm256_mul_ps(a12a23, a31a44_a41a34); 
		
		// calculate a12a24(a31a43 - a41a33) for 8 matrices
		     a12a24 = _mm256_mul_ps(a12, a24);
		     a31a43 = _mm256_mul_ps(a31, a43);
		     a41a33 = _mm256_mul_ps(a41, a33);
		     a31a43_a41a33 = _mm256_sub_ps(a31a43, a41a33);
		     a12a24_D = _mm256_mul_ps(a12a24, a31a43_a41a33); 

		// calculate partial determinant for 8 matrices: (a12a21(a33a44 - a43a34) + 
		// a12a23(a31a44 - a41a34) + a12a24(a31a43 - a41a33)) 
		     a12_D = _mm256_sub_ps(_mm256_sub_ps(a12a23_D, a12a21_D), a12a24_D); 

		/// ******Determinants for a13*****************//
		     __m256	a13a21, a13a22, a13a24, a31a42, a41a32;
		     __m256	a13a21_D, a13a22_D, a13a24_D, a13_D;
		     __m256	a31a42_a41a32;

		// calculate a13a21(a32a44 - a42a34) for 8 matrices
		     a13a21 = _mm256_mul_ps(a13, a21);
		     a13a21_D = _mm256_mul_ps(a13a21, a32a44_a34a42); 
		
		// calculate a13a22(a31a43 - a41a33) for 8 matrices
		     a13a22 = _mm256_mul_ps(a13, a22);
		     a13a22_D = _mm256_mul_ps(a13a22, a31a43_a41a33); 
		
		// calculate a13a24(a31a42 - a41a32) for 8 matrices
		     a13a24 = _mm256_mul_ps(a13, a24);
		     a31a42 = _mm256_mul_ps(a31, a42);
		     a41a32 = _mm256_mul_ps(a41, a32);
		     a31a42_a41a32 = _mm256_sub_ps(a31a42, a41a32);
		     a13a24_D = _mm256_mul_ps(a13a24, a31a42_a41a32); 

		// calculate partial determinant for 8 matrices: (a13a21(a32a44 - a42a34) + 
		// a13a22(a31a43 - a41a33) + a13a24(a31a42 - a41a32)) 
		     a13_D = _mm256_add_ps(_mm256_sub_ps(a13a21_D, a13a22_D), a13a24_D); 

		/// ******Determinants for a14*****************//
		     __m256	a14a21, a14a22, a14a23;
		     __m256	a14a21_D, a14a22_D, a14a23_D, a14_D;

		// calculate a14a21(a32a43 - a42a33) for 8 matrices
		     a14a21 = _mm256_mul_ps(a14, a21);
		     a14a21_D = _mm256_mul_ps(a14a21, a32a43_a42a33); 
		
		// calculate a14a22(a31a43 - a41a33) for 8 matrices
		     a14a22 = _mm256_mul_ps(a14, a22);
		     a14a22_D = _mm256_mul_ps(a14a22, a31a43_a41a33); 
		
		// calculate a14a23(a31a42 - a41a32) for 8 matrices
		     a14a23 = _mm256_mul_ps(a14, a23);
		     a14a23_D = _mm256_mul_ps(a14a23, a31a42_a41a32); 
		
		// calculate partial determinant for 8 matrices: (a14a21(a32a43 - a42a33)) + 
		// a14a22(a31a43 - a41a33) + a14a23(a31a42 - a41a32)) 
		     a14_D = _mm256_sub_ps(_mm256_sub_ps(a14a22_D, a14a21_D), a14a23_D); 

		// Calculate final Determinant value for 8 matrices: addition of 
		// determinants for a11, a12, a13, a14
		    a11_D = _mm256_add_ps(a11_D, a12_D);
		    a13_D = _mm256_add_ps(a13_D, a14_D);
		    a11_D = _mm256_add_ps(a11_D, a13_D); 

		    _mm256_store_ps((float *)output, a11_D);
	
	   }


	   // Unaligned version
              __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      mat8x8_det_u_ymm8r4(const float * __restrict input,
		                          float * __restrict output) {

	              __m256	a11, a12, a13, a14;
	              __m256	a21, a22, a23, a24;
	              __m256	a31, a32, a33, a34;
	              __m256	a41, a42, a43, a44;
	              __m256    ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	              __m256	ymma, ymmb, ymmc, ymmd, ymme, ymmf, ymmg, ymmh;

	


		      //Read the first eight rows of the 8x8 block of floats into ymm registers
		      	//Begin the transpose
		      ymm0 = _mm256_loadu_ps((float *) (input + (0*MAT_SIZE)));	//8 (first 2 rows) from matrix 1 
		      ymm1 = _mm256_loadu_ps((float *) (input + (1*MAT_SIZE)));	//8 (first 2 rows) from matrix 2
		      ymma = _mm256_unpacklo_ps(ymm0, ymm1);	
		      ymmb = _mm256_unpackhi_ps(ymm0, ymm1);	
		      ymm2 = _mm256_loadu_ps((float *) (input + (2*MAT_SIZE)));	//8 (first 2 rows) from matrix 3 
		      ymm3 = _mm256_loadu_ps((float *) (input + (3*MAT_SIZE)));	//8 (first 2 rows) from matrix 4
		      ymmc = _mm256_unpacklo_ps(ymm2, ymm3);	
		      ymmd = _mm256_unpackhi_ps(ymm2, ymm3);	
		      ymm4 = _mm256_loadu_ps((float *) (input + (4*MAT_SIZE)));	//8 (first 2 rows) from matrix 5 
		      ymm5 = _mm256_loadu_ps((float *) (input + (5*MAT_SIZE)));	//8 (first 2 rows) from matrix 6
		      ymme = _mm256_unpacklo_ps(ymm4, ymm5);	
		      ymmf = _mm256_unpackhi_ps(ymm4, ymm5);	
		      ymm6 = _mm256_loadu_ps((float *) (input + (6*MAT_SIZE)));	//8 (first 2 rows) from matrix 7 
		      ymm7 = _mm256_loadu_ps((float *) (input + (7*MAT_SIZE)));	//8 (first 2 rows) from matrix 8 
                      ymmg = _mm256_unpacklo_ps(ymm6, ymm7);	
		      ymmh = _mm256_unpackhi_ps(ymm6, ymm7);	
			
		//Create output rows 0, 1, 4 and 5
		      ymm7 = _mm256_shuffle_ps(ymma, ymmc, 0x4e);	
		      ymm5 = _mm256_blend_ps(ymm7, ymma, 0x33);	
		      ymm3 = _mm256_blend_ps(ymm7, ymmc, 0xcc);	
		
		      ymm6 = _mm256_shuffle_ps(ymme, ymmg, 0x4e);	
		      ymm4 = _mm256_blend_ps(ymm6, ymme, 0x33);	
		      ymm2 = _mm256_blend_ps(ymm6, ymmg, 0xcc);	

		//Last step to create rows 0, 1, 4, and 5
		     a11 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	// 11, 11, 11, 11, 11, 11, 11, 11
		     a12 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	// 12, 12, 12, 12, 12, 12, 12, 12
		     a21 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	// 13, 13, 13, 13, 13, 13, 13, 13
		     a22 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);	// 14, 14, 14, 14, 14, 14, 14, 14

		//Create output rows 2, 3, 6 and 7
		     ymm7 = _mm256_shuffle_ps(ymmb, ymmd, 0x4e);	
		     ymm5 = _mm256_blend_ps(ymm7, ymmb, 0x33);	
		     ymm3 = _mm256_blend_ps(ymm7, ymmd, 0xcc);

		     ymm6 = _mm256_shuffle_ps(ymmf, ymmh, 0x4e);	
		     ymm4 = _mm256_blend_ps(ymm6, ymmf, 0x33);	
		     ymm2 = _mm256_blend_ps(ymm6, ymmh, 0xcc);	

		//Last step to create rows 2, 3, 6, and 7
		     a13 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	
		     a14 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	
		     a23 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	
		     a24 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);

		//Read the first eight rows of the 8x8 block of floats into ymm registers
		     ymm0 = _mm256_loadu_ps((float *) (input + (0*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 1 
		     ymm1 = _mm256_loadu_ps((float *) (input + (1*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 2
		     ymma = _mm256_unpacklo_ps(ymm0, ymm1);	
		     ymmb = _mm256_unpackhi_ps(ymm0, ymm1);	
		     ymm2 = _mm256_loadu_ps((float *) (input + (2*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 3 
		     ymm3 = _mm256_loadu_ps((float *) (input + (3*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 4
		     ymmc = _mm256_unpacklo_ps(ymm2, ymm3);	
		     ymmd = _mm256_unpackhi_ps(ymm2, ymm3);	
		     ymm4 = _mm256_loadu_ps((float *) (input + (4*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 5 
		     ymm5 = _mm256_loadu_ps((float *) (input + (5*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 6
		     ymme = _mm256_unpacklo_ps(ymm4, ymm5);	
		     ymmf = _mm256_unpackhi_ps(ymm4, ymm5);	
		     ymm6 = _mm256_loadu_ps((float *) (input + (6*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 7 
		     ymm7 = _mm256_loadu_ps((float *) (input + (7*MAT_SIZE) + (2*MAT_WIDTH)));	//8 (next 2 rows) from matrix 8 
                     ymmg = _mm256_unpacklo_ps(ymm6, ymm7);	
		     ymmh = _mm256_unpackhi_ps(ymm6, ymm7);	
	
				//Create output rows 0, 1, 4 and 5
		     ymm7 = _mm256_shuffle_ps(ymma, ymmc, 0x4e);	 
		     ymm5 = _mm256_blend_ps(ymm7, ymma, 0x33);	
		     ymm3 = _mm256_blend_ps(ymm7, ymmc, 0xcc);	 
		
		     ymm6 = _mm256_shuffle_ps(ymme, ymmg, 0x4e);	 
		     ymm4 = _mm256_blend_ps(ymm6, ymme, 0x33);	 
		     ymm2 = _mm256_blend_ps(ymm6, ymmg, 0xcc);	 

		//Last step to create rows 0, 1, 4, and 5
		     a31 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	
		     a32 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	
		     a41 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	
		     a42 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);	

		//Create output rows 2, 3, 6 and 7
		     ymm7 = _mm256_shuffle_ps(ymmb, ymmd, 0x4e);	
		     ymm5 = _mm256_blend_ps(ymm7, ymmb, 0x33);	
		     ymm3 = _mm256_blend_ps(ymm7, ymmd, 0xcc);	

		     ymm6 = _mm256_shuffle_ps(ymmf, ymmh, 0x4e);	
		     ymm4 = _mm256_blend_ps(ymm6, ymmf, 0x33);	
		     ymm2 = _mm256_blend_ps(ymm6, ymmh, 0xcc);	

		//Last step to create rows 2, 3, 6, and 7
		     a33 = _mm256_permute2f128_ps(ymm5, ymm4, 0x20);	//33, 33, 33, 33, 33, 33, 33, 33
		     a34 = _mm256_permute2f128_ps(ymm3, ymm2, 0x20);	
		     a43 = _mm256_permute2f128_ps(ymm5, ymm4, 0x31);	
		     a44 = _mm256_permute2f128_ps(ymm3, ymm2, 0x31);	

		/// ******Determinants for a11*****************
		     __m256      a11a22, a33a44, a43a34, a11a23, a32a44, a34a42, a11a24, a32a43, a42a33;
		     __m256	a11a22_D, a11a23_D, a11a24_D, a11_D;
		     __m256	a33a44_a43a34, a32a44_a34a42, a32a43_a42a33;

		// calculate a11a22(a33a44 - a43a34) for 8 matrices
		     a11a22 = _mm256_mul_ps(a11, a22); 
		     a33a44 = _mm256_mul_ps(a33, a44);
		     a43a34 = _mm256_mul_ps(a43, a34);
		     a33a44_a43a34 = _mm256_sub_ps(a33a44, a43a34);
		     a11a22_D = _mm256_mul_ps(a11a22, a33a44_a43a34); 
		
		// calculate a11a23(a32a44 - a34a42) for 8 matrices
		     a11a23 = _mm256_mul_ps(a11, a23);
		     a32a44 = _mm256_mul_ps(a32, a44);
		     a34a42 = _mm256_mul_ps(a34, a42);
		     a32a44_a34a42 = _mm256_sub_ps(a32a44, a34a42);
		     a11a23_D = _mm256_mul_ps(a11a23, a32a44_a34a42); 
		
		// calculate a11a24(a32a43 - a42a33) for 8 matrices
		     a11a24 = _mm256_mul_ps(a11, a24);
		     a32a43 = _mm256_mul_ps(a32, a43);
		     a42a33 = _mm256_mul_ps(a42, a33);
		     a32a43_a42a33 = _mm256_sub_ps(a32a43, a42a33);
		     a11a24_D = _mm256_mul_ps(a11a24, a32a43_a42a33); 

		// calculate partial determinant for 8 matrices: (a11a22(a33a44 - a43a34) + 
		//a11a23(a32a44 - a34a42) + a11a24(a32a43 - a42a33)) 
		     a11_D = _mm256_add_ps(_mm256_sub_ps(a11a22_D, a11a23_D), a11a24_D);

		/// ******Determinants for a12*****************//
		     __m256	a12a21, a12a23, a31a44, a41a34, a12a24, a31a43, a41a33;
		     __m256	a12a21_D, a12a23_D, a12a24_D, a12_D;
		     __m256	a31a44_a41a34, a31a43_a41a33;

		// calculate a12a21(a33a44 - a43a34) for 8 matrices
		     a12a21 = _mm256_mul_ps(a12, a21);
		     a12a21_D = _mm256_mul_ps(a12a21, a33a44_a43a34); 
		
		// calculate a12a23(a31a44 - a41a34) for 8 matrices
		     a12a23 = _mm256_mul_ps(a12, a23);
		     a31a44 = _mm256_mul_ps(a31, a44);
		     a41a34 = _mm256_mul_ps(a41, a34);
		     a31a44_a41a34 = _mm256_sub_ps(a31a44, a41a34);
		     a12a23_D = _mm256_mul_ps(a12a23, a31a44_a41a34); 
		
		// calculate a12a24(a31a43 - a41a33) for 8 matrices
		     a12a24 = _mm256_mul_ps(a12, a24);
		     a31a43 = _mm256_mul_ps(a31, a43);
		     a41a33 = _mm256_mul_ps(a41, a33);
		     a31a43_a41a33 = _mm256_sub_ps(a31a43, a41a33);
		     a12a24_D = _mm256_mul_ps(a12a24, a31a43_a41a33); 

		// calculate partial determinant for 8 matrices: (a12a21(a33a44 - a43a34) + 
		// a12a23(a31a44 - a41a34) + a12a24(a31a43 - a41a33)) 
		     a12_D = _mm256_sub_ps(_mm256_sub_ps(a12a23_D, a12a21_D), a12a24_D); 

		/// ******Determinants for a13*****************//
		     __m256	a13a21, a13a22, a13a24, a31a42, a41a32;
		     __m256	a13a21_D, a13a22_D, a13a24_D, a13_D;
		     __m256	a31a42_a41a32;

		// calculate a13a21(a32a44 - a42a34) for 8 matrices
		     a13a21 = _mm256_mul_ps(a13, a21);
		     a13a21_D = _mm256_mul_ps(a13a21, a32a44_a34a42); 
		
		// calculate a13a22(a31a43 - a41a33) for 8 matrices
		     a13a22 = _mm256_mul_ps(a13, a22);
		     a13a22_D = _mm256_mul_ps(a13a22, a31a43_a41a33); 
		
		// calculate a13a24(a31a42 - a41a32) for 8 matrices
		     a13a24 = _mm256_mul_ps(a13, a24);
		     a31a42 = _mm256_mul_ps(a31, a42);
		     a41a32 = _mm256_mul_ps(a41, a32);
		     a31a42_a41a32 = _mm256_sub_ps(a31a42, a41a32);
		     a13a24_D = _mm256_mul_ps(a13a24, a31a42_a41a32); 

		// calculate partial determinant for 8 matrices: (a13a21(a32a44 - a42a34) + 
		// a13a22(a31a43 - a41a33) + a13a24(a31a42 - a41a32)) 
		     a13_D = _mm256_add_ps(_mm256_sub_ps(a13a21_D, a13a22_D), a13a24_D); 

		/// ******Determinants for a14*****************//
		     __m256	a14a21, a14a22, a14a23;
		     __m256	a14a21_D, a14a22_D, a14a23_D, a14_D;

		// calculate a14a21(a32a43 - a42a33) for 8 matrices
		     a14a21 = _mm256_mul_ps(a14, a21);
		     a14a21_D = _mm256_mul_ps(a14a21, a32a43_a42a33); 
		
		// calculate a14a22(a31a43 - a41a33) for 8 matrices
		     a14a22 = _mm256_mul_ps(a14, a22);
		     a14a22_D = _mm256_mul_ps(a14a22, a31a43_a41a33); 
		
		// calculate a14a23(a31a42 - a41a32) for 8 matrices
		     a14a23 = _mm256_mul_ps(a14, a23);
		     a14a23_D = _mm256_mul_ps(a14a23, a31a42_a41a32); 
		
		// calculate partial determinant for 8 matrices: (a14a21(a32a43 - a42a33)) + 
		// a14a22(a31a43 - a41a33) + a14a23(a31a42 - a41a32)) 
		     a14_D = _mm256_sub_ps(_mm256_sub_ps(a14a22_D, a14a21_D), a14a23_D); 

		// Calculate final Determinant value for 8 matrices: addition of 
		// determinants for a11, a12, a13, a14
		    a11_D = _mm256_add_ps(a11_D, a12_D);
		    a13_D = _mm256_add_ps(a13_D, a14_D);
		    a11_D = _mm256_add_ps(a11_D, a13_D); 

		    _mm256_storeu_ps((float *)output, a11_D);
	
	   }

	   
	            
     }

}

















#endif /*__GMS_MAT8X8_DET_AVX_HPP__*/
