
#ifndef __GMS_AVX512_TRANSPOSITION_16X16_HPP__
#define __GMS_AVX512_TRANSPOSITION_16X16_HPP__

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
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
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
						      __m512 &x15) {

                           register __m512 y00,y01,y02,y03;
			   register __m512 y04,y05,y06,y07;
			   register __m512 y08,y09,y10,y11;
			   register __m512 y12,y13,y14,y15;



		           y00 = _mm512_unpacklo_ps(x00,x01);
			   y01 = _mm512_unpackhi_ps(x00,x01);
			   y02 = _mm512_unpacklo_ps(x02,x03);
			   y03 = _mm512_unpackhi_ps(x02,x03);
			   y04 = _mm512_unpacklo_ps(x04,x05);
			   y05 = _mm512_unpackhi_ps(x04,x05);
			   y06 = _mm512_unpacklo_ps(x06,x07);
			   y07 = _mm512_unpackhi_ps(0x6,x07);
			   y08 = _mm512_unpacklo_ps(x08,x09);
			   y09 = _mm512_unpackhi_ps(x08,x09);
			   y10 = _mm512_unpacklo_ps(x10,x11);
			   y11 = _mm512_unpackhi_ps(x10,x11);
			   y12 = _mm512_unpacklo_ps(x12,x13);
			   y13 = _mm512_unpackhi_ps(x12,x13);
			   y14 = _mm512_unpacklo_ps(x14,x15);
			   y15 = _mm512_unpackhi_ps(x14,x15);

			   x00 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(1,0,1,0));
			   x01 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(3,2,3,2));
			   x02 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(1,0,1,0));
			   x03 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(3,2,3,2));
			   x04 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(1,0,1,0));
			   x05 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(3,2,3,2));
			   x06 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(1,0,1,0));
			   x07 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(3,2,3,2));
			   x08 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(1,0,1,0));
			   x09 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(3,2,3,2));
			   x10 = _mm512_shuffle_ps(y09,y11,_MM_SHUFFLE(1,0,1,0));
			   x11 = _mm512_shuffle_ps(y09,y11,_MM_SHUFFLE(3,2,3,2));
			   x12 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(1,0,1,0));
			   x13 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(3,2,3,2));
			   x14 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(1,0,1,0));
			   x15 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(3,2,3,2));

			   y00 = _mm512_shuffle_f32x4(x00,x04,0x88);
			   y01 = _mm512_shuffle_f32x4(x01,x05,0x88);
			   y02 = _mm512_shuffle_f32x4(x02,x06,0x88);
                           y03 = _mm512_shuffle_f32x4(x03,x07,0x88);
			   y04 = _mm512_shuffle_f32x4(x00,x04,0xdd);
			   y05 = _mm512_shuffle_f32x4(x01,x05,0xdd);
			   y06 = _mm512_shuffle_f32x4(x02,x06,0xdd);
			   y07 = _mm512_shuffle_f32x4(x03,x07,0xdd);
			   y08 = _mm512_shuffle_f32x4(x08,x12,0x88);
			   y09 = _mm512_shuffle_f32x4(x09,x13,0x88);
			   y10 = _mm512_shuffle_f32x4(x10,x14,0x88);
			   y11 = _mm512_shuffle_f32x4(x11,x15,0x88);
			   y12 = _mm512_shuffle_f32x4(x08,x12,0xdd);
			   y13 = _mm512_shuffle_f32x4(x09,x13,0xdd);
			   y14 = _mm512_shuffle_f32x4(x10,x14,0xdd);
			   y15 = _mm512_shuffle_f32x4(x11,x15,0xdd);

			   x00 = _mm512_shuffle_f32x4(y00,y88,0x88);
			   x01 = _mm512_shuffle_f32x4(y01,y09,0x88);
			   x02 = _mm512_shuffle_f32x4(y02,y10,0x88);
			   x03 = _mm512_shuffle_f32x4(y03,y11,0x88);
			   x04 = _mm512_shuffle_f32x4(y04,y12,0x88);
			   x05 = _mm512_shuffle_f32x4(y05,y13,0x88);
			   x06 = _mm512_shuffle_f32x4(y06,y14,0x88);
			   x07 = _mm512_shuffle_f32x4(y07,y15,0x88);
			   x08 = _mm512_shuffle_f32x4(y00,y08,0xdd);
			   x09 = _mm512_shuffle_f32x4(y01,y09,0xdd);
			   x10 = _mm512_shuffle_f32x4(y02,y10,0xdd);
			   x11 = _mm512_shuffle_f32x4(y03,y11,0xdd);
			   x12 = _mm512_shuffle_f32x4(y04,y12,0xdd);
			   x13 = _mm512_shuffle_f32x4(y05,y13,0xdd);
			   x14 = _mm512_shuffle_f32x4(y06,y14,0xdd);
			   x15 = _mm512_shuffle_f32x4(y07,y15,0xdd);

			   
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_u_zmm16r4_16x16_ip(float * __restrict x00,
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
							float * __restrict x15) {

                           register __m512 y00,y01,y02,y03;
			   register __m512 y04,y05,y06,y07;
			   register __m512 y08,y09,y10,y11;
			   register __m512 y12,y13,y14,y15;
			   register __m512 z00,z01,z02,z03;
			   register __m512 z04,z05,z06,z07;
			   register __m512 z08,z09,z10,z11;
			   register __m512 z12,z13,z14,z15;

			   z00 = _mm512_loadu_ps(x00);
			   z01 = _mm512_loadu_ps(x01);
			   y00 = _mm512_unpacklo_ps(z00,z01);
			   y01 = _mm512_unpackhi_ps(z00,z01);
			   z02 = _mm512_loadu_ps(x02);
			   z03 = _mm512_loadu_ps(x03);
			   y02 = _mm512_unpacklo_ps(z02,z03);
			   y03 = _mm512_unpackhi_ps(z02,z03);
			   z04 = _mm512_loadu_ps(x04);
			   z05 = _mm512_loadu_ps(x05);
			   y04 = _mm512_unpacklo_ps(z04,z05);
			   y05 = _mm512_unpackhi_ps(z04,z05);
			   z06 = _mm512_loadu_ps(x06);
			   z07 = _mm512_loadu_ps(x07);
			   y06 = _mm512_unpacklo_ps(z06,z07);
			   y07 = _mm512_unpackhi_ps(z06,z07);
			   z08 = _mm512_loadu_ps(x08);
			   z09 = _mm512_loadu_ps(x09);
			   y08 = _mm512_unpacklo_ps(z08,z09);
			   y09 = _mm512_unpackhi_ps(z08,z09);
			   z10 = _mm512_loadu_ps(x10);
			   z11 = _mm512_loadu_ps(x11);
			   y10 = _mm512_unpacklo_ps(z10,z11);
			   y11 = _mm512_unpackhi_ps(z10,z11);
			   z12 = _mm512_loadu_ps(x12);
			   z13 = _mm512_loadu_ps(x13);
			   y12 = _mm512_unpacklo_ps(z12,z13);
			   y13 = _mm512_unpackhi_ps(z12,z13);
			   z14 = _mm512_loadu_ps(x14);
			   z15 = _mm512_loadu_ps(x15);
			   y14 = _mm512_unpacklo_ps(z14,z15);
			   y15 = _mm512_unpackhi_ps(z14,z15);

			   z00 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(1,0,1,0));
			   z01 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(3,2,3,2));
			   z02 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(1,0,1,0));
			   z03 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(3,2,3,2));
			   z04 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(1,0,1,0));
			   z05 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(3,2,3,2));
			   z06 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(1,0,1,0));
			   z07 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(3,2,3,2));
			   z08 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(1,0,1,0));
			   z09 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(3,2,3,2));
			   z10 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(1,0,1,0));
			   z11 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(3,2,3,2));
			   z12 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(1,0,1,0));
			   z13 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(3,2,3,2));
			   z14 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(1,0,1,0));
			   z15 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(3,2,3,2));

			   y00 = _mm512_shuffle_f32x4(z00,z04,0x88);
			   y01 = _mm512_shuffle_f32x4(z01,z05,0x88);
			   y02 = _mm512_shuffle_f32x4(z02,z06,x088);
			   y03 = _mm512_shuffle_f32x4(z03,z07,0x88);
			   y04 = _mm512_shuffle_f32x4(z00,z04,0xdd);
			   y05 = _mm512_shuffle_f32x4(z01,z05,0xdd);
			   y06 = _mm512_shuffle_f32x4(z02,z06,0xdd);
			   y07 = _mm512_shuffle_f32x4(z03,z07,x0dd);
			   y08 = _mm512_shuffle_f32x4(z08,z12,0x88);
			   _mm512_storeu_ps(x00,_mm512_shuffle_f32x4(y00,y08,0x88);
			   _mm512_storeu_ps(x08,_mm512_shuffle_f32x4(y00,y08,0xdd);
			   y09 = _mm512_shuffle_f32x4(z09,z13,0x88);
			   _mm512_storeu_ps(x09,_mm512_shuffle_f32x4(y01,y09,0xdd);
			   _mm512_storeu_ps(x01,_mm512_shuffle_f32x4(y01,y09,0x88);
			   y10 = _mm512_shuffle_f32x4(z10,z14,0x88);
			   _mm512_storeu_ps(x02,_mm512_shuffle_f32x4(y02,y10,0x88);
			   _mm512_storeu_ps(x10,_mm512_shuffle_f32x4(y02,y10,0xdd);
			   y11 = _mm512_shuffle_f32x4(z11,z15,0x88);
			   _mm512_storeu_ps(x03,_mm512_shuffle_f32x4(y03,y11,0x88);
			   _mm512_storeu_ps(x11,_mm512_shuffle_f32x4(y03,y11,0xdd);
			   y12 = _mm512_shuffle_f32x4(z08,z12,0xdd);
			   _mm512_storeu_ps(x04,_mm512_shuffle_f32x4(y04,y12,0x88);
			   _mm512_storeu_ps(x12,_mm512_shuffle_f32x4(y04,y12,0xdd);
			   y13 = _mm512_shuffle_f32x4(z09,z13,0xdd);
			   _mm512_storeu_ps(x05,_mm512_shuffle_f32x4(y05,y13,0x88);
			   _mm512_storeu_ps(x13,_mm512_shuffle_f32x4(y05,y13,0xdd);
			   y14 = _mm512_shuffle_f32x4(z10,z14,0xdd);
			   _mm512_storeu_ps(x06,_mm512_shuffle_f32x4(y06,y14,0x88);
			   _mm512_storeu_ps(x14,_mm512_shuffle_f32x4(y06,y14,0xdd);
			   y15 = _mm512_shuffle_f32x4(z11,z15,0xdd);
			   _mm512_storeu_ps(x07,_mm512_shuffle_f32x4(y07,y15,0x88);
			   _mm512_storeu_ps(x15,_mm512_shuffle_f32x4(y07,y15,0xdd);

			   
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
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
							float * __restrict __ATTR_ALIGN__(64) x15) {

                           register __m512 y00,y01,y02,y03;
			   register __m512 y04,y05,y06,y07;
			   register __m512 y08,y09,y10,y11;
			   register __m512 y12,y13,y14,y15;
			   register __m512 z00,z01,z02,z03;
			   register __m512 z04,z05,z06,z07;
			   register __m512 z08,z09,z10,z11;
			   register __m512 z12,z13,z14,z15;

			   z00 = _mm512_load_ps(x00);
			   z01 = _mm512_load_ps(x01);
			   y00 = _mm512_unpacklo_ps(z00,z01);
			   y01 = _mm512_unpackhi_ps(z00,z01);
			   z02 = _mm512_load_ps(x02);
			   z03 = _mm512_load_ps(x03);
			   y02 = _mm512_unpacklo_ps(z02,z03);
			   y03 = _mm512_unpackhi_ps(z02,z03);
			   z04 = _mm512_load_ps(x04);
			   z05 = _mm512_load_ps(x05);
			   y04 = _mm512_unpacklo_ps(z04,z05);
			   y05 = _mm512_unpackhi_ps(z04,z05);
			   z06 = _mm512_load_ps(x06);
			   z07 = _mm512_load_ps(x07);
			   y06 = _mm512_unpacklo_ps(z06,z07);
			   y07 = _mm512_unpackhi_ps(z06,z07);
			   z08 = _mm512_load_ps(x08);
			   z09 = _mm512_load_ps(x09);
			   y08 = _mm512_unpacklo_ps(z08,z09);
			   y09 = _mm512_unpackhi_ps(z08,z09);
			   z10 = _mm512_load_ps(x10);
			   z11 = _mm512_load_ps(x11);
			   y10 = _mm512_unpacklo_ps(z10,z11);
			   y11 = _mm512_unpackhi_ps(z10,z11);
			   z12 = _mm512_load_ps(x12);
			   z13 = _mm512_load_ps(x13);
			   y12 = _mm512_unpacklo_ps(z12,z13);
			   y13 = _mm512_unpackhi_ps(z12,z13);
			   z14 = _mm512_load_ps(x14);
			   z15 = _mm512_load_ps(x15);
			   y14 = _mm512_unpacklo_ps(z14,z15);
			   y15 = _mm512_unpackhi_ps(z14,z15);

			   z00 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(1,0,1,0));
			   z01 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(3,2,3,2));
			   z02 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(1,0,1,0));
			   z03 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(3,2,3,2));
			   z04 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(1,0,1,0));
			   z05 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(3,2,3,2));
			   z06 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(1,0,1,0));
			   z07 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(3,2,3,2));
			   z08 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(1,0,1,0));
			   z09 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(3,2,3,2));
			   z10 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(1,0,1,0));
			   z11 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(3,2,3,2));
			   z12 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(1,0,1,0));
			   z13 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(3,2,3,2));
			   z14 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(1,0,1,0));
			   z15 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(3,2,3,2));

			   y00 = _mm512_shuffle_f32x4(z00,z04,0x88);
			   y01 = _mm512_shuffle_f32x4(z01,z05,0x88);
			   y02 = _mm512_shuffle_f32x4(z02,z06,x088);
			   y03 = _mm512_shuffle_f32x4(z03,z07,0x88);
			   y04 = _mm512_shuffle_f32x4(z00,z04,0xdd);
			   y05 = _mm512_shuffle_f32x4(z01,z05,0xdd);
			   y06 = _mm512_shuffle_f32x4(z02,z06,0xdd);
			   y07 = _mm512_shuffle_f32x4(z03,z07,x0dd);
			   y08 = _mm512_shuffle_f32x4(z08,z12,0x88);
			   _mm512_store_ps(x00,_mm512_shuffle_f32x4(y00,y08,0x88);
			   _mm512_store_ps(x08,_mm512_shuffle_f32x4(y00,y08,0xdd);
			   y09 = _mm512_shuffle_f32x4(z09,z13,0x88);
			   _mm512_store_ps(x09,_mm512_shuffle_f32x4(y01,y09,0xdd);
			   _mm512_store_ps(x01,_mm512_shuffle_f32x4(y01,y09,0x88);
			   y10 = _mm512_shuffle_f32x4(z10,z14,0x88);
			   _mm512_store_ps(x02,_mm512_shuffle_f32x4(y02,y10,0x88);
			   _mm512_store_ps(x10,_mm512_shuffle_f32x4(y02,y10,0xdd);
			   y11 = _mm512_shuffle_f32x4(z11,z15,0x88);
			   _mm512_store_ps(x03,_mm512_shuffle_f32x4(y03,y11,0x88);
			   _mm512_store_ps(x11,_mm512_shuffle_f32x4(y03,y11,0xdd);
			   y12 = _mm512_shuffle_f32x4(z08,z12,0xdd);
			   _mm512_store_ps(x04,_mm512_shuffle_f32x4(y04,y12,0x88);
			   _mm512_store_ps(x12,_mm512_shuffle_f32x4(y04,y12,0xdd);
			   y13 = _mm512_shuffle_f32x4(z09,z13,0xdd);
			   _mm512_store_ps(x05,_mm512_shuffle_f32x4(y05,y13,0x88);
			   _mm512_store_ps(x13,_mm512_shuffle_f32x4(y05,y13,0xdd);
			   y14 = _mm512_shuffle_f32x4(z10,z14,0xdd);
			   _mm512_store_ps(x06,_mm512_shuffle_f32x4(y06,y14,0x88);
			   _mm512_store_ps(x14,_mm512_shuffle_f32x4(y06,y14,0xdd);
			   y15 = _mm512_shuffle_f32x4(z11,z15,0xdd);
			   _mm512_store_ps(x07,_mm512_shuffle_f32x4(y07,y15,0x88);
			   _mm512_store_ps(x15,_mm512_shuffle_f32x4(y07,y15,0xdd);

			   
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_u_zmm16r4_16x16_ip(float * __restrict x) {
		                                       
                           register __m512 y00,y01,y02,y03;
			   register __m512 y04,y05,y06,y07;
			   register __m512 y08,y09,y10,y11;
			   register __m512 y12,y13,y14,y15;
			   register __m512 z00,z01,z02,z03;
			   register __m512 z04,z05,z06,z07;
			   register __m512 z08,z09,z10,z11;
			   register __m512 z12,z13,z14,z15;

			   z00 = _mm512_loadu_ps(&x[0*16]);
			   z01 = _mm512_loadu_ps(&x[1*16]);
			   y00 = _mm512_unpacklo_ps(z00,z01);
			   y01 = _mm512_unpackhi_ps(z00,z01);
			   z02 = _mm512_loadu_ps(&x[2*16]);
			   z03 = _mm512_loadu_ps(&x[3*16]);
			   y02 = _mm512_unpacklo_ps(z02,z03);
			   y03 = _mm512_unpackhi_ps(z02,z03);
			   z04 = _mm512_loadu_ps(&x[4*16]);
			   z05 = _mm512_loadu_ps(&x[5*16]);
			   y04 = _mm512_unpacklo_ps(z04,z05);
			   y05 = _mm512_unpackhi_ps(z04,z05);
			   z06 = _mm512_loadu_ps(&x[6*16]);
			   z07 = _mm512_loadu_ps(&x[7*16]);
			   y06 = _mm512_unpacklo_ps(z06,z07);
			   y07 = _mm512_unpackhi_ps(z06,z07);
			   z08 = _mm512_loadu_ps(&x[8*16]);
			   z09 = _mm512_loadu_ps(&x[9*16]);
			   y08 = _mm512_unpacklo_ps(z08,z09);
			   y09 = _mm512_unpackhi_ps(z08,z09);
			   z10 = _mm512_loadu_ps(&x[10*16]);
			   z11 = _mm512_loadu_ps(&x[11*16]);
			   y10 = _mm512_unpacklo_ps(z10,z11);
			   y11 = _mm512_unpackhi_ps(z10,z11);
			   z12 = _mm512_loadu_ps(&x[12*16]);
			   z13 = _mm512_loadu_ps(&x[13*16]);
			   y12 = _mm512_unpacklo_ps(z12,z13);
			   y13 = _mm512_unpackhi_ps(z12,z13);
			   z14 = _mm512_loadu_ps(&x[14*16]);
			   z15 = _mm512_loadu_ps(&x[15*16]);
			   y14 = _mm512_unpacklo_ps(z14,z15);
			   y15 = _mm512_unpackhi_ps(z14,z15);

			   z00 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(1,0,1,0));
			   z01 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(3,2,3,2));
			   z02 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(1,0,1,0));
			   z03 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(3,2,3,2));
			   z04 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(1,0,1,0));
			   z05 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(3,2,3,2));
			   z06 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(1,0,1,0));
			   z07 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(3,2,3,2));
			   z08 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(1,0,1,0));
			   z09 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(3,2,3,2));
			   z10 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(1,0,1,0));
			   z11 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(3,2,3,2));
			   z12 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(1,0,1,0));
			   z13 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(3,2,3,2));
			   z14 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(1,0,1,0));
			   z15 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(3,2,3,2));

			   y00 = _mm512_shuffle_f32x4(z00,z04,0x88);
			   y01 = _mm512_shuffle_f32x4(z01,z05,0x88);
			   y02 = _mm512_shuffle_f32x4(z02,z06,x088);
			   y03 = _mm512_shuffle_f32x4(z03,z07,0x88);
			   y04 = _mm512_shuffle_f32x4(z00,z04,0xdd);
			   y05 = _mm512_shuffle_f32x4(z01,z05,0xdd);
			   y06 = _mm512_shuffle_f32x4(z02,z06,0xdd);
			   y07 = _mm512_shuffle_f32x4(z03,z07,x0dd);
			   y08 = _mm512_shuffle_f32x4(z08,z12,0x88);
			   _mm512_storeu_ps(&x[0*16],_mm512_shuffle_f32x4(y00,y08,0x88);
			   _mm512_storeu_ps(&x[8*16],_mm512_shuffle_f32x4(y00,y08,0xdd);
			   y09 = _mm512_shuffle_f32x4(z09,z13,0x88);
			   _mm512_storeu_ps(&x[9*16],_mm512_shuffle_f32x4(y01,y09,0xdd);
			   _mm512_storeu_ps(&x[1*16],_mm512_shuffle_f32x4(y01,y09,0x88);
			   y10 = _mm512_shuffle_f32x4(z10,z14,0x88);
			   _mm512_storeu_ps(&x[2*16],_mm512_shuffle_f32x4(y02,y10,0x88);
			   _mm512_storeu_ps(&x[10*16],_mm512_shuffle_f32x4(y02,y10,0xdd);
			   y11 = _mm512_shuffle_f32x4(z11,z15,0x88);
			   _mm512_storeu_ps(&x[3*16],_mm512_shuffle_f32x4(y03,y11,0x88);
			   _mm512_storeu_ps(&x[11*16],_mm512_shuffle_f32x4(y03,y11,0xdd);
			   y12 = _mm512_shuffle_f32x4(z08,z12,0xdd);
			   _mm512_storeu_ps(&x[4*16],_mm512_shuffle_f32x4(y04,y12,0x88);
			   _mm512_storeu_ps(&x[12*16],_mm512_shuffle_f32x4(y04,y12,0xdd);
			   y13 = _mm512_shuffle_f32x4(z09,z13,0xdd);
			   _mm512_storeu_ps(&x[5*16],_mm512_shuffle_f32x4(y05,y13,0x88);
			   _mm512_storeu_ps(&x[13*16],_mm512_shuffle_f32x4(y05,y13,0xdd);
			   y14 = _mm512_shuffle_f32x4(z10,z14,0xdd);
			   _mm512_storeu_ps(&x[6*16],_mm512_shuffle_f32x4(y06,y14,0x88);
			   _mm512_storeu_ps(&x[14*16],_mm512_shuffle_f32x4(y06,y14,0xdd);
			   y15 = _mm512_shuffle_f32x4(z11,z15,0xdd);
			   _mm512_storeu_ps(&x[7*16],_mm512_shuffle_f32x4(y07,y15,0x88);
			   _mm512_storeu_ps(&x[15*16],_mm512_shuffle_f32x4(y07,y15,0xdd);

			   
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_a_zmm16r4_16x16_ip(float * __restrict __ATTR_ALIGN__(64) x) {
		                                       
                           register __m512 y00,y01,y02,y03;
			   register __m512 y04,y05,y06,y07;
			   register __m512 y08,y09,y10,y11;
			   register __m512 y12,y13,y14,y15;
			   register __m512 z00,z01,z02,z03;
			   register __m512 z04,z05,z06,z07;
			   register __m512 z08,z09,z10,z11;
			   register __m512 z12,z13,z14,z15;

			   z00 = _mm512_load_ps(&x[0*16]);
			   z01 = _mm512_load_ps(&x[1*16]);
			   y00 = _mm512_unpacklo_ps(z00,z01);
			   y01 = _mm512_unpackhi_ps(z00,z01);
			   z02 = _mm512_load_ps(&x[2*16]);
			   z03 = _mm512_load_ps(&x[3*16]);
			   y02 = _mm512_unpacklo_ps(z02,z03);
			   y03 = _mm512_unpackhi_ps(z02,z03);
			   z04 = _mm512_load_ps(&x[4*16]);
			   z05 = _mm512_load_ps(&x[5*16]);
			   y04 = _mm512_unpacklo_ps(z04,z05);
			   y05 = _mm512_unpackhi_ps(z04,z05);
			   z06 = _mm512_load_ps(&x[6*16]);
			   z07 = _mm512_load_ps(&x[7*16]);
			   y06 = _mm512_unpacklo_ps(z06,z07);
			   y07 = _mm512_unpackhi_ps(z06,z07);
			   z08 = _mm512_load_ps(&x[8*16]);
			   z09 = _mm512_load_ps(&x[9*16]);
			   y08 = _mm512_unpacklo_ps(z08,z09);
			   y09 = _mm512_unpackhi_ps(z08,z09);
			   z10 = _mm512_load_ps(&x[10*16]);
			   z11 = _mm512_load_ps(&x[11*16]);
			   y10 = _mm512_unpacklo_ps(z10,z11);
			   y11 = _mm512_unpackhi_ps(z10,z11);
			   z12 = _mm512_load_ps(&x[12*16]);
			   z13 = _mm512_load_ps(&x[13*16]);
			   y12 = _mm512_unpacklo_ps(z12,z13);
			   y13 = _mm512_unpackhi_ps(z12,z13);
			   z14 = _mm512_load_ps(&x[14*16]);
			   z15 = _mm512_load_ps(&x[15*16]);
			   y14 = _mm512_unpacklo_ps(z14,z15);
			   y15 = _mm512_unpackhi_ps(z14,z15);

			   z00 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(1,0,1,0));
			   z01 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(3,2,3,2));
			   z02 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(1,0,1,0));
			   z03 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(3,2,3,2));
			   z04 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(1,0,1,0));
			   z05 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(3,2,3,2));
			   z06 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(1,0,1,0));
			   z07 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(3,2,3,2));
			   z08 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(1,0,1,0));
			   z09 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(3,2,3,2));
			   z10 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(1,0,1,0));
			   z11 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(3,2,3,2));
			   z12 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(1,0,1,0));
			   z13 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(3,2,3,2));
			   z14 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(1,0,1,0));
			   z15 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(3,2,3,2));

			   y00 = _mm512_shuffle_f32x4(z00,z04,0x88);
			   y01 = _mm512_shuffle_f32x4(z01,z05,0x88);
			   y02 = _mm512_shuffle_f32x4(z02,z06,x088);
			   y03 = _mm512_shuffle_f32x4(z03,z07,0x88);
			   y04 = _mm512_shuffle_f32x4(z00,z04,0xdd);
			   y05 = _mm512_shuffle_f32x4(z01,z05,0xdd);
			   y06 = _mm512_shuffle_f32x4(z02,z06,0xdd);
			   y07 = _mm512_shuffle_f32x4(z03,z07,x0dd);
			   y08 = _mm512_shuffle_f32x4(z08,z12,0x88);
			   _mm512_store_ps(&x[0*16],_mm512_shuffle_f32x4(y00,y08,0x88);
			   _mm512_store_ps(&x[8*16],_mm512_shuffle_f32x4(y00,y08,0xdd);
			   y09 = _mm512_shuffle_f32x4(z09,z13,0x88);
			   _mm512_store_ps(&x[9*16],_mm512_shuffle_f32x4(y01,y09,0xdd);
			   _mm512_store_ps(&x[1*16],_mm512_shuffle_f32x4(y01,y09,0x88);
			   y10 = _mm512_shuffle_f32x4(z10,z14,0x88);
			   _mm512_store_ps(&x[2*16],_mm512_shuffle_f32x4(y02,y10,0x88);
			   _mm512_store_ps(&x[10*16],_mm512_shuffle_f32x4(y02,y10,0xdd);
			   y11 = _mm512_shuffle_f32x4(z11,z15,0x88);
			   _mm512_store_ps(&x[3*16],_mm512_shuffle_f32x4(y03,y11,0x88);
			   _mm512_store_ps(&x[11*16],_mm512_shuffle_f32x4(y03,y11,0xdd);
			   y12 = _mm512_shuffle_f32x4(z08,z12,0xdd);
			   _mm512_store_ps(&x[4*16],_mm512_shuffle_f32x4(y04,y12,0x88);
			   _mm512_store_ps(&x[12*16],_mm512_shuffle_f32x4(y04,y12,0xdd);
			   y13 = _mm512_shuffle_f32x4(z09,z13,0xdd);
			   _mm512_store_ps(&x[5*16],_mm512_shuffle_f32x4(y05,y13,0x88);
			   _mm512_store_ps(&x[13*16],_mm512_shuffle_f32x4(y05,y13,0xdd);
			   y14 = _mm512_shuffle_f32x4(z10,z14,0xdd);
			   _mm512_store_ps(&x[6*16],_mm512_shuffle_f32x4(y06,y14,0x88);
			   _mm512_store_ps(&x[14*16],_mm512_shuffle_f32x4(y06,y14,0xdd);
			   y15 = _mm512_shuffle_f32x4(z11,z15,0xdd);
			   _mm512_store_ps(&x[7*16],_mm512_shuffle_f32x4(y07,y15,0x88);
			   _mm512_store_ps(&x[15*16],_mm512_shuffle_f32x4(y07,y15,0xdd);

			   
		    }


		      __ATTR_REGCALL__
		      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_u_zmm16r4_16x16(float * __restrict x,
		                                     float * __restrict y){
		                                       
                           register __m512 y00,y01,y02,y03;
			   register __m512 y04,y05,y06,y07;
			   register __m512 y08,y09,y10,y11;
			   register __m512 y12,y13,y14,y15;
			   register __m512 z00,z01,z02,z03;
			   register __m512 z04,z05,z06,z07;
			   register __m512 z08,z09,z10,z11;
			   register __m512 z12,z13,z14,z15;

			   z00 = _mm512_loadu_ps(&x[0*16]);
			   z01 = _mm512_loadu_ps(&x[1*16]);
			   y00 = _mm512_unpacklo_ps(z00,z01);
			   y01 = _mm512_unpackhi_ps(z00,z01);
			   z02 = _mm512_loadu_ps(&x[2*16]);
			   z03 = _mm512_loadu_ps(&x[3*16]);
			   y02 = _mm512_unpacklo_ps(z02,z03);
			   y03 = _mm512_unpackhi_ps(z02,z03);
			   z04 = _mm512_loadu_ps(&x[4*16]);
			   z05 = _mm512_loadu_ps(&x[5*16]);
			   y04 = _mm512_unpacklo_ps(z04,z05);
			   y05 = _mm512_unpackhi_ps(z04,z05);
			   z06 = _mm512_loadu_ps(&x[6*16]);
			   z07 = _mm512_loadu_ps(&x[7*16]);
			   y06 = _mm512_unpacklo_ps(z06,z07);
			   y07 = _mm512_unpackhi_ps(z06,z07);
			   z08 = _mm512_loadu_ps(&x[8*16]);
			   z09 = _mm512_loadu_ps(&x[9*16]);
			   y08 = _mm512_unpacklo_ps(z08,z09);
			   y09 = _mm512_unpackhi_ps(z08,z09);
			   z10 = _mm512_loadu_ps(&x[10*16]);
			   z11 = _mm512_loadu_ps(&x[11*16]);
			   y10 = _mm512_unpacklo_ps(z10,z11);
			   y11 = _mm512_unpackhi_ps(z10,z11);
			   z12 = _mm512_loadu_ps(&x[12*16]);
			   z13 = _mm512_loadu_ps(&x[13*16]);
			   y12 = _mm512_unpacklo_ps(z12,z13);
			   y13 = _mm512_unpackhi_ps(z12,z13);
			   z14 = _mm512_loadu_ps(&x[14*16]);
			   z15 = _mm512_loadu_ps(&x[15*16]);
			   y14 = _mm512_unpacklo_ps(z14,z15);
			   y15 = _mm512_unpackhi_ps(z14,z15);

			   z00 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(1,0,1,0));
			   z01 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(3,2,3,2));
			   z02 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(1,0,1,0));
			   z03 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(3,2,3,2));
			   z04 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(1,0,1,0));
			   z05 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(3,2,3,2));
			   z06 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(1,0,1,0));
			   z07 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(3,2,3,2));
			   z08 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(1,0,1,0));
			   z09 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(3,2,3,2));
			   z10 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(1,0,1,0));
			   z11 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(3,2,3,2));
			   z12 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(1,0,1,0));
			   z13 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(3,2,3,2));
			   z14 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(1,0,1,0));
			   z15 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(3,2,3,2));

			   y00 = _mm512_shuffle_f32x4(z00,z04,0x88);
			   y01 = _mm512_shuffle_f32x4(z01,z05,0x88);
			   y02 = _mm512_shuffle_f32x4(z02,z06,x088);
			   y03 = _mm512_shuffle_f32x4(z03,z07,0x88);
			   y04 = _mm512_shuffle_f32x4(z00,z04,0xdd);
			   y05 = _mm512_shuffle_f32x4(z01,z05,0xdd);
			   y06 = _mm512_shuffle_f32x4(z02,z06,0xdd);
			   y07 = _mm512_shuffle_f32x4(z03,z07,x0dd);
			   y08 = _mm512_shuffle_f32x4(z08,z12,0x88);
			   _mm512_storeu_ps(&y[0*16],_mm512_shuffle_f32x4(y00,y08,0x88);
			   _mm512_storeu_ps(&y[8*16],_mm512_shuffle_f32x4(y00,y08,0xdd);
			   y09 = _mm512_shuffle_f32x4(z09,z13,0x88);
			   _mm512_storeu_ps(&y[9*16],_mm512_shuffle_f32x4(y01,y09,0xdd);
			   _mm512_storeu_ps(&y[1*16],_mm512_shuffle_f32x4(y01,y09,0x88);
			   y10 = _mm512_shuffle_f32x4(z10,z14,0x88);
			   _mm512_storeu_ps(&y[2*16],_mm512_shuffle_f32x4(y02,y10,0x88);
			   _mm512_storeu_ps(&y[10*16],_mm512_shuffle_f32x4(y02,y10,0xdd);
			   y11 = _mm512_shuffle_f32x4(z11,z15,0x88);
			   _mm512_storeu_ps(&y[3*16],_mm512_shuffle_f32x4(y03,y11,0x88);
			   _mm512_storeu_ps(&y[11*16],_mm512_shuffle_f32x4(y03,y11,0xdd);
			   y12 = _mm512_shuffle_f32x4(z08,z12,0xdd);
			   _mm512_storeu_ps(&y[4*16],_mm512_shuffle_f32x4(y04,y12,0x88);
			   _mm512_storeu_ps(&y[12*16],_mm512_shuffle_f32x4(y04,y12,0xdd);
			   y13 = _mm512_shuffle_f32x4(z09,z13,0xdd);
			   _mm512_storeu_ps(&y[5*16],_mm512_shuffle_f32x4(y05,y13,0x88);
			   _mm512_storeu_ps(&y[13*16],_mm512_shuffle_f32x4(y05,y13,0xdd);
			   y14 = _mm512_shuffle_f32x4(z10,z14,0xdd);
			   _mm512_storeu_ps(&y[6*16],_mm512_shuffle_f32x4(y06,y14,0x88);
			   _mm512_storeu_ps(&y[14*16],_mm512_shuffle_f32x4(y06,y14,0xdd);
			   y15 = _mm512_shuffle_f32x4(z11,z15,0xdd);
			   _mm512_storeu_ps(&y[7*16],_mm512_shuffle_f32x4(y07,y15,0x88);
			   _mm512_storeu_ps(&y[15*16],_mm512_shuffle_f32x4(y07,y15,0xdd);

			   
		    }


		      __ATTR_REGCALL__
		      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_u_zmm16r4_16x16(float * __restrict __ATTR_ALIGN__(64)  x,
		                                     float * __restrict __ATTR_ALIGN__(64)  y){
		                                       
                           register __m512 y00,y01,y02,y03;
			   register __m512 y04,y05,y06,y07;
			   register __m512 y08,y09,y10,y11;
			   register __m512 y12,y13,y14,y15;
			   register __m512 z00,z01,z02,z03;
			   register __m512 z04,z05,z06,z07;
			   register __m512 z08,z09,z10,z11;
			   register __m512 z12,z13,z14,z15;

			   z00 = _mm512_load_ps(&x[0*16]);
			   z01 = _mm512_load_ps(&x[1*16]);
			   y00 = _mm512_unpacklo_ps(z00,z01);
			   y01 = _mm512_unpackhi_ps(z00,z01);
			   z02 = _mm512_load_ps(&x[2*16]);
			   z03 = _mm512_load_ps(&x[3*16]);
			   y02 = _mm512_unpacklo_ps(z02,z03);
			   y03 = _mm512_unpackhi_ps(z02,z03);
			   z04 = _mm512_load_ps(&x[4*16]);
			   z05 = _mm512_load_ps(&x[5*16]);
			   y04 = _mm512_unpacklo_ps(z04,z05);
			   y05 = _mm512_unpackhi_ps(z04,z05);
			   z06 = _mm512_load_ps(&x[6*16]);
			   z07 = _mm512_load_ps(&x[7*16]);
			   y06 = _mm512_unpacklo_ps(z06,z07);
			   y07 = _mm512_unpackhi_ps(z06,z07);
			   z08 = _mm512_load_ps(&x[8*16]);
			   z09 = _mm512_load_ps(&x[9*16]);
			   y08 = _mm512_unpacklo_ps(z08,z09);
			   y09 = _mm512_unpackhi_ps(z08,z09);
			   z10 = _mm512_load_ps(&x[10*16]);
			   z11 = _mm512_load_ps(&x[11*16]);
			   y10 = _mm512_unpacklo_ps(z10,z11);
			   y11 = _mm512_unpackhi_ps(z10,z11);
			   z12 = _mm512_load_ps(&x[12*16]);
			   z13 = _mm512_load_ps(&x[13*16]);
			   y12 = _mm512_unpacklo_ps(z12,z13);
			   y13 = _mm512_unpackhi_ps(z12,z13);
			   z14 = _mm512_load_ps(&x[14*16]);
			   z15 = _mm512_load_ps(&x[15*16]);
			   y14 = _mm512_unpacklo_ps(z14,z15);
			   y15 = _mm512_unpackhi_ps(z14,z15);

			   z00 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(1,0,1,0));
			   z01 = _mm512_shuffle_ps(y00,y02,_MM_SHUFFLE(3,2,3,2));
			   z02 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(1,0,1,0));
			   z03 = _mm512_shuffle_ps(y01,y03,_MM_SHUFFLE(3,2,3,2));
			   z04 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(1,0,1,0));
			   z05 = _mm512_shuffle_ps(y04,y06,_MM_SHUFFLE(3,2,3,2));
			   z06 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(1,0,1,0));
			   z07 = _mm512_shuffle_ps(y05,y07,_MM_SHUFFLE(3,2,3,2));
			   z08 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(1,0,1,0));
			   z09 = _mm512_shuffle_ps(y08,y10,_MM_SHUFFLE(3,2,3,2));
			   z10 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(1,0,1,0));
			   z11 = _mm512_shuffle_ps(y09,y10,_MM_SHUFFLE(3,2,3,2));
			   z12 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(1,0,1,0));
			   z13 = _mm512_shuffle_ps(y12,y14,_MM_SHUFFLE(3,2,3,2));
			   z14 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(1,0,1,0));
			   z15 = _mm512_shuffle_ps(y13,y15,_MM_SHUFFLE(3,2,3,2));

			   y00 = _mm512_shuffle_f32x4(z00,z04,0x88);
			   y01 = _mm512_shuffle_f32x4(z01,z05,0x88);
			   y02 = _mm512_shuffle_f32x4(z02,z06,x088);
			   y03 = _mm512_shuffle_f32x4(z03,z07,0x88);
			   y04 = _mm512_shuffle_f32x4(z00,z04,0xdd);
			   y05 = _mm512_shuffle_f32x4(z01,z05,0xdd);
			   y06 = _mm512_shuffle_f32x4(z02,z06,0xdd);
			   y07 = _mm512_shuffle_f32x4(z03,z07,x0dd);
			   y08 = _mm512_shuffle_f32x4(z08,z12,0x88);
			   _mm512_store_ps(&y[0*16],_mm512_shuffle_f32x4(y00,y08,0x88);
			   _mm512_store_ps(&y[8*16],_mm512_shuffle_f32x4(y00,y08,0xdd);
			   y09 = _mm512_shuffle_f32x4(z09,z13,0x88);
			   _mm512_store_ps(&y[9*16],_mm512_shuffle_f32x4(y01,y09,0xdd);
			   _mm512_store_ps(&y[1*16],_mm512_shuffle_f32x4(y01,y09,0x88);
			   y10 = _mm512_shuffle_f32x4(z10,z14,0x88);
			   _mm512_store_ps(&y[2*16],_mm512_shuffle_f32x4(y02,y10,0x88);
			   _mm512_store_ps(&y[10*16],_mm512_shuffle_f32x4(y02,y10,0xdd);
			   y11 = _mm512_shuffle_f32x4(z11,z15,0x88);
			   _mm512_store_ps(&y[3*16],_mm512_shuffle_f32x4(y03,y11,0x88);
			   _mm512_store_ps(&y[11*16],_mm512_shuffle_f32x4(y03,y11,0xdd);
			   y12 = _mm512_shuffle_f32x4(z08,z12,0xdd);
			   _mm512_store_ps(&y[4*16],_mm512_shuffle_f32x4(y04,y12,0x88);
			   _mm512_store_ps(&y[12*16],_mm512_shuffle_f32x4(y04,y12,0xdd);
			   y13 = _mm512_shuffle_f32x4(z09,z13,0xdd);
			   _mm512_store_ps(&y[5*16],_mm512_shuffle_f32x4(y05,y13,0x88);
			   _mm512_store_ps(&y[13*16],_mm512_shuffle_f32x4(y05,y13,0xdd);
			   y14 = _mm512_shuffle_f32x4(z10,z14,0xdd);
			   _mm512_store_ps(&y[6*16],_mm512_shuffle_f32x4(y06,y14,0x88);
			   _mm512_store_ps(&y[14*16],_mm512_shuffle_f32x4(y06,y14,0xdd);
			   y15 = _mm512_shuffle_f32x4(z11,z15,0xdd);
			   _mm512_store_ps(&y[7*16],_mm512_shuffle_f32x4(y07,y15,0x88);
			   _mm512_store_ps(&y[15*16],_mm512_shuffle_f32x4(y07,y15,0xdd);

			   
		    }







		    
               
     }

}











#endif /*__GMS_AVX512_TRANSPOSITION_16X16_HPP__*/
