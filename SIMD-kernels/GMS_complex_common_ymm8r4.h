
#ifndef __GMS_COMPLEX_COMMON_YMM8R4_H__
#define __GMS_COMPLEX_COMMON_YMM8R4_H__ 061020191431

namespace file_info {


	const unsigned int GMS_COMPLEX_COMMON_YMM8R4_MAJOR = 1U;

	const unsigned int GMS_COMPLEX_COMMON_YMM8R4_MINOR = 0U;

	const unsigned int GMS_COMPLEX_COMMON_YMM8R4_MICRO = 0U;

	const unsigned int GMS_COMPLEX_COMMON_YMM8R4_FULLVER = 
	 1000U*gGMS_COMPLEX_COMMON_YMM8R4_MAJOR+100U*gGMS_COMPLEX_COMMON_YMM8R4_MINOR+10U*gGMS_COMPLEX_COMMON_YMM8R4_MICRO;

	const char * const GMS_COMPLEX_COMMON_YMM8R4_CREATE_DATE = "06-10-2019 14:31 +00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const GMS_COMPLEX_COMMON_YMM8R4_BUILD_DATE = __DATE__ ":" __TIME__;

	const char * const GMS_COMPLEX_COMMON_YMM8R4_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const GMS_COMPLEX_COMMON_YMM8R4_SYNOPSIS = "Common procedures for AVX complex field (1D) classes.";

}

#include <cstdint>
#include <type_traits>
#include <immintrin.h>
#include "GMS_config.h"


namespace gms {
	namespace math {
			
#if !defined (AVX_COMPLEX_ADDITION)
#define AVX_COMPLEX_ADDITION(out,v1,v2,idx,off) \
	(out) = _mm256_add_ps(_mm256_mul_ps(_mm256_load_ps(&(v1).m_Re[(idx)+(off)]), \
	_mm256_load_ps(&(v2).m_Re[(idx)+(off)])), _mm256_mul_ps(_mm256_load_ps(&(v1).m_Im[(idx)+(off)]), \
	_mm256_load_ps(&(v2).m_Im[(idx)+(off)])));
#endif

#if !defined (AVX_COMPLEX_SUBTRACTION)
#define AVX_COMPLEX_SUBTRACTION(out,v1,v2,idx,off) \
	(out) = _mm256_sub_ps(_mm256_mul_ps(_mm256_load_ps(&(v1).m_Im[(idx)+(off)]), \
	_mm256_load_ps(&(v2).m_Re[(idx)+(off)])), _mm256_mul_ps(_mm512_load_ps(&(v1).m_Re[(idx)+(off)]), \
	_mm512_load_ps(&(v2).m_Im[(idx)+(off)])));
#endif

		// Warning macro parameter v2 must be an exact copy
		// of parameter v1. This should done by calling class (AVX512VComplex1D)
		// Move Constructor.
#if !defined (AVX_COMPLEX_MAGNITUDE)
#define AVX_COMPLEX_MAGNITUDE(out,v1,v2,idx,off) \
	(out) = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_load_ps(&(v1).m_Re[(idx)+(off)]), \
	_mm256_load_ps(&(v2).m_Re[(idx)+(off)])), _mm256_mul_ps(_mm256_load_ps(&(v1).m_Im[(idx)+(off)]), \
	_mm256_load_ps(&(v2).m_Im[(idx)+(off)]))));
#endif

#if !defined (AVX_HORIZONTAL_ADD)
#define AVX_HORIZONTAL_ADD(sum,ymmx,tmp) \
	(tmp) = _mm256_add_ps((ymmx), _mm256_permute_ps((ymmx),0x05)); \
	(sum) = _mm_cvtsd_f64(_mm_add_sd(_mm256_castpd256_ps128((tmp)),_mm256_extractf128_ps((tmp),1)));
#endif

			
			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cnormalize_prod(AVXVField1D &out,
					       const AVXVField1D &v1,
					       const AVXVField1D &v2,
					       const bool do_nt_store) {
				if (v1.m_nsize != v2.m_nsize) {return;}
				int32_t i;
				if (do_nt_store) {
					for (i = 0; i != ROUND_TO_EIGHT(v1.m_nsize, 8); i += 8) {
						const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
						const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
						const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
						const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
						const __m256 re_part( _mm256_sub_ps(
									 _mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2,ymm3)));
						const __m256 im_part(_mm256_add_ps(
									 _mm256_mul_ps(ymm2, ymm0), _mm256_mul_ps(ymm0,ymm3)));
						const __m256 sqrt_term1(_mm256_mul_ps(re_part,re_part));
						const __m256 sqrt_term2(_mm256_mul_ps(im_part,im_part));
						const __m256 mag_term(_mm256_sqrt_ps(_mm256_add_ps(sqrt_term1,sqrt_term2)));
						_mm256_stream_ps(&out.m_Re[i], _mm256_div_ps(re_part,mag_part));
						_mm256_stream_ps(&out.m_Im[i], _mm256_div_ps(im_part,mag_term));

					}
					_mm_sfence();
					// Warning remainder is cached upon store.
					for (; i != out.m_nsize; ++i) {
						const float re = (v1.m_Re[i] * v2.m_Re[i]) - (v1.m_Im[i] * v2.m_Im[i]);
						const float im = (v1.m_Im[i] * v2.m_Re[i]) + (v1.m_Re[i] * v2.m_Im[i]);
						const float mag = std::sqrt(re*re + im*im);
						out.m_Re[i] = re / mag;
						out.m_Im[i] = im / mag;
					}
				}
				else {
					for (i = 0; i != ROUND_TO_EIGHT(v1.m_nsize, 8); i += 8) {
						const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
						const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
						const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
						const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
						const __m256 re_part(_mm256_sub_ps(
							_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2, ymm3)));
						const __m256 im_part(_mm256_add_ps(
							_mm256_mul_ps(ymm2, ymm0), _mm256_mul_ps(ymm0, ymm3)));
						const __m256 sqrt_term1(_mm256_mul_ps(re_part, re_part));
						const __m256 sqrt_term2(_mm256_mul_ps(im_part, im_part));
						const __m256 mag_term(_mm256_sqrt_ps(_mm256_add_ps(sqrt_term1, sqrt_term2)));
						_mm256_store_ps(&out.m_Re[i], _mm256_div_ps(re_part, mag_part));
						_mm256_store_ps(&out.m_Im[i], _mm256_div_ps(im_part, mag_term));

					}
					for (; i != out.m_nsize; ++i) {
						const float re = (v1.m_Re[i] * v2.m_Re[i]) - (v1.m_Im[i] * v2.m_Im[i]);
						const float im = (v1.m_Im[i] * v2.m_Re[i]) + (v1.m_Re[i] * v2.m_Im[i]);
						const float mag = std::sqrt(re*re + im*im);
						out.m_Re[i] = re / mag;
						out.m_Im[i] = im / mag;
					}
				}
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmean_prod(std::complex<float> &mean,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2) {
				if (v1.m_nsize != v2.m_nsize) { return;}


			    __ATTR_ALIGN__(64) struct {
                                       	float sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					int32_t i; char pad[44];
			    }ca;

				for (ca.i = 0; i != ROUND_TO_EIGHT(v1.m_nsize, 8); ca.i += 8) {
					const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
					const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
					const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
					const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
					const __m256 re_part(_mm256_sub_ps(_mm256_mul_ps(ymm0,ymm1),
										  _mm256_mul_ps(ymm2,ymm3)));
					const __m256 im_part(_mm256_add_ps(_mm256_mul_ps(ymm2,ymm3),
										 _mm256_mul_ps(ymm0,ymm3)));
					const __m256 t1(_mm256_add_ps(re_part, _mm256_permute_ps(re_part,0x05)));
					ca.sumre = (float)(_mm_cvtsd_f64(_mm_add_ss(_mm256_castps256_ps128(t1), _mm256_extractf128_ps(t1,1))));
					ca.accre += ca.sumre;
					const __m256 t2(_mm256_add_ps(im_part, _mm256_permute_ps(im_part,0x05)));
					ca.sumim = (float)(_mm_cvtsd_f64(_mm_add_ss(_mm256_castps256_ps128(t2), _mm256_extractf128_ps(t2,1))));
					ca.accim += ca.sumim;
				}
				for (; ca.i != v1.m_nsize; ++i) {
					ca.accre += (v1.m_Re[ca.i] * v2.m_Re[ca.i]) - (v1.m_Im[ca.i] * v2.m_Im[ca.i]);
					ca.accim += (v1.m_Im[ca.i] * v2.m_Re[ca.i]) + (v1.m_Re[ca.i] * v2.m_Im[ca.i]);
				}
				mean._Val[0] = ca.accre /= static_cast<float>(v1.m_nsize);
				mean._Val[1] = ca.accim /= static_cast<float>(v1.m_nsize);
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmean_quot(std::complex<float> &mean,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2) {
				if (v1.m_nsize != v2.m_nsize) { return;}

				__ATTR_ALIGN__(64) struct {
                                        float sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{ 0.0 },
					t1{ 0.0 }, t2{0.0};
					int32_t i;
					char pad[12];
				}ca;

				for (ca.i = 0; i != ROUND_TO_EIGHT(v1.m_nsize, 8); ca.i += 8) {
					const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
					const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
					const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
					const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
					const __m256 re_part(_mm256_sub_ps(_mm256_mul_ps(ymm0,ymm1),
										              _mm256_mul_ps(ymm2,ymm3)));
					const __m256 im_part(_mm256_add_ps(_mm256_mul_ps(ymm2,ymm1),
													  _mm256_mul_ps(ymm0,ymm3)));
					const __m256 den_term(_mm256_add_ps(_mm256_mul_ps(ymm1,ymm1),
													   _mm256_mul_ps(ymm3,ymm3)));
					const __m256 re_quot(_mm256_div_ps(re_part,den_term));
					const __m256 tmp1(_mm256_add_ps(re_quot, _mm256_permute_ps(re_quot,0x05)));
					ca.sumre = _mm_cvtsd_f64(_mm_add_ss(_mm256_castps256_ps128(tmp1), _mm256_extractf128_ps(tmp1,1)));
					ca.accre += ca.sumre;
					const __m256 im_quot(_mm256_div_ps(im_part,den_term));
					const __m256 tmp2(_mm256_add_ps(im_quot, _mm256_permute_ps(im_quot,0x05)));
					ca.sumim = _mm_cvtsd_f64(_mm_add_sd(_mm256_castps256_ps128(tmp2), _mm256_extractf128_ps(tmp2,1)));
					ca.accre += ca.sumim;
				}
				for (; ca.i != v1.m_nsize; ++ca.i) {
					const float re = (v1.m_Re[ca.i] * v2.m_Re[ca.i]) - (v1.m_Im[ca.i] * v2.m_Im[ca.i]);
					const float im = (v1.m_Im[ca.i] * v2.m_Re[ca.i]) + (v1.m_Re[ca.i] * v2.m_Im[ca.i]);
					const float den = (v2.m_Re[ca.i] * v2.m_Re[ca.i]) + (v2.m_Im[ca.i] * v2.m_Im[ca.i]);
					ca.t1 += re / den;
					ca.t2 += im / den;
				}
				ca.accre += ca.t1;
				ca.accim += ca.t2;
				mean._Val[0] = ca.accre / static_cast<float>(v1.m_nsize);
				mean._Val[1] = ca.accim / static_cast<float>(v1.m_nsize);
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cconj_prod(AVXVField1D &out,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2,
					  const bool do_nt_store) {
				if (v1.m_nsize != v2.m_nsize) { return;}
				int32_t i;
				if (do_nt_store) {
					for (i = 0; i != ROUND_TO_EIGHT(out.m_nsize, 8); i += 8) {
						const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
						const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
						const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
						const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
						_mm256_stream_ps(&out.m_Re[i], _mm256_add_ps(
							_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2,ymm3)));
						_mm256_stream_ps(&out.m_Im[i], _mm256_sub_ps(
							_mm256_mul_ps(ymm2, ymm1), _mm256_mul_ps(ymm0,ymm3)));
					}
					_mm_sfence();
					for (; i != v1.m_nsize; ++i) { // Cache storing remainder.
						out.m_Re[i] = (v1.m_Re[i] * v2.m_Re[i]) + (v1.m_Im[i] * v2.m_Im[i]);
						out.m_Im[i] = (v1.m_Im[i] * v2.m_Re[i]) - (v1.m_Re[i] * v2.m_Im[i]);
					}
				}
				else {
					for (i = 0; i != ROUND_TO_EIGHT(out.m_nsize, 8); i += 8) {
						const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
						const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
						const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
						const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
						_mm256_store_ps(&out.m_Re[i], _mm256_add_ps(
							_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2, ymm3)));
						_mm256_store_ps(&out.m_Im[i], _mm256_sub_ps(
							_mm256_mul_ps(ymm2, ymm1), _mm256_mul_ps(ymm0, ymm3)));
					}
					for (; i != v1.m_nsize; ++i) { 
						out.m_Re[i] = (v1.m_Re[i] * v2.m_Re[i]) + (v1.m_Im[i] * v2.m_Im[i]);
						out.m_Im[i] = (v1.m_Im[i] * v2.m_Re[i]) - (v1.m_Re[i] * v2.m_Im[i]);
					}
				}
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cnorm_conjprod(AVXVField1D &out,
					      const AVXVField1D &v1,
					      const AVXVField1D &v2,
					      const bool do_nt_store) {
				if (v1.m_nsize != v2.m_nsize) { return;}
				int32_t i;
				if (do_nt_store) {
					for (i = 0; i != ROUND_TO_EIGHT(out.m_nsize, 8); i += 8) {
						const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
						const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
						const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
						const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
						const __m256 re_part(_mm256_add_ps(_mm256_mul_ps(ymm0,ymm1),
									  _mm256_mul_ps(ymm2,ymm3)));
						const __m256 im_part(_mm256_sub_ps(_mm256_mul_ps(ymm2,ymm1),
									 _mm256_mul_ps(ymm0,ymm3)));
						const __m256 mag_c1(_mm256_mul_ps(re_part,re_part));
						const __m256 mag_c2(_mm256_mul_ps(im_part,im_part));
						const __m256 vcmag(_mm256_sqrt_ps(_mm256_add_ps(mag_c1,mag_c2)));
						_mm256_stream_ps(&out.m_Re[i], _mm256_div_ps(re_part,vcmag));
						_mm256_stream_ps(&out.m_Im[i], _mm256_div_ps(im_part,vcmag));
					}
					_mm_sfence();
					for (; i != out.m_nsize; ++i) { // Cache storing remainder.
						const float re = (v1.m_Re[i] * v2.m_Re[i]) - (v1.m_Im[i] * v2.m_Im[i]);
						const float im = (v1.m_Im[i] * v2.m_Re[i]) + (v1.m_Re[i] * v2.m_Im[i]);
						const float mag = std::sqrt(re*re + im*im);
						out.m_Re[i] = re / mag;
						out.m_Im[i] = im / mag;
					}
				}
				else {
					for (i = 0; i != ROUND_TO_EIGHT(out.m_nsize, 8); i += 8) {
						const __m256 ymm0(_mm256_load_ps(&v1.m_Re[i]));
						const __m256 ymm1(_mm256_load_ps(&v2.m_Re[i]));
						const __m256 ymm2(_mm256_load_ps(&v1.m_Im[i]));
						const __m256 ymm3(_mm256_load_ps(&v2.m_Im[i]));
						const __m256 re_part(_mm256_add_ps(_mm256_mul_ps(ymm0, ymm1),
							_mm256_mul_ps(ymm2, ymm3)));
						const __m256 im_part(_mm256_sub_ps(_mm256_mul_ps(ymm2, ymm1),
							_mm256_mul_ps(ymm0, ymm3)));
						const __m256 mag_c1(_mm256_mul_ps(re_part, re_part));
						const __m256 mag_c2(_mm256_mul_ps(im_part, im_part));
						const __m256 vcmag(_mm256_sqrt_ps(_mm256_add_ps(mag_c1, mag_c2)));
						_mm256_store_ps(&out.m_Re[i], _mm256_div_ps(re_part, vcmag));
						_mm256_store_ps(&out.m_Im[i], _mm256_div_ps(im_part, vcmag));
					}
					for (; i != out.m_nsize; ++i) { 
						const float re = (v1.m_Re[i] * v2.m_Re[i]) - (v1.m_Im[i] * v2.m_Im[i]);
						const float im = (v1.m_Im[i] * v2.m_Re[i]) + (v1.m_Re[i] * v2.m_Im[i]);
						const float mag = std::sqrt(re*re + im*im);
						out.m_Re[i] = re / mag;
						out.m_Im[i] = im / mag;
					}
				}
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmean_conjprod(std::complex<float> &mean,
					      const AVXVField1D &v1,
					      const AVXVField1D &v2) {
				if (v1.m_nsize != v2.m_nsize) { return;}

				__ATTR_ALIGN__(64) struct {
                                        float sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					int32_t i;
					char pad[44];
				}ca;

				__m256 re = _mm256_setzero_ps(), tmp1 = _mm256_setzero_ps();
				__m256 im = _mm256_setzero_ps(), tmp2 = _mm256_setzero_ps();
				for (ca.i = 0; i != ROUND_TO_EIGHT(v1.m_nsize, 8); i += 8) {
					AVX_COMPLEX_ADDITION(re,v1,v2,ca.i,0)
					AVX_COMPLEX_SUBTRACTION(im,v1,v2,ca.i,0)
					AVX_HORIZONTAL_ADD(ca.sumre,re,tmp1)
					ca.accre += ca.sum_re;
					AVX_HORIZONTAL_ADD(ca.sumim,im,tmp2)
					ca.accim += ca.sumim;
				}	
				for (; i != v1.m_nsize; ++i) {
					const float re = (v1.m_Re[i] * v2.m_Re[i]) + (v1.m_Im[i] * v2.m_Im[i]);
					const float im = (v1.m_Im[i] * v2.m_Re[i]) - (v1.m_Re[i] * v2.m_Im[i]);
					ca.accre += re;
					ca.accim += im;
				}
				mean._Val[0] = ca.accre / static_cast<float>(v1.m_nsize);
				mean._Val[1] = ca.accim / static_cast<float>(v1.m_nsize);
			}	
				
			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_arith_mean(std::complex<float> &mean,
					  const AVXVField1D &v) {

				__ATTR_ALIGN__(64) struct {
                                        float sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					int32_t i;
					char pad[44];
				}ca;
				__attribute__((align(64))) struct {
                                        __m256 re_part,im_part,tmp1,tmp2;
				}ca2;

				ca2.tmp1 = _mm256_setzero_ps();
				ca2.tmp2 = _mm256_setzero_ps();
				for (ca.i = 0; i != ROUND_TO_EIGHT(v.m_nsize, 8); ca.i += 8) {
					ca2.re_part = _mm256_load_ps(&v.m_Re[i]);
					AVX_HORIZONTAL_ADD(ca.sumre,ca2.re_part,ca2.tmp1)
					ca.accre += ca.sumre;
					ca2.im_part = _mm256_load_ps(&v.m_Im[i]);
					AVX_HORIZONTAL_ADD(ca.sumim,ca2.im_part,tmp2)
					ca.accim += ca.sumim;
				}
				for (; ca.i != v.m_nsize; ++ca.i) {
					ca.accre += v.m_Re[i];
					ca.accim += v.m_Im[i];
				}
				mean._Val[0] = ca.accre / static_cast<float>(v.m_nsize);
				mean._Val[1] = ca.accim / static_cast<float>(v.m_nsize);
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cnormalize(AVXVField1D &out,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2,
					  const bool do_nt_store) {
				if (v1.m_nsize != v2.m_nsize) { return; }
				__m256 cvmag(_mm256_setzero_ps());
				int32_t i;
				if (do_nt_store) {

					for (i = 0, i != ROUND_TO_EIGHT(out.m_nsize,8); i += 8) {
						AVX_COMPLEX_MAGNITUDE(cvmag, v1, v2, i, 0)
							const __m256 re_part(_mm256_load_ps(&v1.m_Re[i]));
						_mm256_stream_ps(&out.m_Re[i], _mm256_div_ps(re_part, cvmag));
						const __m256 im_part(_mm256_load_ps(&v1.m_Im[i]));
						_mm256_stream_ps(&out.m_Im[i], _mm256_div_ps(im_part, cvmag));
					}
					_mm_sfence();
					for (; i != v1.m_nsize; ++i) {
						const float cmag = std::sqrt(v1.m_Re[i] * v2.m_Re[i] + v1.m_Im[i] * v2.m_Im[i]);
						out.m_Re[i] = v1.m_Re[i] / cmag;
						out.m_Im[i] = v1.m_Im[i] / cmag;
					}
				}
				else {
					for (i = 0, i != ROUND_TO_EIGHT(out.m_nsize, 8); i += 8) {
						AVX_COMPLEX_MAGNITUDE(cvmag, v, cv, i, 0)
							const __m256 re_part(_mm256_load_ps(&v1.m_Re[i]));
						_mm256_store_ps(&out.m_Re[i], _mm256_div_ps(re_part, cvmag));
						    const __m256 im_part(_mm256_load_ps(&v1.m_Im[i]));
						_mm256_store_ps(&out.m_Im[i], _mm256_div_ps(im_part, cvmag));
					}
					for (; i != v1.m_nsize; ++i) {
						const float cmag = std::sqrt(v1.m_Re[i] * v2.m_Re[i] + v1.m_Im[i] * v2.m_Im[i]);
						out.m_Re[i] = v1.m_Re[i] / cmag;
						out.m_Im[i] = v1.m_Im[i] / cmag;
					}
				}
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmagnitude(AVXVField1D &out,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2) {
				if (v1.m_nsize != v2.m_nsize) { return;}
				__m256 vcmag(_mm256_setzero_ps());
				int32_t i;
				for (i = 0; i != ROUND_TO_EIGHT(out.m_nsize, 8); i += 8) {
					    AVX_COMPLEX_MAGNITUDE(vcmag,v1,v2,i,0)
						_mm256_store_ps(&out.m_Re[i], vcmag);
				}
				for (; i != out.m_nsize; ++i) {
					out.m_Re[i] = std::sqrt(v1.m_Re[i] * v2.m_Re[i] + v1.m_Im[i] * v2.m_Im[i]);
				}
			}
			
	}
}





#endif /*__GMS_COMPLEX_COMMON_YMM8R4_H__*/
