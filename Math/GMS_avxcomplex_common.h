
#ifndef __GMS_AVXCOMPLEX_COMMON_H__
#define __GMS_AVXCOMPLEX_COMMON_H__

namespace file_info {
#if defined _WIN64  
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

	const unsigned int gGMS_AVXCOMPLEX_COMMON_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_AVXCOMPLEX_COMMON_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_AVXCOMPLEX_COMMON_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_AVXCOMPLEX_COMMON_FULLVER = 
	 1000U*gGMS_AVXCOMPLEX_COMMON_MAJOR+100U*gGMS_AVXCOMPLEX_COMMON_MINOR+10U*gGMS_AVXCOMPLEX_COMMON_MICRO;

	const char * const pgGMS_AVXCOMPLEX_COMMON_CREATE_DATE = "06-10-2019 14:31 +00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_AVXCOMPLEX_COMMON_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_AVXCOMPLEX_COMMON_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_AVXCOMPLEX_COMMON_SYNOPSIS = "Common procedures for AVX complex field (1D) classes.";

}

#include <cstdint>
#include <complex>
#include <type_traits>
#include <immintrin.h>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#endif

namespace gms {
	namespace math {
			
#if !defined (AVX_COMPLEX_ADDITION)
#define AVX_COMPLEX_ADDITION(out,v1,v2,idx,off) \
	(out) = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&(v1).data.m_Re[(idx)+(off)]), \
	_mm256_load_pd(&(v2).data.m_Re[(idx)+(off)])), _mm256_mul_pd(_mm256_load_pd(&(v1).data.m_Im[(idx)+(off)]), \
	_mm256_load_pd(&(v2).data.m_Im[(idx)+(off)])));
#endif

#if !defined (AVX_COMPLEX_SUBTRACTION)
#define AVX_COMPLEX_SUBTRACTION(out,v1,v2,idx,off) \
	(out) = _mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&(v1).data.m_Im[(idx)+(off)]), \
	_mm256_load_pd(&(v2).data.m_Re[(idx)+(off)])), _mm256_mul_pd(_mm512_load_pd(&(v1).data.m_Re[(idx)+(off)]), \
	_mm512_load_pd(&(v2).data.m_Im[(idx)+(off)])));
#endif

		// Warning macro parameter v2 must be an exact copy
		// of parameter v1. This should done by calling class (AVX512VComplex1D)
		// Move Constructor.
#if !defined (AVX_COMPLEX_MAGNITUDE)
#define AVX_COMPLEX_MAGNITUDE(out,v1,v2,idx,off) \
	(out) = _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&(v1).data.m_Re[(idx)+(off)]), \
	_mm256_load_pd(&(v2).data.m_Re[(idx)+(off)])), _mm256_mul_pd(_mm256_load_pd(&(v1).data.m_Im[(idx)+(off)]), \
	_mm256_load_pd(&(v2).data.m_Im[(idx)+(off)]))));
#endif

#if !defined (AVX_HORIZONTAL_ADD)
#define AVX_HORIZONTAL_ADD(sum,ymmx,tmp) \
	(tmp) = _mm256_add_pd((ymmx), _mm256_permute_pd((ymmx),0x05)); \
	(sum) = _mm_cvtsd_f64(_mm_add_sd(_mm256_castpd256_pd128((tmp)),_mm256_extractf128_pd((tmp),1)));
#endif

			
			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cnormalize_prod(AVXVField1D &out,
					       const AVXVField1D &v1,
					       const AVXVField1D &v2,
					       const bool do_nt_store) {
				if (v1.data.m_nsize != v2.data.m_nsize) {return;}
				int32_t i;
				if (do_nt_store) {
					for (i = 0; i != ROUND_TO_FOUR(v1.data.m_nsize, 4); i += 4) {
						const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
						const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
						const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
						const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
						const __m256d re_part( _mm256_sub_pd(
									 _mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2,ymm3)));
						const __m256d im_part(_mm256_add_pd(
									 _mm256_mul_pd(ymm2, ymm0), _mm256_mul_pd(ymm0,ymm3)));
						const __m256d sqrt_term1(_mm256_mul_pd(re_part,re_part));
						const __m256d sqrt_term2(_mm256_mul_pd(im_part,im_part));
						const __m256d mag_term(_mm256_sqrt_pd(_mm256_add_pd(sqrt_term1,sqrt_term2)));
						_mm256_stream_pd(&out.data.m_Re[i], _mm256_div_pd(re_part,mag_part));
						_mm256_stream_pd(&out.data.m_Im[i], _mm256_div_pd(im_part,mag_term));

					}
					_mm_sfence();
					// Warning remainder is cached upon store.
					for (; i != out.data.m_nsize; ++i) {
						const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) - (v1.data.m_Im[i] * v2.data.m_Im[i]);
						const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) + (v1.data.m_Re[i] * v2.data.m_Im[i]);
						const double mag = std::sqrt(re*re + im*im);
						out.data.m_Re[i] = re / mag;
						out.data.m_Im[i] = im / mag;
					}
				}
				else {
					for (i = 0; i != ROUND_TO_FOUR(v1.data.m_nsize, 4); i += 4) {
						const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
						const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
						const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
						const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
						const __m256d re_part(_mm256_sub_pd(
							_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm3)));
						const __m256d im_part(_mm256_add_pd(
							_mm256_mul_pd(ymm2, ymm0), _mm256_mul_pd(ymm0, ymm3)));
						const __m256d sqrt_term1(_mm256_mul_pd(re_part, re_part));
						const __m256d sqrt_term2(_mm256_mul_pd(im_part, im_part));
						const __m256d mag_term(_mm256_sqrt_pd(_mm256_add_pd(sqrt_term1, sqrt_term2)));
						_mm256_store_pd(&out.data.m_Re[i], _mm256_div_pd(re_part, mag_part));
						_mm256_store_pd(&out.data.m_Im[i], _mm256_div_pd(im_part, mag_term));

					}
					for (; i != out.data.m_nsize; ++i) {
						const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) - (v1.data.m_Im[i] * v2.data.m_Im[i]);
						const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) + (v1.data.m_Re[i] * v2.data.m_Im[i]);
						const double mag = std::sqrt(re*re + im*im);
						out.data.m_Re[i] = re / mag;
						out.data.m_Im[i] = im / mag;
					}
				}
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmean_prod(std::complex<double> &mean,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2) {
				if (v1.data.m_nsize != v2.data.m_nsize) { return;}
#if defined _WIN64
				__declspec(align(64)) struct {
					double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					int32_t i; char pad[4];
				}ca;
#elif defined __linux
			    __attribute__((align(64))) struct {
                                       	double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					int32_t i; char pad[4];
			    }ca;
#endif
				for (ca.i = 0; i != ROUND_TO_FOUR(v1.data.m_nsize, 4); ca.i += 4) {
					const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
					const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
					const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
					const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
					const __m256d re_part(_mm256_sub_pd(_mm256_mul_pd(ymm0,ymm1),
										  _mm256_mul_pd(ymm2,ymm3)));
					const __m256d im_part(_mm256_add_pd(_mm256_mul_pd(ymm2,ymm3),
										 _mm256_mul_pd(ymm0,ymm3)));
					const __m256d t1(_mm256_add_pd(re_part, _mm256_permute_pd(re_part,0x05)));
					ca.sumre = _mm_cvtsd_f64(_mm_add_sd(_mm256_castpd256_pd128(t1), _mm256_extractf128_pd(t1,1)));
					ca.accre += ca.sumre;
					const __m256d t2(_mm256_add_pd(im_part, _mm256_permute_pd(im_part,0x05)));
					ca.sumim = _mm_cvtsd_f64(_mm_add_sd(_mm256_castpd256_pd128(t2), _mm256_extractf128_pd(t2,1)));
					ca.accim += ca.sumim;
				}
				for (; ca.i != v1.data.m_nsize; ++i) {
					ca.accre += (v1.data.m_Re[ca.i] * v2.data.m_Re[ca.i]) - (v1.data.m_Im[ca.i] * v2.data.m_Im[ca.i]);
					ca.accim += (v1.data.m_Im[ca.i] * v2.data.m_Re[ca.i]) + (v1.data.m_Re[ca.i] * v2.data.m_Im[ca.i]);
				}
				mean._Val[0] = ca.accre /= static_cast<double>(v1.data.m_nsize);
				mean._Val[1] = ca.accim /= static_cast<double>(v1.data.m_nsize);
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmean_quot(std::complex<double> &mean,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2) {
				if (v1.data.m_nsize != v2.data.m_nsize) { return;}
#if defined _WIN64
				__declspec(align(64)) struct {
					double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{ 0.0 },
					t1{ 0.0 }, t2{0.0};
					int32_t i;
					char pad[4];
				}ca;
#elif defined __linux
				__attribute__(align(64))) struct {
                                        double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{ 0.0 },
					t1{ 0.0 }, t2{0.0};
					int32_t i;
					char pad[4];
				}ca;
#endif
				for (ca.i = 0; i != ROUND_TO_FOUR(v1.data.m_nsize, 4); ca.i += 4) {
					const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
					const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
					const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
					const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
					const __m256d re_part(_mm256_sub_pd(_mm256_mul_pd(ymm0,ymm1),
										              _mm256_mul_pd(ymm2,ymm3)));
					const __m256d im_part(_mm256_add_pd(_mm256_mul_pd(ymm2,ymm1),
													  _mm256_mul_pd(ymm0,ymm3)));
					const __m256d den_term(_mm256_add_pd(_mm256_mul_pd(ymm1,ymm1),
													   _mm256_mul_pd(ymm3,ymm3)));
					const __m256d re_quot(_mm256_div_pd(re_part,den_term));
					const __m256d tmp1(_mm256_add_pd(re_quot, _mm256_permute_pd(re_quot,0x05)));
					ca.sumre = _mm_cvtsd_f64(_mm_add_sd(_mm256_castpd256_pd128(tmp1), _mm256_extractf128_pd(tmp1,1)));
					ca.accre += ca.sumre;
					const __m256d im_quot(_mm256_div_pd(im_part,den_term));
					const __m256d tmp2(_mm256_add_pd(im_quot, _mm256_permute_pd(im_quot,0x05)));
					ca.sumim = _mm_cvtsd_f64(_mm_add_sd(_mm256_castpd256_pd128(tmp2), _mm256_extractf128_pd(tmp2,1)));
					ca.accre += ca.sumim;
				}
				for (; ca.i != v1.data.m_nsize; ++ca.i) {
					const double re = (v1.data.m_Re[ca.i] * v2.data.m_Re[ca.i]) - (v1.data.m_Im[ca.i] * v2.data.m_Im[ca.i]);
					const double im = (v1.data.m_Im[ca.i] * v2.data.m_Re[ca.i]) + (v1.data.m_Re[ca.i] * v2.data.m_Im[ca.i]);
					const double den = (v2.data.m_Re[ca.i] * v2.data.m_Re[ca.i]) + (v2.data.m_Im[ca.i] * v2.data.m_Im[ca.i]);
					ca.t1 += re / den;
					ca.t2 += im / den;
				}
				ca.accre += ca.t1;
				ca.accim += ca.t2;
				mean._Val[0] = ca.accre / static_cast<double>(v1.data.m_nsize);
				mean._Val[1] = ca.accim / static_cast<double>(v1.data.m_nsize);
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cconj_prod(AVXVField1D &out,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2,
					  const bool do_nt_store) {
				if (v1.data.m_nsize != v2.data.m_nsize) { return;}
				int32_t i;
				if (do_nt_store) {
					for (i = 0; i != ROUND_TO_FOUR(out.data.m_nsize, 4); i += 4) {
						const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
						const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
						const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
						const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
						_mm256_stream_pd(&out.data.m_Re[i], _mm256_add_pd(
							_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2,ymm3)));
						_mm256_stream_pd(&out.data.m_Im[i], _mm256_sub_pd(
							_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0,ymm3)));
					}
					for (; i != v1.data.m_nsize; ++i) { // Cache storing remainder.
						out.data.m_Re[i] = (v1.data.m_Re[i] * v2.data.m_Re[i]) + (v1.data.m_Im[i] * v2.data.m_Im[i]);
						out.data.m_Im[i] = (v1.data.m_Im[i] * v2.data.m_Re[i]) - (v1.data.m_Re[i] * v2.data.m_Im[i]);
					}
				}
				else {
					for (i = 0; i != ROUND_TO_FOUR(out.data.m_nsize, 4); i += 4) {
						const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
						const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
						const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
						const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
						_mm256_store_pd(&out.data.m_Re[i], _mm256_add_pd(
							_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm3)));
						_mm256_store_pd(&out.data.m_Im[i], _mm256_sub_pd(
							_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0, ymm3)));
					}
					for (; i != v1.data.m_nsize; ++i) { 
						out.data.m_Re[i] = (v1.data.m_Re[i] * v2.data.m_Re[i]) + (v1.data.m_Im[i] * v2.data.m_Im[i]);
						out.data.m_Im[i] = (v1.data.m_Im[i] * v2.data.m_Re[i]) - (v1.data.m_Re[i] * v2.data.m_Im[i]);
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
				if (v1.data.m_nsize != v2.data.m_nsize) { return;}
				int32_t i;
				if (do_nt_store) {
					for (i = 0; i != ROUND_TO_FOUR(out.data.m_nsize, 4); i += 4) {
						const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
						const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
						const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
						const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
						const __m256d re_part(_mm256_add_pd(_mm256_mul_pd(ymm0,ymm1),
									  _mm256_mul_pd(ymm2,ymm3)));
						const __m256d im_part(_mm256_sub_pd(_mm256_mul_pd(ymm2,ymm1),
									 _mm256_mul_pd(ymm0,ymm3)));
						const __m256d mag_c1(_mm256_mul_pd(re_part,re_part));
						const __m256d mag_c2(_mm256_mul_pd(im_part,im_part));
						const __m256d vcmag(_mm256_sqrt_pd(_mm256_add_pd(mag_c1,mag_c2)));
						_mm256_stream_pd(&out.data.m_Re[i], _mm256_div_pd(re_part,vcmag));
						_mm256_stream_pd(&out.data.m_Im[i], _mm256_div_pd(im_part,vcmag));
					}
					for (; i != out.data.m_nsize; ++i) { // Cache storing remainder.
						const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) - (v1.data.m_Im[i] * v2.data.m_Im[i]);
						const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) + (v1.data.m_Re[i] * v2.data.m_Im[i]);
						const double mag = std::sqrt(re*re + im*im);
						out.data.m_Re[i] = re / mag;
						out.data.m_Im[i] = im / mag;
					}
				}
				else {
					for (i = 0; i != ROUND_TO_FOUR(out.data.m_nsize, 4); i += 4) {
						const __m256d ymm0(_mm256_load_pd(&v1.data.m_Re[i]));
						const __m256d ymm1(_mm256_load_pd(&v2.data.m_Re[i]));
						const __m256d ymm2(_mm256_load_pd(&v1.data.m_Im[i]));
						const __m256d ymm3(_mm256_load_pd(&v2.data.m_Im[i]));
						const __m256d re_part(_mm256_add_pd(_mm256_mul_pd(ymm0, ymm1),
							_mm256_mul_pd(ymm2, ymm3)));
						const __m256d im_part(_mm256_sub_pd(_mm256_mul_pd(ymm2, ymm1),
							_mm256_mul_pd(ymm0, ymm3)));
						const __m256d mag_c1(_mm256_mul_pd(re_part, re_part));
						const __m256d mag_c2(_mm256_mul_pd(im_part, im_part));
						const __m256d vcmag(_mm256_sqrt_pd(_mm256_add_pd(mag_c1, mag_c2)));
						_mm256_store_pd(&out.data.m_Re[i], _mm256_div_pd(re_part, vcmag));
						_mm256_store_pd(&out.data.m_Im[i], _mm256_div_pd(im_part, vcmag));
					}
					for (; i != out.data.m_nsize; ++i) { 
						const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) - (v1.data.m_Im[i] * v2.data.m_Im[i]);
						const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) + (v1.data.m_Re[i] * v2.data.m_Im[i]);
						const double mag = std::sqrt(re*re + im*im);
						out.data.m_Re[i] = re / mag;
						out.data.m_Im[i] = im / mag;
					}
				}
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmean_conjprod(std::complex<double> &mean,
					      const AVXVField1D &v1,
					      const AVXVField1D &v2) {
				if (v1.data.m_nsize != v2.data.m_nsize) { return;}
#if defined _WIN64
				__declspec(align(64)) struct {
					double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					
					int32_t i;
					char pad[4];
				}ca;
#elif defined __linux
				__attribute__((align(64))) struct {
                                        double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					
					int32_t i;
					char pad[4];
				}ca;
#endif
				__m256d re = _mm256_setzero_pd(), tmp1 = _mm256_setzero_pd();
				__m256d im = _mm256_setzero_pd(), tmp2 = _mm256_setzero_pd();
				for (ca.i = 0; i != ROUND_TO_FOUR(v1.data.m_nsize, 4); i += 4) {
					AVX_COMPLEX_ADDITION(re,v1,v2,ca.i,0)
					AVX_COMPLEX_SUBTRACTION(im,v1,v2,ca.i,0)
					AVX_HORIZONTAL_ADD(ca.sumre,re,tmp1)
					ca.accre += ca.sum_re;
					AVX_HORIZONTAL_ADD(ca.sumim,im,tmp2)
					ca.accim += ca.sumim;
				}	
				for (; i != v1.data.m_nsize; ++i) {
					const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) + (v1.data.m_Im[i] * v2.data.m_Im[i]);
					const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) - (v1.data.m_Re[i] * v2.data.m_Im[i]);
					ca.accre += re;
					ca.accim += im;
				}
				mean._Val[0] = ca.accre / static_cast<double>(v1.data.m_nsize);
				mean._Val[1] = ca.accim / static_cast<double>(v1.data.m_nsize);
			}	
				
			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_arith_mean(std::complex<double> &mean,
					  const AVXVField1D &v) {
#if defined _WIN64
				__declspec(align(64)) struct {
					double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					int32_t i;
					char pad[4];
				}ca;
				__declspec(align(64)) struct {
					__m256d re_part,im_part,tmp1,tmp2;
				}ca2;
#elif defined __linux
				__attribute__((align(64))) struct {
                                        double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
					int32_t i;
					char pad[4];
				}ca;
				__attribute__((align(64))) struct {
                                        __m256d re_part,im_part,tmp1,tmp2;
				}ca2;
#endif
				ca2.tmp1 = _mm256_setzero_pd();
				ca2.tmp2 = _mm256_setzero_pd();
				for (ca.i = 0; i != ROUND_TO_FOUR(v.data.m_nsize, 4); ca.i += 4) {
					ca2.re_part = _mm256_load_pd(&v.data.m_Re[i]);
					AVX_HORIZONTAL_ADD(ca.sumre,ca2.re_part,ca2.tmp1)
					ca.accre += ca.sumre;
					ca2.im_part = _mm256_load_pd(&v.data.m_Im[i]);
					AVX_HORIZONTAL_ADD(ca.sumim,ca2.im_part,tmp2)
					ca.accim += ca.sumim;
				}
				for (; ca.i != v.data.m_nsize; ++ca.i) {
					ca.accre += v.data.m_Re[i];
					ca.accim += v.data.m_Im[i];
				}
				mean._Val[0] = ca.accre / static_cast<double>(v.data.m_nsize);
				mean._Val[1] = ca.accim / static_cast<double>(v.data.m_nsize);
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cnormalize(AVXVField1D &out,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2,
					  const bool do_nt_store) {
				if (v1.data.m_nsize != v2.data.m_nsize) { return; }
				__m256d cvmag(_mm256_setzero_pd());
				int32_t i;
				if (do_nt_store) {

					for (i = 0, i != ROUND_TO_FOUR(out.data.m_nsize,4); i += 4) {
						AVX_COMPLEX_MAGNITUDE(cvmag, v1, v2, i, 0)
							const __m256d re_part(_mm256_load_pd(&v1.data.m_Re[i]));
						_mm256_stream_pd(&out.data.m_Re[i], _mm256_div_pd(re_part, cvmag));
						const __m256d im_part(_mm256_load_pd(&v1.data.m_Im[i]));
						_mm256_stream_pd(&out.data.m_Im[i], _mm256_div_pd(im_part, cvmag));
					}
					for (; i != v1.data.m_nsize; ++i) {
						const double cmag = std::sqrt(v1.data.m_Re[i] * v2.data.m_Re[i] + v1.data.m_Im[i] * v2.data.m_Im[i]);
						out.data.m_Re[i] = v1.data.m_Re[i] / cmag;
						out.data.m_Im[i] = v1.data.m_Im[i] / cmag;
					}
				}
				else {
					for (i = 0, i != ROUND_TO_FOUR(out.data.m_nsize, 4); i += 4) {
						AVX_COMPLEX_MAGNITUDE(cvmag, v, cv, i, 0)
							const __m256d re_part(_mm256_load_pd(&v1.data.m_Re[i]));
						_mm256_store_pd(&out.data.m_Re[i], _mm256_div_pd(re_part, cvmag));
						    const __m256d im_part(_mm256_load_pd(&v1.data.m_Im[i]));
						_mm256_store_pd(&out.data.m_Im[i], _mm256_div_pd(im_part, cvmag));
					}
					for (; i != v1.data.m_nsize; ++i) {
						const double cmag = std::sqrt(v1.data.m_Re[i] * v2.data.m_Re[i] + v1.data.m_Im[i] * v2.data.m_Im[i]);
						out.data.m_Re[i] = v1.data.m_Re[i] / cmag;
						out.data.m_Im[i] = v1.data.m_Im[i] / cmag;
					}
				}
			}

			template<class AVXVField1D,
			typename = std::enable_if<
			std::is_class<AVXVField1D>::value,void>::type>
			avx256_cmagnitude(AVXVField1D &out,
					  const AVXVField1D &v1,
					  const AVXVField1D &v2) {
				if (v1.data.m_nsize != v2.data.m_nsize) { return;}
				__m256d vcmag(_mm256_setzero_pd());
				int32_t i;
				for (i = 0; i != ROUND_TO_FOUR(out.data.m_nsize, 4); i += 4) {
					    AVX_COMPLEX_MAGNITUDE(vcmag,v1,v2,i,0)
						_mm256_store_pd(&out.data.m_Re[i], vcmag);
				}
				for (; i != out.data.m_nsize; ++i) {
					out.data.m_Re[i] = std::sqrt(v1.data.m_Re[i] * v2.data.m_Re[i] + v1.data.m_Im[i] * v2.data.m_Im[i]);
				}
			}
			
	}
}





#endif /*__GMS_AVXCOMPLEX_COMMON_H__*/
