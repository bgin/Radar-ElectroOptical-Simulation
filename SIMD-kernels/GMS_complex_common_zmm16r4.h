
#ifndef __GMS_COMPLEX_COMMON_ZMM16R4_H__
#define __GMS_COMPLEX_COMMON_ZMM16R4_H__ 121020191814

namespace file_info {


	const unsigned int GMS_COMPLEX_COMMON_ZMM16R4_MAJOR = 1U;

	const unsigned int GMS_COMPLEX_COMMON_ZMM16R4_MINOR = 0U;

	const unsigned int GMS_COMPLEX_COMMON_ZMM16R4_MICRO = 0U;

	const unsigned int GMS_COMPLEX_COMMON_ZMM16R4_FULLVER = 
		1000U*GMS_COMPLEX_COMMON_ZMM16R4_MAJOR + 100U*GMS_COMPLEX_COMMON_ZMM16R4_MINOR + 10U*GMS_COMPLEX_COMMON_ZMM16R4_MICRO;

	const char * const GMS_COMPLEX_COMMON_ZMM16R4_CREATE_DATE = "12-10-2019 18:14 +00200 (SAT 12 OCT 2019 GMT+2)";

	const char * const GMS_COMPLEX_COMMON_ZMM16R4_BUILD_DATE = __DATE__ ":" __TIME__;

	const char * const GMS_COMPLEX_COMMON_ZMM16R4_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const GMS_COMPLEX_COMMON_ZMM16R4_SYNOPSIS = "Common procedures for AVX512 complex field (1D) classes.";

}


#include <cstdint>
#include <type_traits>
#include <immintrin.h>
#include <complex>
#include "GMS_config.h"

namespace gms {
	namespace math {

#if !defined (AVX512_COMPLEX_ADDITION)
#define AVX512_COMPLEX_ADDITION(out,v1,v2,idx,off) \
	(out) = _mm512_add_ps(_mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Re[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Re[(idx)+(off)])), _mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Im[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Im[(idx)+(off)])));
#endif

#if !defined (AVX512_COMPLEX_SUBTRACTION)
#define AVX512_COMPLEX_SUBTRACTION(out,v1,v2,idx,off) \
	(out) = _mm512_sub_ps(_mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Im[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Re[(idx)+(off)])), _mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Re[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Im[(idx)+(off)])));
#endif

		// Warning macro parameter v2 must be an exact copy
		// of parameter v1. This should done by calling class (AVX512VComplex1D)
		// Move Constructor.
#if !defined (AVX512_COMPLEX_MAGNITUDE)
#define AVX512_COMPLEX_MAGNITUDE(out,v1,v2,idx,off) \
	(out) = _mm512_sqrt_ps(_mm512_add_ps(_mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Re[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Re[(idx)+(off)])), _mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Im[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Im[(idx)+(off)]))));
#endif

			template<class AVX512VField1D,
				    typename = std::enable_if<
					 std::is_class<AVX512VField1D>::value,void>::type>
					  avx512_cnormalize_prod(AVX512VField1D &out,
								 const AVX512VField1D &v1,
					                         const AVX512VField1D &v2,
					                         const bool do_nt_stream) {
				 if (v1.data.m_nsize != v2.data.m_nsize) {return;}
				 int64_t i;
				 if (do_nt_stream) {// Incurs branch penalty here circa. 50% upon first encounters.
					 for (i = 0; i != ROUND_TO_SIXTEEN(out.data.m_nsize, 16); i += 16) {

						 const __m512 zmm0 = _mm512_load_ps(&v1.data.m_Re[i]);
						 const __m512 zmm1 = _mm512_load_ps(&v2.data.m_Re[i]);
						 const __m512 zmm2 = _mm512_load_ps(&v1.data.m_Im[i]);
						 const __m512 zmm3 = _mm512_load_ps(&v2.data.m_Im[i]);

						 const __m512 re_part = _mm512_sub_ps(
							 _mm512_mul_ps(zmm0, zmm1), _mm512_mul_ps(zmm2,zmm3));
						 const __m512 im_part = _mm512_add_ps(
							 _mm512_mul_ps(zmm2, zmm1), _mm512_mul_ps(zmm0,zmm3));
						 const __m512 sqrt_term1 = _mm512_mul_ps(re_part,re_part);
						 const __m512 sqrt_term2 = _mm512_mul_ps(im_part,im_part);
						 const __m512 mag_term = _mm512_sqrt_ps(_mm512_add_ps(sqrt_term1,sqrt_term2));
						 _mm512_stream_ps(&out.data.m_Re[i], _mm512_div_ps(re_part,mag_term));
						 _mm512_stream_ps(&out.data.m_Im[i], _mm512_div_ps(im_part,mag_term));
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
					    
					 for (i = 0; i != ROUND_TO_SIXTEEN(out.data.m_nsize, 16); i += 16) {

						 const __m512 zmm0 = _mm512_load_ps(&v1.data.m_Re[i]);
						 const __m512 zmm1 = _mm512_load_ps(&v2.data.m_Re[i]);
						 const __m512 zmm2 = _mm512_load_ps(&v1.data.m_Im[i]);
						 const __m512 zmm2 = _mm512_load_ps(&v2.data.m_Im[i]);

						 const __m512 re_part = _mm512_sub_ps(
							 _mm512_mul_ps(zmm0, zmm1), _mm512_mul_ps(zmm2, zmm3));
						 const __m512 im_part = _mm512_add_ps(
							 _mm512_mul_ps(zmm2, zmm1), _mm512_mul_ps(zmm0, zmm3));
						 const __m512 sqrt_term1 = _mm512_mul_ps(re_part, re_part);
						 const __m512 sqrt_term2 = _mm512_mul_ps(im_part, im_part);
						 const __m512 mag_term = _mm512_sqrt_ps(_mm512_add_ps(sqrt_term1, sqrt_term2));
						 _mm512_store_ps(&out.data.m_Re[i], _mm512_div_ps(re_part, mag_term));
						 _mm512_store_ps(&out.data.m_Im[i], _mm512_div_ps(im_part, mag_term));
					 }

					 for (; i != out.data.m_nsize; ++i) {
						 const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) - (v1.data.m_Im[i] * v2.data.m_Im[i]);
						 const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) + (v1.data.m_Re[i] * v2.data.m_Im[i]);
						 const double mag = std::sqrt((re*re) + (im*im));
						 out.data.m_Re[i] = re / mag;
						 out.data.m_Im[i] = im / mag;
					 }
			  }
		}

		template<class AVX512VField1D,
			     typename = std::enable_if<
				 std::is_class<AVX512VField1D>::value,void>::type>
				 avx512_cmean_prod(std::complex<double> &mean,
						   const AVX512VField1D &v1,
						   const AVX512VField1D &v2) {
					 if (v1.data.m_nsize != v2.data.m_nsize) { return;}
				   
				 	__ATTR_ALIGN__(64) struct {
                                                 double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{0.0};
						 int64_t i;
					 } ca;
				 
					 for (ca.i = 0; ca.i != ROUND_TO_SIXTEEN(v1.data.m_nsize, 16); ca.i += 16) {
						 
						 const __m512 zmm0(_mm512_load_ps(&v1.data.m_Re[i]));
						 const __m512 zmm1(_mm512_load_ps(&v2.data.m_Re[i]));
						 const __m512 zmm2(_mm512_load_ps(&v1.data.m_Im[i]));
						 const __m512 zmm3(_mm512_load_ps(&v2.data.m_Im[i]));
						 const __m512 re_part(_mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1)
											   _mm512_mul_ps(zmm2,zmm3)));
						 const __m512 im_part(_mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
											   _mm512_mul_ps(zmm0,zmm3)));
						 ca.sumre = _mm512_reduce_ps(re_part);
						 ca.accre += ca.sumre;
						 ca.sumim = _mm512_reduce_ps(im_part);
						 ca.accim += ca.sumim;
				}
					 for (; ca.i != v1.data.m_nsize; ++i) {
						   ca.accre += (v1.data.m_Re[ca.i] * v2.data.m_Re[ca.i]) - (v1.data.m_Im[ca.i] * v2.data.m_Im[ca.i]);
						   ca.accim += (v1.data.m_Im[ca.i] * v2.data.m_Re[ca.i]) + (v1.data.m_Re[ca.i] * v2.data.m_Im[ca.i]);
					 }
					 mean._Val[0] = ca.accre /= static_cast<double>(v1.data.m_nsize);
					 mean._Val[1] = ca.accim /= static_cast<double>(v1.data.m_nsize);
		 }

		 template<class AVX512VField1D,
				  typename = std::enable_if<
				  std::is_class<AVX512VField1D>::value,void>::type>
				  avx512_cmean_quot(std::complex<double> &mean,
						    const AVX512VField1D &v1,
						    const AVX512VField1D &v2) {
				  if (v1.data.m_nsize != v2.data.m_nsize) { return;}
			    
				  __ATTR_ALIGN__(64) struct {
					  double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{ 0.0 },
					  t1{ 0.0 }, t2{0.0};
					  int64_t i;
				  }ca;
			  
				  for (ca.i = 0; ca.i != ROUND_TO_SIXTEEN(v1.data.m_nsize, 16); ca.i += 16) {
						
						const __m512 zmm0(_mm512_load_ps(&v1.data.m_Re[i]));
						const __m512 zmm1(_mm512_load_ps(&v2.data.m_Re[i]));
						const __m512 zmm2(_mm512_load_ps(&v1.data.m_Im[i]));
						const __m512 zmm3(_mm512_load_ps(&v2.data.m_Im[i]));
						const __m512 re_part(_mm512_sub_ps(_mm512_mul_ps(zmm0, zmm1),
											  _mm512_mul_ps(zmm2, zmm3)));
						const __m512 im_part(_mm512_add_ps(_mm512_mul_ps(zmm2, zmm1),
											  _mm512_mul_ps(zmm0, zmm3)));
						const __m512 den_term(_mm512_add_ps(_mm512_mul_ps(zmm1, zmm1),
											  _mm512_mul_ps(zmm3, zmm3)));
						const __m512 re_quot(_mm512_div_ps(re_part,den_term));
						ca.sumre = _mm512_reduce_ps(re_quot);
						ca.accre += ca.sumre;
						const __m512 im_quot(_mm512_div_ps(im_part,den_term));
						ca.sumim = _mm512_reduce_ps(im_quot);
						ca.accim += ca.sumim;
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

		template<class AVX512VField1D,
			    typename = std::enable_if<
				std::is_class<AVX512VField1D>::value,void>::type>
				avx512_cconj_prod(AVX512VField1D &out,
						  const AVX512VField1D &v1,
						  const AVX512VField1D &v2,
						  const bool do_nt_store) {
				if (v1.data.m_nsize != v1.data.m_nsize) {return;}
				int64_t i;
				if (do_nt_store) { // Perform streaming-store

					for (i = 0; i != ROUND_TO_SIXTEEN(out.data.m_nsize, 16); i += 16) {
						 
						const __m512 zmm0(_mm512_load_ps(&v1.data.m_Re[i]));
						const __m512 zmm1(_mm512_load_ps(&v2.data.m_Re[i]));
						const __m512 zmm2(_mm512_load_ps(&v1.data.m_Im[i]));
						const __m512 zmm3(_mm512_load_ps(&v2.data.m_Im[i]));
						_mm512_stream_ps(&out.data.m_Re[i], _mm512_add_ps(
							            _mm512_mul_ps(zmm0, zmm1), _mm512_mul_ps(zmm2, zmm3)));
						_mm512_stream_ps(&out.data.m_Im[i], _mm512_sub_ps(
									    _mm512_mul_ps(zmm2, zmm1), _mm512_mul_ps(zmm0, zmm3)));
					}

					for (; i != v1.data.m_nsize; ++i) { // Cache storing remainder.
						out.data.m_Re[i] = (v1.data.m_Re[i] * v2.data.m_Re[i]) + (v1.data.m_Im[i] * v2.data.m_Im[i]);
						out.data.m_Im[i] = (v1.data.m_Im[i] * v2.data.m_Re[i]) - (v1.data.m_Re[i] * v2.data.m_Im[i]);
					}
				}
				else {
					
					for (i = 0; i != ROUND_TO_SIXTEEN(out.data.m_nsize, 16); i += 16) {

						const __m512 zmm0(_mm512_load_ps(&v1.data.m_Re[i]));
						const __m512 zmm1(_mm512_load_ps(&v2.data.m_Re[i]));
						const __m512 zmm2(_mm512_load_ps(&v1.data.m_Im[i]));
						const __m512 zmm3(_mm512_load_ps(&v2.data.m_Im[i]));
						_mm512_store_ps(&out.data.m_Re[i], _mm512_add_ps(
							_mm512_mul_ps(zmm0, zmm1), _mm512_mul_ps(zmm2, zmm3)));
						_mm512_store_ps(&out.data.m_Im[i], _mm512_sub_ps(
							_mm512_mul_ps(zmm2, zmm1), _mm512_mul_ps(zmm0, zmm3)));
					}

					for (; i != v1.data.m_nsize; ++i) { 
						out.data.m_Re[i] = (v1.data.m_Re[i] * v2.data.m_Re[i]) + (v1.data.m_Im[i] * v2.data.m_Im[i]);
						out.data.m_Im[i] = (v1.data.m_Im[i] * v2.data.m_Re[i]) - (v1.data.m_Re[i] * v2.data.m_Im[i]);
					}
				}
	        }

			template<class AVX512VField1D,
					 typename = std::enable_if<
					 std::is_class<AVX512VField1D>::value,void>::type>
					 avx512_cnorm_conjprod(AVX512VField1D &out,
							       const AVX512VField1D &v1,
							       const AVX512VField1D &v2,
							       const bool do_nt_store) {
					 if (v1.data.m_nsize != v2.data.m_nsize) { return;}
					 int64_t i;
					 if (do_nt_store) {
							
						 for (i = 0; i != ROUND_TO_SIXTEEN(out.data.m_nsize, 16); i += 16) {
							
							 const __m512 zmm0(_mm512_load_ps(&v1.data.m_Re[i]));
							 const __m512 zmm1(_mm512_load_ps(&v2.data.m_Re[i]));
							 const __m512 zmm2(_mm512_load_ps(&v1.data.m_Im[i]));
							 const __m512 zmm3(_mm512_load_ps(&v2.data.m_Im[i]));
							 const __m512 re_part(_mm512_add_ps(_mm512_mul_ps(zmm0, zmm1),
													   _mm512_mul_ps(zmm2, zmm3)));
							 const __m512 im_part(_mm512_sub_ps(_mm512_mul_ps(zmm2, zmm1),
													   _mm512_mul_ps(zmm0, zmm3)));
							 const __m512 mag_c1(_mm512_mul_ps(re_part, re_part));
							 const __m512 mag_c2(_mm512_mul_ps(im_part, im_part));
							 const __m512 vcmag(_mm512_sqrt_ps(_mm512_add_ps(mag_c1, mag_c2)));
							 _mm512_stream_ps(&out.data.m_Re[i], _mm512_div_ps(re_part,vcmag));
							 _mm512_stream_ps(&out.data.m_Im[i], _mm512_div_ps(im_part,vcmag));
						 }

						 for (; i != out.data.m_nsize; ++i) { // Cache storing remainder.
							 const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) - (v1.data.m_Im[i] * v2.data.m_Im[i]);
							 const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) + (v1.data.m_Re[i] * v2.data.m_Im[i]);
							 const double mag = std::sqrt((re*re) + (im*im));
							 out.data.m_Re[i] = re / mag;
							 out.data.m_Im[i] = im / mag;
						 }
					 }
					 else {
							
						 for (i = 0LL; i != ROUND_TO_SIXTEEN(out.data.m_nsize, 8LL); i += 8LL) {

							 const __m512 zmm0(_mm512_load_ps(&v1.data.m_Re[i]));
							 const __m512 zmm1(_mm512_load_ps(&v2.data.m_Re[i]));
							 const __m512 zmm2(_mm512_load_ps(&v1.data.m_Im[i]));
							 const __m512 zmm3(_mm512_load_ps(&v2.data.m_Im[i]));
							 const __m512 re_part(_mm512_add_ps(_mm512_mul_ps(zmm0, zmm1),
								 _mm512_mul_ps(zmm2, zmm3)));
							 const __m512 im_part(_mm512_sub_ps(_mm512_mul_ps(zmm2, zmm1),
								 _mm512_mul_ps(zmm0, zmm3)));
							 const __m512 mag_c1(_mm512_mul_ps(re_part, re_part));
							 const __m512 mag_c2(_mm512_mul_ps(im_part, im_part));
							 const __m512 vcmag(_mm512_sqrt_ps(_mm512_add_ps(mag_c1, mag_c2)));
							 _mm512_store_ps(&out.data.m_Re[i], _mm512_div_ps(re_part, vcmag));
							 _mm512_store_ps(&out.data.m_Im[i], _mm512_div_ps(im_part, vcmag));
						 }

						 for (; i != out.data.m_nsize; ++i) { 
							 const double re = (v1.data.m_Re[i] * v2.data.m_Re[i]) - (v1.data.m_Im[i] * v2.data.m_Im[i]);
							 const double im = (v1.data.m_Im[i] * v2.data.m_Re[i]) + (v1.data.m_Re[i] * v2.data.m_Im[i]);
							 const double mag = std::sqrt((re*re) + (im*im));
							 out.data.m_Re[i] = re / mag;
							 out.data.m_Im[i] = im / mag;
						 }
					 }

			    }
		  
		  template<class AVX512VField1D,
				   typename = std::enable_if<
				   std::is_class<AVX512VField1D>::value,void>::type>
				   avx512_cmean_conjprod(std::complex<double> &mean,
							 const AVX512VField1D &v1,
							 const AVX512VField1D &v2) {
					   if (v1.data.m_nsize != v2.data.m_nsize) { return;}
				    
					__ATTR_ALIGN__(64) struct {
                                                    double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{ 0.0 };
						    int64_t i;
					   }ca;
					   __m512 re = _mm512_setzero_ps();
					   __m512 im = _mm512_setzero_ps();
					   for (ca.i = 0; ca.i != ROUND_TO_SIXTEEN(v1.data.m_nsize, 16); i += 16) {
						    
						   AVX512_COMPLEX_ADDITION(re,v1,v2,ca.i,0)
						   AVX512_COMPLEX_SUBTRACTION(im,v1,v2,ca.i,0)
						   ca.sumre = _mm512_reduce_ps(re);
						   ca.accre += ca.sumre;
						   ca.sumim = _mm512_reduce_ps(im);
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

			
			template<class AVX512VField1D, 
			         typename = std::enable_if<
					 std::is_class<AVX512VField1D>::value,void>::type>
					 avx512_arith_mean(std::complex<double> &mean,
							   const AVX512VField1D &v) {
			                 
						 __ATTR_ALIGN__(64) struct {
                                                         double sumre{ 0.0 }, sumim{ 0.0 }, accre{ 0.0 }, accim{ 0.0 };
							 int64_t i;
						 }ca;
					 
						 for (ca.i = 0; ca.i != ROUND_TO_SIXTEEN(v.data.m_nsize, 16); ca.i += 16) {
							 const __m512 re_part(_mm512_load_ps(&v.data.m_Re[i]));
							 ca.sumre = _mm512_reduce_ps(re_part);
							 ca.accre += ca.sumre;
							 const __m512 im_part(_mm512_load_ps(&v.data.m_Im[i]));
							 ca.sumim = _mm512_reduce_ps(im_part);
							 ca.accim += ca.sumim;
						 }
						 for (; ca.i != v.data.m_nsize; ++ca.i) {
							 ca.accre += v.data.m_Re[i];
							 ca.accim += v.data.m_Im[i];
						 }
						 mean._Val[0] = ca.accre / static_cast<double>(v.data.m_nsize);
						 mean._Val[1] = ca.accim / static_cast<double>(v.data.m_nsize);
			}

			template<class AVX512VField1D,
					 typename = std::enable_if<
								std::is_class<AVX512VField1D>::value,void>::type>
								avx512_cnormalize(AVX512VField1D &out,
										  const AVX512VField1D &v,
										  const AVX512VField1D &cv,
										  const bool do_nt_store) {
						if (v.data.m_nsize != cv.data.m_nsize) { return;}
						__m512 cvmag(_mm512_setzero_ps());
						int64_t i;
						if (do_nt_store) {
								
							for (i = 0, i != ROUND_TO_SIXTEEN(out.data.m_nsize, 16); i += 16) {
								AVX512_COMPLEX_MAGNITUDE(cvmag,v,cv,i,0)
								const __m512 re_part(_mm512_load_ps(&v.data.m_Re[i]));
								_mm512_stream_ps(&out.data.m_Re[i], _mm512_div_ps(re_part,cvmag));
								const __m512 im_part(_mm512_load_ps(&v.data.m_Im[i]));
								_mm512_stream_ps(&out.data.m_Im[i], _mm512_div_ps(im_part,cvmag));
							}
							for (; i != v.data.m_nsize; ++i) {
								const double cmag = std::sqrt(v.data.m_Re[i] * cv.data.m_Re[i] + v.data.m_Im[i] * cv.data.m_Im[i]);
								out.data.m_Re[i] = v.data.m_Re[i] / cmag;
								out.data.m_Im[i] = v.data.m_Im[i] / cmag;
							}
						}
						else {
							for (i = 0LL, i != ROUND_TO_SIXTEEN(out.data.m_nsize, 8LL); i += 8LL) {
								AVX512_COMPLEX_MAGNITUDE(cvmag, v, cv, i, 0LL)
									const __m512 re_part(_mm512_load_ps(&v.data.m_Re[i]));
								_mm512_store_ps(&out.data.m_Re[i], _mm512_div_ps(re_part, cvmag));
								const __m512 im_part(_mm512_load_ps(&v.data.m_Im[i]));
								_mm512_store_ps(&out.data.m_Im[i], _mm512_div_ps(im_part, cvmag));
							}
							for (; i != v.data.m_nsize; ++i) {
								const double cmag = std::sqrt(v.data.m_Re[i] * cv.data.m_Re[i] + v.data.m_Im[i] * cv.data.m_Im[i]);
								out.data.m_Re[i] = v.data.m_Re[i] / cmag;
								out.data.m_Im[i] = v.data.m_Im[i] / cmag;
							}
						}
				}

				template<class AVX512VField1D,
						 typename = std::enable_if<
						 std::is_class<AVX512VField1D>::value,void>::type>
						void avx512_cmagnitude(double * __restrict vmag, 
								       const AVX512VField1D &v,
								       const AVX512VField1D &cv) {
						if (v.data.m_nsize != cv.data.m_nsize) { return;}
						__m512 vcmag(_mm512_setzero_ps());
						int64_t i;
						for (i = 0; i != ROUND_TO_SIXTEEN(v.data.m_nsize, 16); i += 16) {

							AVX512_COMPLEX_MAGNITUDE(vcmag, v, cv, i, 0)
								_mm512_store_ps(&vmag[i], vcmag);
						}
						for (; i != v.data.m_nsize; ++i)
							vmag[i] = std::sqrt(v.data.m_Re[i] * cv.data.m_Re[i] + v.data.m_Im[i] * cv.data.m_Im[i]);
			   }

			

	}
}

 


#endif /*__GMS_COMPLEX_COMMON_ZMM16R4_H__*/
