
#include <iostream>
#include <iomanip>
#include <limits>
#include <cfenv>
#include <cmath>
#include "GMS_helpers.h"




bool  check_null_ptr_2args(const float* __restrict src_ptr,
                           float* __restrict dst_ptr,
			   const char * func_name) {

	if ((src_ptr == NULL) || (dst_ptr == NULL)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in src_ptr=" << std::hex << "0x" << src_ptr  << "\n"
			<< "address held in dst_ptr=" << std::hex << "0x" << dst_ptr  << "\n"
			<< "address   of    src_ptr=" << std::hex << "0x" << &src_ptr << "\n"
			<< "address   of    dst_ptr=" << std::hex << "0x" << &dst_ptr << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
		  return (true);
		
	}
	return (false);
}


bool  check_null_ptr_2args(const double* __restrict src_ptr,
                           double* __restrict dst_ptr,
			   const char * func_name) {

	if ((src_ptr == NULL) || (dst_ptr == NULL)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in  src_ptr=" << std::hex << "0x" << src_ptr << "\n"
			<< "address held in  dst_ptr=" << std::hex << "0x" << dst_ptr << "\n"
			<< "address   of     src_ptr="  << std::hex << "0x" << &src_ptr << "\n"
			<< "address   of     dst_ptr="  << std::hex << "0x" << &dst_ptr << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
		  return (true);
	}
	return (false);
}


bool  check_null_ptr_3args(const float* __restrict src1_ptr,
			   const float* src2_ptr,
	                   float* __restrict dst_ptr,
			   const char * func_name       ) {

	if ((src1_ptr == NULL) || (src2_ptr == NULL) || (dst_ptr == NULL)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in  src1_ptr=" << std::hex << "0x" << src1_ptr << "\n"
			<< "address held in  src2_ptr=" << std::hex << "0x" << src2_ptr << "\n"
			<< "address held in  dst_ptr =" << std::hex << "0x" << dst_ptr  << "\n"
			<< "address   of     src1_ptr=" << std::hex << "0x" << &src1_ptr << "\n"
			<< "address   of     src2_ptr=" << std::hex << "0x" << &src2_ptr << "\n"
			<< "address   of     dst2_ptr=" << std::hex << "0x" << &dst_ptr  << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
		  return (true);
	}
	return (false);
}


bool  check_null_ptr_3args(const double* __restrict src1_ptr,
			   const double* __restrict src2_ptr,
	                   double* __restrict dst_ptr,
			   const char * func_name)   {

	if ((src1_ptr == NULL) || (src2_ptr == NULL) || (dst_ptr == NULL)) {
#if (GMS_DEBUG_ON) == 1    
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in  src1_ptr=" << std::hex << "0x" << src1_ptr << "\n"
			<< "address held in  src2_ptr=" << std::hex << "0x" << src2_ptr << "\n"
			<< "address held in  dst_ptr =" << std::hex << "0x" << dst_ptr << "\n"
			<< "address   of     src1_ptr=" << std::hex << "0x" << &src1_ptr << "\n"
			<< "address   of     src2_ptr=" << std::hex << "0x" << &src2_ptr << "\n"
			<< "address   of     dst2_ptr=" << std::hex << "0x" << &dst_ptr << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
		  return (true);
	}
	return (false);
}


bool  check_null_ptr_2args(const int* __restrict src_ptr,
			   int* __restrict dst_ptr,
			   const char * func_name) {

	if ((src_ptr == NULL) || (dst_ptr == NULL)){
#if (GMS_DEBUG_ON) == 1	
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in  src_ptr=" << std::hex << "0x" << src_ptr << "\n"
			<< "address held in  dst_ptr=" << std::hex << "0x" << dst_ptr << "\n"
			<< "address   of     src_ptr=" << std::hex << "0x" << &src_ptr << "\n"
			<< "address   of     dst_ptr=" << std::hex << "0x" << &dst_ptr << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
                    return (true);
	}
       return (false);
}


bool check_null_ptr_2args(const long long* __restrict src_ptr,
	                   long long* __restrict dst_ptr,
			   const char * func_name) {

	if ((src_ptr == NULL) || (dst_ptr == NULL)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in  src_ptr=" << std::hex << "0x" << src_ptr << "\n"
			<< "address held in  dst_ptr=" << std::hex << "0x" << dst_ptr << "\n"
			<< "address   of     src_ptr=" << std::hex << "0x" << &src_ptr << "\n"
			<< "address   of     dst_ptr=" << std::hex << "0x" << &dst_ptr << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
                     return (true);
	}
       return (false);
}


bool  check_null_ptr_3args(const int* __restrict src1_ptr,
                           const int* __restrict src2_ptr,
			   int* __restrict dst_ptr,
			   const char * func_name              ) {

	if ((src1_ptr == NULL) || (src2_ptr == NULL) || (dst_ptr == NULL)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in  src1_ptr=" << std::hex << "0x" << src1_ptr << "\n"
			<< "address held in  src2_ptr=" << std::hex << "0x" << src2_ptr << "\n"
			<< "address held in  dst_ptr =" << std::hex << "0x" << dst_ptr << "\n"
			<< "address   of     src1_ptr=" << std::hex << "0x" << &src1_ptr << "\n"
			<< "address   of     src2_ptr=" << std::hex << "0x" << &src2_ptr << "\n"
			<< "address   of     dst2_ptr=" << std::hex << "0x" << &dst_ptr << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
                     return (true);
	}
       return (false);
}


bool  check_null_ptr_3args(const long long* __restrict src1_ptr,
                           const long long* __restrict src2_ptr,
			   long long* __restrict dst_ptr,
			   const char * func_name                    ) {

	if ((src1_ptr == NULL) || (src2_ptr == NULL) || (dst_ptr == NULL)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Invalid Pointer in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "address held in  src1_ptr=" << std::hex << "0x" << src1_ptr << "\n"
			<< "address held in  src2_ptr=" << std::hex << "0x" << src2_ptr << "\n"
			<< "address held in  dst_ptr =" << std::hex << "0x" << dst_ptr << "\n"
			<< "address   of     src1_ptr=" << std::hex << "0x" << &src1_ptr << "\n"
			<< "address   of     src2_ptr=" << std::hex << "0x" << &src2_ptr << "\n"
			<< "address   of     dst2_ptr=" << std::hex << "0x" << &dst_ptr << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
                        return (true);
	}
       return (false);
}


bool  check_arrays_size( const int src_len,
                         const int dst_len,
			 const char * func_name) {

	if ((src_len != dst_len) || (src_len <= 0 || dst_len <= 0)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Array size mismatch in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n"
			<< "value of src_len=" << src_len << "\n"
			<< "value of dst_len=" << dst_len << "\n"
			<< "*****ERROR-DETAILS***** \n"
#endif
                        return (true);
	}
       return (false);
}








void  dump_timing_stats(const unsigned __int64 start_clock, const int src_len,
	                const int unrolled_iters,const int dst_len, const int elem_size,
	                const unsigned __int64 stop_clock, const int tsc_aux, const char * intrin_name) {
	
	std::cout << "Crude approximation of " << intrin_name << "\n"
		<< "--------------- Dumping Stats ---------------- \n"
		<< "----------------------------------------------- \n"
		<< "Loop iterations:          " << src_len << ".\n"
		<< "Unrolled loop iterations: " <<  unrolled_iters << ".\n"
		<< "src array size:           " << static_cast<double>(src_len * elem_size) / 1024.0 << "KiB.\n"
		<< "dst array size:           " << static_cast<double>(dst_len * elem_size) / 1024.0 << "KiB.\n"
		<< "IA32_TSC_AUX MSR:         " << std::hex << "0x" << tsc_aux << ".\n"
		<< "TSC on entry:             " << start_clock << " TSC-cycles. \n"
		<< "TSC on exit:              " << stop_clock << " TSC-cycles.  \n"
		<< "TSC Delta:                " << stop_clock - start_clock << "cycles. \n"
		<< "TSC per loop cycle:       " << (static_cast<double>(stop_clock - start_clock) / (src_len / elem_size)) << " TSC-cycles.\n"
		<< "--------------- End of Dump ------------------- \n";

}


void  is_domain_f32_lt_zero(const float* __restrict input,
			    const int length,
			    int* err) {

	 __m256 v_zmask(_mm256_setzero_ps());
	 const __m256 v_zerof32(_mm256_set1_ps(0.F));
	 constexpr int f32_off{8};
         *err = 0;
	 for (int i{ 0 }; i != length; i += f32_off) {
		  
		 v_zmask = _mm256_cmp_ps(v_zerof32, _mm256_loadu_ps(&input[i]), _CMP_LT_OQ);
		 if (_mm256_testc_ps(v_zmask, _mm256_setzero_ps())) {
			 *err = 1;
			 break;
		 }
	 }

}


void  is_domain_f32_lt_zero(const float* __restrict in1,
                            const float* __restrict in2,
	                    const int length,
			    int* err ) {

	__m256 v_zmask1(_mm256_setzero_ps());
	__m256 v_zmask2(_mm256_setzero_ps());
	const __m256 v_zerof32(_mm256_set1_ps(0.F));
	constexpr int f32_off{8};
        *err = 0;  

	for (int i{ 0 }; i != length; i += f32_off) {

		v_zmask1 = _mm256_cmp_ps(v_zerof32, _mm256_loadu_ps(&in1[i]), _CMP_LT_OQ);
		v_zmask2 = _mm256_cmp_ps(v_zerof32, _mm256_loadu_ps(&in2[i]), _CMP_LT_OQ);
		if ((_mm256_testc_ps(v_zmask1, _mm256_setzero_ps())) ||
			                (_mm256_testc_ps(v_zmask2, _mm256_setzero_ps()))) {
				*err = 1;
				break;
		}
	}
}


void  is_domain_f64_lt_zero(const double* __restrict in,
	                    const int length,
			   int* err) {

	__m256d v_zmask(_mm256_setzero_pd());
	const __m256d v_zerof64(_mm256_set1_pd(0.0));
	constexpr int f64_off{4};
        *err = 0;

	for (int i{ 0 }; i != length; i += f64_off) {

		v_zmask = _mm256_cmp_pd(v_zerof64, _mm256_loadu_pd(&in[i]), _CMP_LT_OQ);
		if (_mm256_testc_pd(v_zmask, _mm256_setzero_pd())) {
			*err = 1;
			break;
		}
	}
}


void  is_domain_f64_lt_zero(const double* __restrict in1,
                            const double* __restrict in2,
			    const int length,
			    int* err                ) {
	
	__m256d v_zmask1(_mm256_setzero_pd());
	__m256d v_zmask2(_mm256_setzero_pd());
	const __m256d v_zerof64(_mm256_set1_pd(0.0));
	constexpr int f64_off{4};
        *err = 0;

	for (int i{ 0 }; i != length; i += f64_off) {

		v_zmask1 = _mm256_cmp_pd(v_zerof64, _mm256_loadu_pd(&in1[i]), _CMP_LT_OQ);
		v_zmask2 = _mm256_cmp_pd(v_zerof64, _mm256_loadu_pd(&in2[i]), _CMP_LT_OQ);
		if ((_mm256_testc_pd(v_zmask1, _mm256_setzero_pd())) ||
			               (_mm256_testc_pd(v_zmask2, _mm256_setzero_pd()))) {
				*err = 1;
				break;
		}
	}
}


void  is_domain_f32_bt_ones(const float* __restrict in,
	                    const int length,
			    int* err) {

	__m256 v_maskgt(_mm256_setzero_ps());
	__m256 v_masklt(_mm256_setzero_ps());
	const __m256 v_f32one(_mm256_set1_ps(1.F));
	const __m256 v_f32negone(_mm256_set1_ps(-1.F));
	constexpr int f32_off{8};
        *err = 0;

	for (int i{ 0 }; i != length; i += f32_off) {

		v_maskgt = _mm256_cmp_ps(v_f32one, _mm256_loadu_ps(&in[i]),_CMP_GT_OQ);
		v_masklt = _mm256_cmp_ps(v_f32negone, _mm256_loadu_ps(&in[i]), _CMP_LT_OQ);
		if ((_mm256_testc_ps(v_maskgt, _mm256_setzero_ps())) ||
			            (_mm256_testc_ps(v_masklt, _mm256_setzero_ps()))) {
			 *err = 1;
			 break;
		}
	}
}


void  is_domain_f32_bt_ones(const float* __restrict in1,
                            const float* __restrict in2,
			    const int length,
			    int* err             ) {
	
	__m256 v_maskgt1(_mm256_setzero_ps());
	__m256 v_maskgt2(_mm256_setzero_ps());
	__m256 v_masklt1(_mm256_setzero_ps());
	__m256 v_masklt2(_mm256_setzero_ps());
	const __m256 v_f32one(_mm256_set1_ps(1.F));
	const __m256 v_f32negone(_mm256_set1_ps(-1.F));
	constexpr int f32_off{ 8 };
        *err = 0;

	for (int i{ 0 }; i != length; i += f32_off) {

		v_maskgt1 = _mm256_cmp_ps(v_f32one,    _mm256_loadu_ps(&in1[i]), _CMP_GT_OQ);
		v_maskgt2 = _mm256_cmp_ps(v_f32one,    _mm256_loadu_ps(&in2[i]), _CMP_GT_OQ);
		v_masklt1 = _mm256_cmp_ps(v_f32negone, _mm256_loadu_ps(&in1[i]), _CMP_LT_OQ);
		v_masklt2 = _mm256_cmp_ps(v_f32negone, _mm256_loadu_ps(&in2[i]), _CMP_LT_OQ);
		if ((_mm256_testc_ps(v_maskgt1, _mm256_setzero_ps())) ||
			          (_mm256_testc_ps(v_masklt1, _mm256_setzero_ps()))) {
			*err = 1;
			break;
		}

		if ((_mm256_testc_ps(v_maskgt2, _mm256_setzero_ps())) || 
		              (_mm256_testc_ps(v_masklt2, _mm256_setzero_ps()))) {
			*err = 1;
			break;
		}
	}
}


void  is_domain_f64_bt_ones(const double* __restrict in,
			   const int length,
			   int* err) {

	__m256d v_maskgt(_mm256_setzero_pd());
	__m256d v_masklt(_mm256_setzero_pd());
	const __m256d v_fone(_mm256_set1_pd(1.0));
	const __m256d v_fnegone(_mm256_set1_pd(-1.0));
	constexpr int f64_off{4};
        *err = 0;

	for (int i{ 0 }; i != length; i += f64_off) {

		v_maskgt = _mm256_cmp_pd(v_fone,    _mm256_loadu_pd(&in[i]), _CMP_GT_OQ);
		v_masklt = _mm256_cmp_pd(v_fnegone, _mm256_loadu_pd(&in[i]), _CMP_LT_OQ);
		if ((_mm256_testc_pd(v_maskgt, _mm256_setzero_pd())) || 
		              (_mm256_testc_pd(v_masklt, _mm256_setzero_pd()))) {
			 *err = 1;
			 break;
		}
	}
}


void  is_domain_f64_bt_ones(const double* __restrict in1,
                            const double* __restrict in2,
			    const int length,
			    int* err) {

	__m256d v_maskgt1(_mm256_setzero_pd());
	__m256d v_maskgt2(_mm256_setzero_pd());
	__m256d v_masklt1(_mm256_setzero_pd());
	__m256d v_masklt2(_mm256_setzero_pd());
	const __m256d v_fone(_mm256_set1_pd(1.0));
	const __m256d v_fnegone(_mm256_set1_pd(-1.0));
	constexpr int f64_off{4};
        *err = 1

	for (int i{ 0 }; i != length; i += f64_off) {

		v_maskgt1 = _mm256_cmp_pd(v_fone,    _mm256_loadu_pd(&in1[i]), _CMP_GT_OQ);
		v_maskgt2 = _mm256_cmp_pd(v_fone,    _mm256_loadu_pd(&in2[i]), _CMP_GT_OQ);
		v_masklt1 = _mm256_cmp_pd(v_fnegone, _mm256_loadu_pd(&in1[i]), _CMP_LT_OQ);
		v_masklt2 = _mm256_cmp_pd(v_fnegone, _mm256_loadu_pd(&in2[i]), _CMP_LT_OQ);
		if ((_mm256_testc_pd(v_maskgt1, _mm256_setzero_pd())) ||
			(_mm256_testc_pd(v_masklt1, _mm256_setzero_pd()))) {
			*err = 1;
			break;
		}

		if ((_mm256_testc_pd(v_maskgt2, _mm256_setzero_pd())) ||
			(_mm256_testc_pd(v_masklt2, _mm256_setzero_pd()))) {
			*err = 1;
			break;
		}
	}
}


void  is_domain_f32_ne_zero(const float* __restrict in1,
                            const float* __restrict in2,
			    const int  length,
			    int* err) {

	__m256 v_xmaskz(_mm256_setzero_ps());
	__m256 v_ymaskz(_mm256_setzero_ps());
	const __m256 v_f32zero(_mm256_set1_ps(0.F));
	constexpr int f32_off{8};
        *err = 0;

	for (int i{ 0 }; i != length; i += f32_off) {

		v_xmaskz = _mm256_cmp_ps(v_f32zero, _mm256_loadu_ps(&in1[i]), _CMP_EQ_OQ);
		v_ymaskz = _mm256_cmp_ps(v_f32zero, _mm256_loadu_ps(&in2[i]), _CMP_EQ_OQ);
		if ((_mm256_testc_ps(v_xmaskz, _mm256_setzero_ps())) || 
		                    (_mm256_testc_ps(v_ymaskz, _mm256_setzero_ps()))) {

			*err = 1;
			break;
		}
	}
}


void  is_domain_f64_ne_zero(const double* __restrict in1,
                            const double* __restrict in2,
			    const int length,
			    int* err) {

	__m256d v_xmaskz(_mm256_setzero_pd());
	__m256d v_ymaskz(_mm256_setzero_pd());
	const __m256d v_f64zero(_mm256_set1_pd(0.0));
	constexpr int f64_off{4};
        *err = 0;

	for (int i{ 0 }; i != length; i += f64_off) {

		v_xmaskz = _mm256_cmp_pd(v_f64zero, _mm256_loadu_pd(&in1[i]), _CMP_EQ_OQ);
		v_ymaskz = _mm256_cmp_pd(v_f64zero, _mm256_loadu_pd(&in2[i]), _CMP_EQ_OQ);
		if ((_mm256_testc_pd(v_xmaskz, _mm256_setzero_pd())) ||
			                (_mm256_testc_pd(v_ymaskz, _mm256_setzero_pd()))) {
			
			*err = 1;
			break;
		}
	}
}


void  is_domain_f32_ne_zero(const float* __restrict in,
			    const int length,
			    int* err) {

	__m256 v_xmaskz(_mm256_setzero_ps());
	const __m256 v_f32zero(_mm256_set1_ps(0.F));
	constexpr int f32_off{8};
        *err = 0;

	for (int i{ 0 }; i != length; i += f32_off) {

		v_xmaskz = _mm256_cmp_ps(v_f32zero, _mm256_loadu_ps(&in[i]), _CMP_EQ_OQ);
		if ((_mm256_testc_ps(v_xmaskz, _mm256_setzero_ps()))) {

			*err = 1;
			break;
		}
	}
}


void  is_domain_f64_ne_zero(const double* __restrict in,
			    const int length,
			    int* err) {

	__m256d v_xmaskz(_mm256_setzero_pd());
	const __m256d v_f64zero(_mm256_set1_pd(0.0));
	constexpr int f64_off{4};
        *err = 0;

	for (int i{ 0 }; i != length; i += f64_off) {

		v_xmaskz = _mm256_cmp_pd(v_f64zero, _mm256_loadu_pd(&in[i]), _CMP_EQ_OQ);
		if ((_mm256_testc_pd(v_xmaskz, _mm256_setzero_pd()))) {

			*err = 1;
			break;
		}
	}
}


void  is_domain_f32_ne_negone(const float* __restrict in,
			      const int length,
			      int* err) {

	__m256 v_xmaskz(_mm256_setzero_ps());
	const __m256 v_f32negone(_mm256_set1_ps(-1.F));
	constexpr int f32_off{8};
        *err = 0;

	for (int i{ 0 }; i != length; i += f32_off) {

		v_xmaskz = _mm256_cmp_ps(v_f32negone, _mm256_loadu_ps(&in[i]), _CMP_LT_OQ);
		if ((_mm256_testc_ps(v_xmaskz, _mm256_setzero_ps()))) {

			*err = 1;
			break;
		}
	}
}


void  is_domain_f64_ne_negone(const double* __restrict in,
			      const int length,
			      int* err) {

	__m256d v_xmaskz(_mm256_setzero_pd());
	const __m256d v_f64negone(_mm256_set1_pd(-1.0L));
	constexpr int f64_off{4};
        *err = 0;

	for (int i{ 0 }; i != length; i += f64_off) {

		v_xmaskz = _mm256_cmp_pd(v_f64negone, _mm256_loadu_pd(&in[i]), _CMP_LT_OQ);
		if ((_mm256_testc_pd(v_xmaskz, _mm256_setzero_pd()))) {

			*err = 1;
			break;
		}
	}
}


void  is_denormal_f32_present(const float* __restrict in,
			      const int length,
			      int* err) {
	
	using namespace std;
	__m256 v_xmaskz(_mm256_setzero_ps());
	__m256 v_ymaskz(_mm256_setzero_ps());
	const __m256 v_zero(_mm256_set1_ps(0.F));
	const __m256 v_minf32(_mm256_set1_ps(numeric_limits<float>::min()));
	constexpr unsigned int abs_mask{ 0x7FFFFFFF };
	constexpr int f32_off{8};
        *err = 0;

	for (int i{ 0 }; i != length; i += f32_off) {
		
		v_ymaskz = _mm256_cmp_ps(v_zero, _mm256_loadu_ps(&in[i]), _CMP_NEQ_OQ);
		__m256 v_abs = _mm256_and_ps(_mm256_loadu_ps(&in[i]), _mm256_set1_ps(abs_mask));
		v_xmaskz = _mm256_cmp_ps(v_abs,v_minf32, _CMP_LT_OQ);
		if ((_mm256_testc_ps(v_ymaskz, _mm256_setzero_ps())) && 
		                    (_mm256_testc_ps(v_xmaskz, _mm256_setzero_ps()))) {

			*err = 1;
			break;
		}
	}
}


void  is_denormal_f64_present(const double* __restrict in,
			      const int length,
			      int* err) {

	  using namespace std;
	  __m256d v_xmaskz(_mm256_setzero_pd());
	  __m256d v_ymaskz(_mm256_setzero_pd());
	  const __m256d v_zero(_mm256_set1_pd(0.0));
	  const __m256d v_minf64(_mm256_set1_pd(numeric_limits<double>::min()));
	  constexpr unsigned long long abs_mask{ 0x7FFFFFFFFFFFFFFF };
	  constexpr int f64_off{4};
          *err = 0;

	  for (int i{ 0 }; i != length; i += f64_off) {

		  v_ymaskz = _mm256_cmp_pd(v_zero, _mm256_loadu_pd(&in[i]), _CMP_NEQ_OQ);
		  __m256d v_abs = _mm256_and_pd(_mm256_loadu_pd(&in[i]), _mm256_set1_pd(abs_mask));
		  v_xmaskz = _mm256_cmp_pd(v_abs,v_minf64, _CMP_LT_OQ);
		  if ((_mm256_testc_pd(v_ymaskz, _mm256_setzero_pd())) && 
		                      (_mm256_testc_pd(v_xmaskz, _mm256_setzero_pd()))) {
			 
			 *err = 1;
			 break;
		  }
	  }
}


void  check_denormf32_raise_except(const float* __restrict in,
				   const int length,
				   int* err) {

	constexpr int one{1};
	is_denormal_f32_present(in,length,err);
#pragma STD FENV_ACCESS ON
	if (*err == one) {
		if (fetestexcept(FE_DENORMAL) != 0) {
			feclearexcept(FE_DENORMAL); // clear previous exception.
			feraiseexcept(FE_DENORMAL);
		}
		else {
			feraiseexcept(FE_DENORMAL);
		}
	}
}


void  check_denormf64_raise_except(const double* __restrict in,
	                           const int length,
				   int* err) {

	constexpr int one{1};
	is_denormal_f64_present(in,length,err);
#pragma STD FENV_ACCESS ON
	if (*err == one) {
		if (fetestexcept(FE_DENORMAL) != 0) {
			feclearexcept(FE_DENORMAL);
			feraiseexcept(FE_DENORMAL);
		}
		else {
			feraiseexcept(FE_DENORMAL);
		}
	}

}


void  is_domain_f32_nan_inf(const float* __restrict in,
	                   const int length,
			   int* err) {

	using namespace std;

	for (int i{ 0 }; i != length; ++i) {

		if (isinf(in[i]) || isnan(in[i])) {
			
			*err = 1;
			break;
		}
	}
}


void  is_domain_f64_nan_inf(const double* __restrict in,
		            const int length,
			    int* err) {

	using namespace std;

	for (int i{ 0 }; i != length; ++i) {

		if (isinf(in[i]) || isnan(in[i])) {

		    *err = 1;
			break;
		}
	}
}


int   clear_all_feexcept() {

 int err{-9999};
#pragma STD FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_ALL_EXCEPT);
		if(err != 0) {
#if (GMS_DEBUG_ON) == 1
			std::cerr << "****FATAL ERROR****\n"
				      << "'feclearexcept' failed with error value: " << err << "\n"
					  << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";
#endif					  
			return err;
		}

		return err;
	}
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  clear_fe_denormal() {

	int err{-9999};
#pragma STD FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_DENORMAL);
		if(err != 0) {
#if (GMS_DEBUG_ON) == 1
			std::cerr << "*****FATAL-ERROR*****\n"
				<< "'feclearexcept' failed with error value: " << err << "\n"
				<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";
#endif
		    return err;
		}

		return err;
	}
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  clear_fe_inexact() {

	int err{-9999};
#pragma STD FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_INEXACT);
		if(err != 0) {
#if (GMS_DEBUG_ON) == 1
			std::cerr << "*****FATAL-ERROR*****\n"
				<< "'feclearexcept' failed with error value: " << err << "\n"
				<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";
#endif
			return err;
	  }

	  return err;
   }
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  clear_fe_invalid() {

	int err{-9999};
#pragma STD FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_INVALID);
		if(err != 0) {
#if (GMS_DEBUG_ON) == 1
			std::cerr << "*****FATAL-ERROR*****\n"
				<< "'feclearexcept' failed with error value: " << err << "\n"
				<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";
#endif
			  return err;
		}

		return err;
	}
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  clear_fe_divbyzero() {

	int err{-9999};
#pragma STD FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_DIVBYZERO);
		if (err != 0) {
#if (GMS_DEBUG_ON) == 1
			std::cerr << "*****FATAL-ERROR*****\n"
				<< "'feclearexcept' failed with error value: " << err << "\n"
				<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";
#endif
			return err;
		}

		return err;
	}
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  clear_fe_overflow() {

	int err{-9999};
#pragma STD FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_OVERFLOW);
		if (err != 0) {
#if (GMS_DEBUG_ON) == 1
			std::cerr << "*****FATAL-ERROR*****\n"
				<< "'feclearexcept' failed with error value: " << err << "\n"
				<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";
#endif
			return err;
		}

		return err;
	}

#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT'"
#endif

}


int  clear_fe_underflow() {

	int err{-9999};
#pragma STD FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_UNDERFLOW);
		if (err != 0) {
#if (GMS_DEBUG_ON) == 1
			std::cerr << "*****FATAL-ERROR*****\n"
				<< "'feclearexcept' failed with error value: " << err << "\n"
				<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";
#endif
			 return err;
        }

		return err;
  }

#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  test_fe_excepts(_In_ const int excepts) {

	int val{-9999};
#pragma STDC FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {
		
		val = fetestexcept(excepts);
		if (val == FE_ALL_EXCEPT) {
			std::cout << "'fetestexcept' detected, that all fe-exceptions are set!!\n";
			return val;
		}
		else {
			std::cout << "'fetestexcept' detected, that some fe-exceptions are set, but not all!!\n";
			return val;
		}
	}
	return val; // return unmodified val containing -9999 i.e. error.

#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  test_fe_invalid(_In_ const int except) {

	int val{-9999};
#pragma STDC FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {

		val = fetestexcept(except);
		if(val == FE_INVALID) {
			std::cout << "'fetestexcept' detected  FE_INVALID exception!!\n"
				      << "exception numeric type = " << val << ".\n";
			return val;
		}
		else {
			std::cout << "'fetestexcept' failed to detected FE_INVALID exception!!\n"
				      << "exception numeric type = " << val << ".\n";
			return val;
		}
    }

	return val; // return unmodified val containing -9999 i.e. error.
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  test_fe_inexact(_In_ const int except) {

	int val{-9999};
#pragma STDC FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {

		val = fetestexcept(except);
		if (val == FE_INEXACT) {
			std::cout << "'fetestexcept' detected FE_INEXACT exception!!\n"
					  << "exception numeric type = " << val << ".\n";
			return val;
		}
		else {
			std::cout << "'fetestexcept' failed to detected FE_INEXACT exception!!\n"
				      << "exception numeric type = " << val << ".\n";
			return val;
		}
	}
	return val; // return unmodified val containing -9999 i.e. error.
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' '"
#endif

}


int  test_fe_divbyzero(_In_ const int except) {

	int val{-9999};
#pragma STDC FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {

		val = fetestexcept(except);
		if(val == FE_DIVBYZERO) {
			std::cout << "'fetestexcept' detected FE_DIVBYZERO exception!!\n"
				      << "exception numeric type = " << val << ".\n";
			return val;
		}
		else {
			std::cout << "'fetestxcept' failed to detect FE_DIVBYZERO exception!!\n"
				      << "exception numeric type = " << val << ".\n";
			return val;
		}
	}
	return val; // return unmodified val containing -9999 i.e. error.
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  test_fe_unnormal(_In_ const int except) {

	int val{-9999};
#pragma STDC FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {

		val = fetestexcept(except);
		if(val == FE_UNNORMAL) {
			std::cout << "'fetestexcept' detected FE_UNNORMAL exception!!\n"
				      << "exception numeric type = " << val << ".\n";
			return val;
		}
		else {
			std::cout << "'fetestexcept' failed to detected FE_UNNORMAL exception!!\n"
				<< "exception numeric type = " << val << ".\n";
			return val;
		}
	}
	return val; // return unmodified val containing -9999 i.e. error.
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  test_fe_underflow(_In_ const int except) {

	int val{-9999};
#pragma STDC FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {

		val = fetestexcept(except);
		if(val == FE_UNDERFLOW) {
			std::cout << "'fetestexcept' detected FE_UNDERFLOW exception!!\n"
				      << "exception numeric type = " << val << ".\n";
			return val;
		}
		else {
			std::cout << "'fetestexcept' failed to detect FE_UNDERFLOW exception!!\n"
				<< "exception numeric type = " << val << ".\n";
			return val;
		}
    }
	return val; // return unmodified val containing -9999 i.e. error.
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif

}


int  test_fe_overflow(_In_ const int except) {

	int val{-9999};
#pragma STDC FENV_ACCESS ON
#if defined math_errhandling && defined MATH_ERREXCEPT
	if (math_errhandling & MATH_ERREXCEPT) {

		val = fetestexcept(except);
		if (val == FE_OVERFLOW) {
			std::cout << "'fetestexcept' detected FE_OVERFLOW exception!!\n"
				<< "exception numeric type = " << val << ".\n";
			return val;
		}
		else {
			std::cout << "'fetestexcept' failed to detect FE_OVERFLOW exception!!\n"
				<< "exception numeric type = " << val << ".\n";
			return val;
		}
	}
	return val; // return unmodified val containing -9999 i.e. error.
#else
#error "Undefined: 'math_errhandling and MATH_ERREXCEPT' "
#endif
}


