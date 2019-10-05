
#include <iostream>
#include <iomanip>
#include "LAM_perf_collector_qpc.h"
#include "../GMS_malloc.h"
#include "../Math/GMS_constants.h"

#if (USE_MKL) == 1 && defined (GMS_COMPILED_BY_ICC)

		#if !defined (PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP)
		#define PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP(status,data1,data2,desc) \
			do {																  \
					if ((status) != 0) {										  \
					    _mm_free((data1)); _mm_free((data2));				      \
					  DftiFreeDescriptor(&(desc));								  \
					  return (false);											  \
				}																  \
			} while (0);
		 #endif

		 #if !defined (PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP2)
		 #define PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP2(status,data1,data2,desc1,desc2) \
			 do {																		   \
					if ((status) != 0) {												   \
					_mm_free((data1)); _mm_free((data2));								   \
					DftiFreeDescriptor(&(desc1));									       \
					DftiFreeDescriptor(&(desc2));									       \
					return (false);														   \
				}																		   \
			} while (0);
		#endif

#endif

gms::system
::PerfCollectorQPC
::PerfCollectorQPC()
:
m_Iscleared{},
m_nsamples{},
m_func_name{},
m_file_name{},
m_date{},
m_time{},
m_nvalid_values{},
m_start_values{},
m_stop_values{},
m_delta_values{},
m_pc_values{}
{}

gms::system
::PerfCollectorQPC
::PerfCollectorQPC(_In_ const PerfCollectorQPC &x)
:
m_Iscleared{ x.m_Iscleared },
m_nsamples{ x.m_nsamples },
m_func_name{ x.m_func_name },
m_file_name{ x.m_file_name },
m_date{ x.m_date },
m_time{ x.m_time },
m_nvalid_values{ x.m_nvalid_values },
m_start_values{ x.m_start_values },
m_stop_values{ x.m_stop_values },
m_delta_values{ x.m_delta_values },
m_pc_values{ x.m_pc_values }
{}

gms::system
::PerfCollectorQPC
::PerfCollectorQPC(_In_ PerfCollectorQPC &&x)
:
m_Iscleared{ std::move(x.m_Iscleared) },
m_nsamples{ std::move(x.m_nsamples) },
m_func_name{ std::move(x.m_func_name) },
m_file_name{ std::move(x.m_file_name) },
m_date{ std::move(x.m_date) },
m_time{ std::move(x.m_time) },
m_nvalid_values{ std::move(x.m_nvalid_values) },
m_start_values{ std::move(x.m_start_values) },
m_stop_values{ std::move(x.m_stop_values) },
m_delta_values{ std::move(x.m_delta_values) },
m_pc_values{ std::move(x.m_pc_values) }
{}

gms::system::PerfCollectorQPC &
gms::system::PerfCollectorQPC
::operator=(_In_ const PerfCollectorQPC &x) {
	if (this == &x) return (*this);
			m_Iscleared = x.m_Iscleared;
			m_nsamples  = x.m_nsamples;
			m_func_name.operator=(x.m_func_name);
			m_file_name.operator=(x.m_file_name);
			m_date.operator=(x.m_date);
			m_time.operator=(x.m_time);
			m_nvalid_values.operator=(x.m_nvalid_values);
			m_start_values.operator=(x.m_start_values);
			m_stop_values.operator=(x.m_stop_values);
			m_delta_values.operator=(x.m_delta_values);
			m_pc_values.operator=(x.m_pc_values);
	 return (*this);
}

gms::system::PerfCollectorQPC &
gms::system::PerfCollectorQPC
::operator=(_In_ PerfCollectorQPC &&x) {
	if (this == &x) return (*this);
			m_Iscleared = std::move(x.m_Iscleared);
			m_nsamples  = std::move(x.m_nsamples);
			m_func_name.operator=(std::move(x.m_func_name));
			m_file_name.operator=(std::move(x.m_file_name));
			m_date.operator=(std::move(x.m_date));
			m_time.operator=(std::move(x.m_time));
			m_nvalid_values.operator=(std::move(x.m_nvalid_values));
			m_start_values.operator=(std::move(x.m_start_values));
			m_stop_values.operator=(std::move(x.m_stop_values));
			m_delta_values.operator=(std::move(x.m_delta_values));
			m_pc_values.operator=(std::move(x.m_pc_values));
	 return (*this);
}

bool
gms::system
::PerfCollectorQPC::start() {
	bool bOk{};
	LARGE_INTEGER start = {0ULL};
	bOk = ::QueryPerformanceCounter(&start);
	if (!bOk) { // On error do not push_back an invalid value.
		std::printf("QueryPerformanceCounter -- failed with an error: 0x%x\n", ::GetLastError());
		return (false);
	}
	m_start_values.push_back(start); // May throw here
	return (true);
}

bool
gms::system
::PerfCollectorQPC::stop() {
	
	bool bOk1{}, bOk2{};
	LARGE_INTEGER stop = {0Ui64}, cfreq = {0Ui64};
	bOk1 = ::QueryPerformanceCounter(&stop);
	if (!bOk1) {
		std::printf("QueryPerformanceCounter -- failed with an error: 0x%x\n", ::GetLastError());
		return (false);
	}
	
	bOk2 = ::QueryPerformanceFrequency(&cfreq);
	if (!bOk2) {
		std::printf("QueryPerformanceFrequency -- failed with an error: 0x%x\n", ::GetLastError());
		return (false);
	}
	m_stop_values.push_back(stop);
	m_pc_values.push_back(cfreq);
	return (true);
}

void
gms::system
::PerfCollectorQPC::clear_all() {
	using namespace gms::math::constants;
	if (!m_Iscleared) {	
		if (m_nvalid_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_nvalid_values.size(); ++i) {
				m_nvalid_values.operator[](i) = false;
			}
		}
		if (m_start_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
				m_start_values.operator[](i).QuadPart = 0Ui64;
			}
		}
		if (m_stop_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_stop_values.size(); ++i) {
				m_stop_values.operator[](i).QuadPart = 0Ui64;
			}
		}
		if (m_delta_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
				m_delta_values.operator[](i) = dinf;
			}
		}
		if (m_pc_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_pc_values.size(); ++i) {
				m_pc_values.operator[](i).QuadPart = 0Ui64;
			}
		}
		m_Iscleared = true;
	}
}

bool
gms::system
::PerfCollectorQPC
::compute_delta() {
	if (!check_data_size()              || 
	    m_start_values.size() <= lowest ||
		m_stop_values.size()  <= lowest ||
		m_pc_values.size()    <= lowest) {
		return (false);
	}
	ULONGLONG res{};
	double tmp{},dfreq{};
	constexpr double mhz = 1000.0;
	for (std::size_t i = 0Ui64; i != m_stop_values.size(); ++i) {
		if ((res = m_stop_values[i].QuadPart - m_start_values[i].QuadPart) > 0Ui64 ) {
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
			tmp = __binary64_from_uint64(m_pc_values[i].QuadPart);
			dfreq = __binary64_div_binary64_binary64(tmp,mhz);
			m_delta_values.push_back( __binary64_div_binary64_binary64(
									  __binary64_from_uint64(res),dfreq));
			m_nvalid_values.push_back(true);
			++m_nsamples;
#else
			tmp = static_cast<double>(m_pc_values[i].QuadPart);
			dfreq = tmp / mhz;
			m_delta_values.push_back(static_cast<double>(res) / dfreq);
			m_nvalid_values.push_back(true);
			++m_nsamples;
#endif			  
		}
		else {
			m_delta_values.push_back(0.0);
			m_nvalid_values.push_back(false);
		}
	}
	return (true);
}

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
bool
gms::system
::PerfCollectorQPC
::compute_stats(_Out_ double &mean,
				_Out_ double &adev,
				_Out_ double &sdev,
				_Out_ double &skew,
				_Out_ double &kurt) {
	using namespace gms::math::constants;
	if (m_Iscleared || m_delta_values.size() < lo_bound) {
		    mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	std::size_t vlen{};
	__declspec(align(64)) struct {
		double len{}, sum{}, var{}, t{}, t2{}, t3{},
		t4{}, t5{}, isdev{};
	}r8_loc;
	vlen = m_delta_values.size();
	r8_loc.len = __binary64_from_uint64(vlen);
	for(std::size_t i = 0Ui64; i != vlen; ++i) {
		r8_loc.sum = 
			__binary64_add_binary64_binary64(r8_loc.sum,m_delta_values[i]);
	}
	mean = __binary64_div_binary64_binary64(r8_loc.sum,r8_loc.len);
	//  Compute average deviation and variance
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (__binary64_quiet_not_equal_binary64(0.0, m_delta_values[i])) {
			r8_loc.t = __binary64_abs(
							__binary64_sub_binary64_binary64(m_delta_values[i],mean));
			adev = __binary64_add_binary64_binary64(adev,r8_loc.t);
			r8_loc.t2 = __binary64_mul_binary64_binary64(
				             __binary64_sub_binary64_binary64(m_delta_values[i],mean),
				             __binary64_sub_binary64_binary64(m_delta_values[i],mean));
			r8_loc.var = __binary64_add_binary64_binary64(r8_loc.var,r8_loc.t2);
		}  
	}
	adev = __binary64_div_binary64_binary64(adev,r8_loc.len);
	if (__binary64_quiet_less_binary64(r8_loc.var, 0.0)) {
		std::printf(" PerfCollectorQPC::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		    mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	r8_loc.var = __binary64_div_binary64_binary64(r8_loc.var,r8_loc.len);
	sdev = __binary64_sqrt_binary64(r8_loc.var);
	r8_loc.isdev = __binary64_div_binary64_binary64(1.0,sdev);
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (__binary64_quiet_not_equal_binary64(0.0, m_delta_values[i])) {
			r8_loc.t3 = __binary64_mul_binary64_binary64(
				              __binary64_sub_binary64_binary64(m_delta_values[i],mean),r8_loc.isdev);
			r8_loc.t4 = __binary64_mul_binary64_binary64(r8_loc.t3,  
							  __binary64_mul_binary64_binary64(r8_loc.t3,r8_loc.t3));
			skew = __binary64_add_binary64_binary64(skew,r8_loc.t4);
			r8_loc.t5 = __binary64_mul_binary64_binary64(r8_loc.t4,r8_loc.t3);
			kurt = __binary64_add_binary64_binary64(kurt,r8_loc.t5);
		}
	}
	skew = __binary64_div_binary64_binary64(skew,r8_loc.len);
	kurt = __binary64_div_binary64_binary64(kurt,
		__binary64_sub_binary64_binary64(r8_loc.len, 3.0));
	if (__binary64_quiet_less_equal_binary64(kurt, 1.0)) {
		std::printf(" PerfCollectorQPC::compute_stats: Invalid kurtosis: %.15f\n", kurt);
		mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	return (true);
}

#else

bool
gms::system
::PerfCollectorQPC::
compute_stats(_Inout_ double &mean,
			  _Inout_ double &adev,
			  _Inout_ double &sdev,
			  _Inout_ double &skew,
			  _Inout_ double &kurt) {
	using namespace gms::math::constants;
	if (m_Iscleared || m_delta_values.size() < lo_bound) {
		    mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	uint64_t vlen{};
	__declspec(align(64)) struct{
		double len{}, sum{}, var{}, t{}, t2{}, t3{}, tmp{},
		t4{}, t5{}, isdev{}, fracp{}, ct2{};
	}r8_loc;
	vlen = m_delta_values.size();
	r8_loc.len = static_cast<double>(vlen);
	for(std::size_t i = 0Ui64; i != vlen; ++i) {
		sum += m_delta_values[i];
	}
	mean = sum / r8_loc.len;
	//  Compute average deviation and variance
	for (size_t i = 0Ui64; i != vlen; ++i) { // <-- no auto-vectorization (no countable loop)
		if (0.0 != m_delta_values.operator[](i)) {
			
			r8_loc.t = std::abs(m_delta_values[i] - mean); //potential catastrophic cancellation if(tmp - mean) both are very close to each other.
			adev += r8_loc.t; // <-- here potential overflow
			r8_loc.t2 = (m_delta_values[i] - mean) * (m_delta_values[i] - mean);
			r8_loc.ct2 = r8_loc.t2; // Checks if catastrophic cancellation has occurred
			r8_loc.fracp = r8_loc.t2 - static_cast<uint64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" PerfCollectorQPC::compute_stats: Losing a significand digits: %.16f\n", r8_loc.fracp);
			}
			r8_loc.var += r8_loc.t2; // potential overflow
		}
	}
	adev /= r8_loc.len;
	if (r8_loc.var < 0.0) {
		std::printf(" PerfCollectorQPC::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		return (false);
	}
	r8_loc.var /= r8_loc.len;
	sdev = std::sqrt(r8_loc.var);
	r8_loc.isdev = 1.0 / sdev;
	r8_loc.fracp = -1.0;
	for (size_t i = 0Ui64; i != vlen; ++i) {
		if (0.0 != m_delta_values.operator[](i)) {
			
			r8_loc.t3 = (m_delta_values[i] - mean) * r8_loc.isdev; // <-- catastrophic cancellation here
			r8_loc.ct2 = r8_loc.t3;
			r8_loc.fracp = r8_loc.t3 - static_cast<int64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" PerfCollectorQPC::compute_stats: Losing a significant digits: %.16f\n", r8_loc.fracp);
			}

			r8_loc.t4 = r8_loc.t3*r8_loc.t3*r8_loc.t3; // <-- // Potential overflow?
			skew += r8_loc.t4; // <-- // Potential overflow?
			r8_loc.t5 = r8_loc.t4 * r8_loc.t3;
			kurt += r8_loc.t5; // <-- // Potential overflow?
		}
	}
	skew /= r8_loc.len;
	kurt /= r8_loc.len - 3.0;
	if (kurt < 1.0) {
		std::printf(" PerfCollectorQPC::compute_stats: Invalid kurtosis: %.15f\n", kurt);
		return (false);
	}
	return (true);
}

#endif

bool
gms::system
::PerfCollectorQPC
::correlation_set(_In_ PerfCollectorQPC &set1,
				  _In_ PerfCollectorQPC &set2,
				  _Out_ double * __restrict correlated_set,
				  _In_ int32_t * __restrict ip,
				  _In_ double * __restrict w,
				  _In_ const int32_t n) {
	if (set1.m_delta_values.size() != set2.m_delta_values.size() ||
		static_cast<size_t>(n) != set1.m_delta_values.size()) {
		return (false);
	}
#if (USE_MKL) == 1
	DFTI_DESCRIPTOR_HANDLE ddh1, ddh2;
	MKL_LONG status{};
#endif
	double t{}, invals{};
	int32_t snvals{};
	auto set_size = set1.m_delta_values.size();
	typedef double * __restrict __declspec(align_value(64)) aligned_r8ptr;
	aligned_r8ptr dataset1 = gms::common::gms_edmalloca(set_size,align64B);
	aligned_r8ptr dataset2 = gms::common::gms_edmalloca(set_size,align64B);
	snvals = n << 1;
	invals = 1.0 / static_cast<double>(snvals);
	// MKL_DFTI
#if (USE_MKL) == 1
	status = DftiCreateDescriptor(&ddh1, DFTI_DOUBLE, DFTI_REAL, 1, set_size);
	PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP(status, dataset1, dataset2, ddh1)
	status = DftiCommitDescriptor(ddh1);
	PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP(status, dataset1, dataset2, ddh1)
	status = DftiComputeForward(ddh1, dataset1);
	PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP(status, dataset1, dataset2, ddh1)
	status = DftiCreateDescriptor(&ddh2, DFTI_DOUBLE, DFTI_REAL, 1, set_size);
	PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	status = DftiCommitDescriptor(ddh2);
	PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	status = DftiComputeForward(ddh2, dataset2);
	PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	for (int32_t i = 0; i != n - 1; i += 2) {
		t = dataset1[i];
		dataset1[i] = (dataset1[i] * dataset2[i] + dataset1[i + 1] * dataset2[i + 1]) * invals;
		dataset1[i + 1] = (dataset1[i + 1] * dataset2[i] - t * dataset2[i + 1]) * invals;
	}
	dataset1[0] = dataset1[0] * dataset2[0] * invals;
	dataset1[1] = dataset1[1] * dataset2[1] * invals;
	status = DftiComputeBackward(ddh1, dataset1);
	PERF_COLLECTOR_QPC_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	memcpy(&correlated_set[0], &dataset1[0], set_size);
	_mm_free(dataset1); _mm_free(dataset2);
	DftiFreeDescriptor(&ddh1); DftiFreeDescriptor(&ddh2);
	return (true);
#else
	ip[0] = 0;
	rdft(n, 1, &dataset1[0], &ip[0], &w[0]);
	rdft(n, 1, &dataset2[0], &ip[0], &w[0]);
	for (int32_t i = 1; i != n - 1; i += 2) {
		t = dataset1[i];
		dataset1[i] = (dataset1[i] * dataset2[i] + dataset1[i + 1] * dataset2[i + 1]) * invals;
		dataset1[i + 1] = (dataset1[i + 1] * dataset2[i] - t * dataset2[i + 1]) * invals;
	}
	dataset1[0] = dataset1[0] * dataset2[0] * invals;
	dataset1[1] = dataset1[1] * dataset2[1] * invals;
	rdft(n, -1, &dataset1[0], &ip[0], &w[0]);
	memcpy(&correlated_set[0], &dataset1[0], set_size);

	if (dataset1) _mm_free(dataset1);
	if (dataset2) _mm_free(dataset2);
	return (true);
#endif
}

void
gms::system
::PerfCollectorQPC::print() const {
	std::cout << " Dumping state of: " << typeid(*this).raw_name() << "\n";
	std::cout << "======================================================================================\n";
	std::cout << " Reset state of collector: " << std::boolalpha << m_Iscleared << "\n"
		<< "       Number of valid samples:  " << m_nsamples << "\n"
		<< "       Collected at function:    " << m_func_name.data() << "\n"
		<< "       Collected in file:        " << m_file_name.data() << "\n"
		<< "       Collected at date:        " << m_date.data() << "\n"
		<< "       Collected at time:        " << m_time.data() << "\n"
		<< "       ======================================================================================\n";
	std::cout << "  Sample state  :  start_value   :    stop_value   :    delta    :    counter          \n";
	for (std::size_t i = 0Ui64; i != m_stop_values.size(); ++i) {
		std::cout << std::setw(4)  << m_nvalid_values[i] <<
			         std::setw(8)  << m_start_values[i].QuadPart   <<
					 std::setw(12) << m_start_values[i].QuadPart   <<
					 std::setw(16) << std::setprecision(15) << m_delta_values[i] <<
					 std::setw(20) << m_pc_values[i].QuadPart << "\n";
	}
	std::cout << "======================================================================================\n"
		      << "===============================	End  of dump  ======================================\n";
}


// Helpers
bool
gms::system
::PerfCollectorQPC
::check_data_size() const {
	if (!m_Iscleared) {
		if (m_stop_values.size() == m_start_values.size() &&
			m_pc_values.size() == m_stop_values.size()) {
			return (true);
		}
		else {
			return (false);
		}
	}
	else return (false);
}

