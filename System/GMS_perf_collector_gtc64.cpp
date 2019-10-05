
#include <iostream>
#include <iomanip>
#include "GMS_perf_collector_gtc64.h"
#include "../GMS_malloc.h"

#include "../Math/GMS_constants.h"

#if (USE_MKL) == 1 && defined (GMS_COMPILED_BY_ICC)

		#if !defined (PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP)
		#define PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP(status,data1,data2,desc)  \
	       do {														             \
		        if ((status) != 0) {									         \
		            _mm_free(data1), _mm_free(data2), DftiFreeDescriptor(&desc); \
		             return (false);											 \
	             }																 \
	      } while (0);
      #endif

	  #if !defined (PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP2)
	  #define PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP2(status,data1,data2,desc1,desc2)  \
		do {																		   \
				if ((status) != 0) {												   \
				_mm_free(data1); _mm_free(data2);									   \
				DftiFreeDescriptor(&desc1);											   \
				DftiFreeDescriptor(&desc2);											   \
				return (false);														   \
			}																		   \
		} while (0);
	  #endif

#endif

gms::system
::PerfCollectorGTC64
::PerfCollectorGTC64()
:
m_Iscleared{},
m_nsamples{},
m_funcname{},
m_filename{},
m_date{},
m_time{},
m_nsuccess(),
m_start_values(),
m_stop_values(),
m_delta_values()
{}

gms::system
::PerfCollectorGTC64
::PerfCollectorGTC64(_In_ const PerfCollectorGTC64 &x)
:
m_Iscleared{ x.m_Iscleared },
m_nsamples{ x.m_nsamples },
m_funcname{ x.m_funcname },
m_filename{ x.m_filename },
m_date{ x.m_date },
m_time{x.m_time},
m_nsuccess(x.m_nsuccess),
m_start_values(x.m_start_values),
m_stop_values(x.m_stop_values),
m_delta_values(x.m_delta_values)
{}

gms::system
::PerfCollectorGTC64
::PerfCollectorGTC64(_In_ PerfCollectorGTC64 &&x)
:
m_Iscleared{ x.m_Iscleared },
m_nsamples{ x.m_nsamples },
m_funcname{ std::move(x.m_funcname )},
m_filename{ std::move(x.m_filename) },
m_date{ std::move(x.m_date) },
m_time{ std::move(x.m_time) },
m_nsuccess(std::move(x.m_nsuccess)),
m_start_values(std::move(x.m_start_values)),
m_stop_values(std::move(x.m_stop_values)),
m_delta_values(std::move(x.m_delta_values))
{}

gms::system::PerfCollectorGTC64 &
gms::system::PerfCollectorGTC64
::operator=(_In_ const PerfCollectorGTC64 &x) {
	if (this == &x) return (*this);
		m_Iscleared = x.m_Iscleared,
		m_nsamples = x.m_nsamples,
		m_funcname.operator=(x.m_funcname),
		m_filename.operator=(x.m_filename),
		m_date.operator=(x.m_date),
		m_time.operator=(x.m_time),
		m_nsuccess.operator=(x.m_nsuccess),
		m_start_values.operator=(x.m_start_values),
		m_stop_values.operator=(x.m_stop_values),
		m_delta_values.operator=(x.m_delta_values);
	 return (*this);
}

gms::system::PerfCollectorGTC64 &
gms::system::PerfCollectorGTC64
::operator=(_In_ PerfCollectorGTC64 &&x) {
	if (this == &x) return (*this);
	    m_Iscleared = std::move(x.m_Iscleared),
		m_nsamples = std::move(x.m_nsamples),
		m_funcname.operator=(std::move(x.m_funcname)),
		m_filename.operator=(std::move(x.m_filename)),
		m_date.operator=(std::move(x.m_date)),
		m_time.operator=(std::move(x.m_time)),
		m_nsuccess.operator=(std::move(x.m_nsuccess)),
		m_start_values.operator=(std::move(x.m_start_values)),
		m_stop_values.operator=(std::move(x.m_stop_values)),
		m_delta_values.operator=(std::move(x.m_delta_values));
	return (*this);
}

void
gms::system
::PerfCollectorGTC64::start() {
	m_start_values.push_back(::GetTickCount64());
}

void
gms::system
::PerfCollectorGTC64::stop() {
	m_stop_values.push_back(::GetTickCount64());
}

void
gms::system
::PerfCollectorGTC64::clear_all(){
	if (!m_Iscleared) {
	    
		if (m_start_values.size() >= lowest) {
			for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
				m_start_values.operator[](i) = zero;
			}
		}
		if (m_stop_values.size() >= lowest) {
			for (std::size_t i = 0Ui64; i != m_stop_values.size(); ++i) {
				m_stop_values.operator[](i) = zero;
			}
		}
		if (m_delta_values.size() >= lowest) {
			for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
				m_delta_values.operator[](i) = zero;
			}
		}
		m_Iscleared = true;
	}
}

bool
gms::system
::PerfCollectorGTC64
::compute_delta() {
	if (check_data_size() == false){
	    return (false);
	}
	uint64_t res{};
	for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
		if ((res = m_stop_values[i] - m_start_values[i]) > zero) {
			m_delta_values.push_back(res);
			m_nsuccess.push_back(true);
			++m_nsamples;
		}
		else {
			m_delta_values.push_back(zero);
			m_nsuccess.push_back(false);
		}
	}
	return (true);
}

#if (USE_ACCURATE_IEEE754_2008_FP) == 1

bool 
gms::system::PerfCollectorGTC64
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
	uint64_t s{}, prev{}, vlen{};
	__declspec(align(64)) struct {
		double len{},  var{}, t{}, t2{}, t3{}, tmp{},
		t4{}, t5{}, isdev{}; 
	}r8_loc;
	vlen = m_delta_values.size();
	r8_loc.len = __binary64_from_uint64(vlen);
	// Compute mean guard against an overflow
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (s < prev) return (false);
		prev = s;
		s += m_delta_values.operator[](i);
	}
	mean = __binary64_from_uint64(s) / r8_loc.len;
	//  Compute average deviation and variance
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (zero != m_delta_values[i]) {
			r8_loc.tmp = __binary64_from_uint64(m_delta_values[i]);
			r8_loc.t = __binary64_abs(__binary64_sub_binary64_binary64(r8_loc.tmp,mean));
			adev = __binary64_add_binary64_binary64(adev,r8_loc.t);
			r8_loc.t2 = __binary64_mul_binary64_binary64(
							     __binary64_sub_binary64_binary64(r8_loc.tmp,mean),
								 __binary64_sub_binary64_binary64(r8_loc.tmp,mean));
			r8_loc.var = __binary64_add_binary64_binary64(r8_loc.var,r8_loc.t2);
		}
	}
	adev = __binary64_div_binary64_binary64(adev,r8_loc.len);
	if (__binary64_quiet_less_binary64(r8_loc.var, 0.0)) {
		std::printf(" PerfCollectorGTC64::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		mean = dinf, adev = dinf,
		sdev = dinf, skew = dinf,
		kurt = dinf;
		return (false);
	}
	r8_loc.var = __binary64_div_binary64_binary64(r8_loc.var,r8_loc.len);
	sdev = __binary64_sqrt_binary64(r8_loc.var);
	r8_loc.isdev = __binary64_div_binary64_binary64(1.0,sdev);
	r8_loc.tmp = 0.0;
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (zero != m_delta_values[i]) {
			r8_loc.tmp = __binary64_from_uint64(m_delta_values[i]);
			r8_loc.t3 =  __binary64_mul_binary64_binary64(
			                   __binary64_sub_binary64_binary64(r8_loc.tmp,mean),r8_loc.isdev);
			r8_loc.t4 = __binary64_mul_binary64_binary64(r8_loc.t3, 
				                        __binary64_mul_binary64_binary64(r8_loc.t3,r8_loc.t3));
			skew = __binary64_add_binary64_binary64(skew,r8_loc.t4);
			r8_loc.t5 = __binary64_mul_binary64_binary64(r8_loc.t4,r8_loc.t3);
			kurt = __binary64_add_binary64_binary64(kurt,r8_loc.t5);
		}
	}
	skew = __binary64_div_binary64_binary64(skew,r8_loc.len);
	kurt = __binary64_div_binary64_binary64(kurt,
				      __binary64_sub_binary64_binary64(r8_loc.len,3.0));
	if (__binary64_quiet_less_equal_binary64(kurt, 1.0)) {
		std::printf(" PerfCollectorGTC64::compute_stats: Invalid kurtosis: %.15f\n", kurt);
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
::PerfCollectorGTC64
::compute_stats(_Inout_ double &mean,
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
	uint64_t s{}, prev{};
	__declspec(align(64)) struct{
		double len{}, sum{}, var{}, t{}, t2{}, t3{}, tmp{},
		t4{}, t5{}, isdev{}, fracp{}, ct2{};
	}r8_loc;
	r8_loc.len = static_cast<double>(m_delta_values.size());
	// Compute mean (guard against an overflow)
	for (size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		if (s < prev) return (false);
		prev = s;
		s += m_delta_values.operator[](i);
	}
	mean = static_cast<double>(s) / r8_loc.len;
	//  Compute average deviation and variance
	for (size_t i = 0Ui64; i != m_delta_values.size(); ++i) { // <-- no auto-vectorization (no countable loop)
		if (zero != m_delta_values.operator[](i)) {
			r8_loc.tmp = static_cast<double>(m_delta_values[i]);
			r8_loc.t   = std::abs(r8_loc.tmp - mean); //potential catastrophic cancellation if(tmp - mean) both are very close to each other.
			adev += r8_loc.t; // <-- here potential overflow
			r8_loc.t2 = (r8_loc.tmp - mean) * (r8_loc.tmp - mean);
			r8_loc.ct2 = r8_loc.t2; // Checks if catastrophic cancellation has occurred
			r8_loc.fracp = r8_loc.t2 - static_cast<uint64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" PerfCollectorGTC64::compute_stats: Losing a significand digits: %.16f\n", r8_loc.fracp);
			}
			r8_loc.var += r8_loc.t2; // potential overflow
		}
	}
	adev /= r8_loc.len;
	if (r8_loc.var < 0.0) {
		std::printf(" PerfCollectorGTC64::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		return (false);
	}
	r8_loc.var /= r8_loc.len;
	sdev = std::sqrt(r8_loc.var);
	r8_loc.isdev = 1.0 / sdev;
	r8_loc.fracp = -1.0;
	r8_loc.tmp = 0.0;
	for (size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		if (m_delta_values.operator[](i) != zero) {
			r8_loc.tmp = static_cast<double>(m_delta_values[i]);
			r8_loc.t3 = (r8_loc.tmp - mean) * r8_loc.isdev; // <-- catastrophic cancellation here
			r8_loc.ct2 = r8_loc.t3;
			r8_loc.fracp = r8_loc.t3 - static_cast<int64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" PerfCollectorGTC64::compute_stats: Losing a significant digits: %.16f\n", r8_loc.fracp);
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
		std::printf(" PerfCollectorGTC64::compute_stats: Invalid kurtosis: %.15f\n", kurt);
		return (false);
	}
	return (true);
}

#endif

bool
gms::system::PerfCollectorGTC64
::correlation_set(_In_ PerfCollectorGTC64 &set1,
				  _In_ PerfCollectorGTC64 &set2,
				  _Out_ double * __restrict correlated_set,
				  _In_ int32_t * __restrict ip,
				  _In_ double * __restrict w,
				  _In_ const int32_t n) {
	if (set1.m_delta_values.size() != set2.m_delta_values.size() ||
		static_cast<size_t>(n) != set1.m_delta_values.size()) {
		return (false);
	}
#if (USE_MKL) == 1
	DFTI_DESCRIPTOR_HANDLE ddh1,ddh2;
	MKL_LONG status{};
#endif
	double t{}, invals{};
	int32_t snvals{};
	auto set_size = set1.m_delta_values.size();
	typedef double * __restrict __declspec(align_value(64)) aligned_r8ptr;
	aligned_r8ptr dataset1 = gms::common::gms_edmalloca(set_size,align64B);
	aligned_r8ptr dataset2 = gms::common::gms_edmalloca(set_size,align64B);
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
	 set1.cvrt_to_double(&dataset1[0],set_size);
	 set2.cvrt_to_double(&dataset2[0],set_size);
#else
	 set1.cvrt_to_double(&dataset1[0],set_size);
	set2.cvrt_to_double(&dataset2[0], set_size);
#endif
	
	snvals = n << 1;
	invals = 1.0 / static_cast<double>(snvals);
	// MKL DFTI
#if (USE_MKL) == 1
	status = DftiCreateDescriptor(&ddh1,DFTI_DOUBLE,DFTI_REAL,1,set_size);
	PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiCommitDescriptor(ddh1);
	PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiComputeForward(ddh1,dataset1);
	PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiCreateDescriptor(&ddh2,DFTI_DOUBLE,DFTI_REAL,1,set_size);
	PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP2(status,dataset1,dataset2,ddh1,ddh2)
	status = DftiCommitDescriptor(ddh2);
	PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2,ddh1, ddh2)
	status = DftiComputeForward(ddh2,dataset2);
	PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP2(status,dataset1,dataset2,ddh1,ddh2)
	for (int32_t i = 0; i != n - 1; i += 2) {
		t = dataset1[i];
		dataset1[i] = (dataset1[i] * dataset2[i] + dataset1[i + 1] * dataset2[i + 1]) * invals;
		dataset1[i + 1] = (dataset1[i + 1] * dataset2[i] - t * dataset2[i + 1]) * invals;
	}
	dataset1[0] = dataset1[0] * dataset2[0] * invals;
	dataset1[1] = dataset1[1] * dataset2[1] * invals;
	status = DftiComputeBackward(ddh1,dataset1);
	PERF_COLLECTOR_GTC64_DFTI_FAIL_CLEANUP2(status,dataset1,dataset2,ddh1,ddh2)
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
::PerfCollectorGTC64::print() const {
	std::cout << " Dumping state of: " << typeid(*this).raw_name() << "\n";
	std::cout << "======================================================================================\n";
	std::cout << " Reset state of collector: " << std::boolalpha << m_Iscleared << "\n"
		<< "       Number of valid samples:  " << m_nsamples << "\n"
		<< "       Collected at function:    " << m_funcname.data() << "\n"
		<< "       Collected in file:        " << m_filename.data() << "\n"
		<< "       Collected at date:        " << m_date.data() << "\n"
		<< "       Collected at time:        " << m_time.data() << "\n"
		<< "       ======================================================================================\n";
	std::cout << "  Sample_state  :    start_value  :   stop_value    :  delta_value    \n";
	for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
		std::cout << std::setw(4)  << m_nsuccess[i]          <<
			         std::setw(8)  << m_start_values[i]      <<
					 std::setw(12) << m_stop_values[i]       <<
					 std::setw(16) << m_delta_values[i] << "\n";
 	}				 
	std::cout << "======================================================================================\n"
		      << "===============================	End  of dump  ======================================\n";		
}

// Helpers
bool
gms::system
::PerfCollectorGTC64
::check_data_size() {
	if (m_Iscleared == false) {
		if (m_start_values.size() == m_stop_values.size() &&
			m_start_values.size() == m_delta_values.size()) {
			return (true);
		}
		else {
			return (false);
		}
	}
	else return (false);
}

void
gms::system
::PerfCollectorGTC64
::cvrt_to_double(_Inout_ double * __restrict data,
			    _In_ const std::size_t data_len) {
	
	for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		data[i] = static_cast<double>(m_delta_values[i]);
	}
	
}
