
#include <iostream>
#include <iomanip>
#include <cstdint>
#include "LAM_gtc64_timer.h"
#include "../GMS_indices.h"
#include "../Math/GMS_constants.h"


gms::system
::GTC64Timer
::GTC64Timer(_In_ const std::size_t nruns,
			 _In_ const std::size_t datum_size,
			 _In_ const fptr tfunc,
			 _In_ const fptr wfunc,
			 _In_ const char * pfuname,
			 _In_ const char * pfilename)
:
m_nruns{ nruns },
m_datum_size{ datum_size },
m_ptfunc( tfunc ),
m_ptwrap( wfunc),
m_func_name{ pfuname },
m_file_name{pfilename},
m_timing_values(0ULL, m_nruns * m_datum_size),
m_overhead_values(0ULL, m_nruns * m_datum_size),
m_delta_values(0ULL, m_nruns * m_datum_size)
{}


gms::system
::GTC64Timer
::GTC64Timer(_In_ const GTC64Timer &x)
:
m_nruns{ x.m_nruns },
m_datum_size{ x.m_datum_size },
m_ptfunc( x.m_ptfunc ),
m_ptwrap(x.m_ptwrap),
m_func_name{ x.m_func_name },
m_file_name{ x.m_file_name },
m_timing_values{ x.m_timing_values },
m_overhead_values{ x.m_overhead_values },
m_delta_values{ x.m_delta_values }
{}


gms::system
::GTC64Timer
::GTC64Timer(_In_ GTC64Timer &&x)
:
m_nruns{ x.m_nruns },
m_datum_size{ x.m_datum_size },
m_ptfunc( x.m_ptfunc ),
m_ptwrap( x.m_ptwrap ),
m_func_name{ std::move(x.m_func_name) },
m_file_name{ std::move(x.m_file_name) },
m_timing_values{ std::move(x.m_timing_values) },
m_overhead_values{ std::move(x.m_overhead_values) },
m_delta_values{ std::move(x.m_delta_values) }
{}


gms::system::GTC64Timer
::~GTC64Timer() {
	if (m_ptfunc != nullptr) m_ptfunc = nullptr;
	if (m_ptwrap != nullptr) m_ptwrap = nullptr;
}

gms::system::GTC64Timer &
gms::system::GTC64Timer
::operator=(_In_ const GTC64Timer &x) {
	if (this == &x) return (*this);

	    m_nruns      = x.m_nruns;
	    m_datum_size = x.m_datum_size;
	    m_ptfunc     =  x.m_ptfunc;
		m_ptwrap     = x.m_ptwrap;
	    m_func_name.operator=(x.m_func_name);
	    m_file_name.operator=(x.m_file_name);
		m_timing_values.operator=(x.m_timing_values);
		m_overhead_values.operator=(x.m_overhead_values);
		m_delta_values.operator=(x.m_delta_values);
	   
	return (*this);
}

gms::system::GTC64Timer &
gms::system::GTC64Timer
::operator=(_In_ GTC64Timer &&x) {
	if (this == &x) return (*this);

		m_nruns      = x.m_nruns;
		m_datum_size = x.m_datum_size;
		m_ptfunc     = x.m_ptfunc;
		m_ptwrap     = x.m_ptwrap;
		m_func_name.operator=(std::move(x.m_func_name));
		m_file_name.operator=(std::move(x.m_file_name));
		m_timing_values.operator=(std::move(x.m_timing_values));
		m_overhead_values.operator=(std::move(x.m_overhead_values));
		m_delta_values.operator=(std::move(x.m_delta_values));
		
	return (*this);
}

void 
gms::system::GTC64Timer::run_benchmark() {

	volatile ULONGLONG STARTF{0Ui64}, STOPF{0Ui64},
	                   STARTW{0Ui64}, STOPW{0Ui64};
   //  Wrapper function test measurement loop 
	for (std::size_t i = 0Ui64; i != m_nruns; ++i) {
		for (std::size_t j = 0Ui64; j != m_datum_size; ++j) {

			STARTW = ::GetTickCount64();
			m_ptwrap(); // indirect call through register
			STOPW = ::GetTickCount64();

			if (STOPW - STARTW <= 0Ui64) {
				std::printf("*** Error *** in wrapper timing loop: invalid timing at i:%llu,j:%llu, STARTW=%llu, STOPW=%llu, DELTAW=%llu\n",
				            i,j,STOPW,STARTW,(STOPW-STARTW));
				m_overhead_values.operator[](Ix2D(i,m_datum_size,j)) = 0Ui64;
			}
			else {
				m_overhead_values.operator[](Ix2D(i, m_datum_size, j)) = (STOPW-STARTW);
			}
		}
	}
	// Now timed function test measurement loop
	for (std::size_t i = 0Ui64; i != m_nruns; ++i) {
		for (std::size_t j = 0Ui64; j != m_datum_size; ++j) {

		   
			 STARTF = ::GetTickCount64();
			 m_ptfunc(); // indirect call through register
			 STOPF = ::GetTickCount64();
			
			 if ((STOPF - STARTF) <= 0Ui64) {
				 std::printf(" *** Error *** in timed function loop invalid read out at: i:%llu,j:%llu, STARTF=%llu, STOPF=%llu, DELTAF=%llu\n",
					 i, j, STARTF, STOPF, (STOPF-STARTF));
				 m_timing_values.operator[](Ix2D(i,m_datum_size,j)) = 0Ui64;
				
			 }
			 else {
				 m_timing_values.operator[](Ix2D(i,m_datum_size,j)) = STOPF - STARTF;
				
			 }
		}
	}

	// Compute delta (do not skip over potentially present zero values)
	for (std::size_t i = 0Ui64; i != m_nruns; ++i) {
		for (std::size_t j = 0Ui64; j != m_datum_size; ++j) {
			if (m_timing_values[i] - m_overhead_values[i] > 0Ui64) {
				m_delta_values[i] = m_timing_values[i] - m_overhead_values[i];
			}
			else {
				m_delta_values[i] = 0Ui64;
			}
		}
	}
}

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
	
bool 
gms::system
::GTC64Timer
::compute_stats(_Out_ double &mean,
			    _Out_ double &adev,
			    _Out_ double &sdev,
				_Out_ double &skew,
				_Out_ double &kurt) {
	using namespace gms::math::constants;
	if (m_delta_values.size() <= lo_bound) {
		mean = dinf, adev = dinf,
		sdev = dinf, skew = dinf,
		kurt = dinf;
		return (false);
	}
	ULONGLONG s{ 0Ui64 }, prev{ 0Ui64 };
	std::size_t vlen{};
	__declspec(align(64)) struct {
		double len{}, sum{}, var{}, t{}, t2{}, t3{},
		t4{}, t5{}, isdev{}, tmp{};
	}r8_loc;
	vlen = m_delta_values.size();
	r8_loc.len = __binary64_from_uint64(vlen);
	// Compute mean guard against an overflow
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (s < prev) return (false);
		prev = s;
		s += m_delta_values[i];
	}
	mean =  __binary64_from_uint64(s) / r8_loc.len;
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
		std::printf(" GTC64Timer::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
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
			r8_loc.t3 = __binary64_mul_binary64_binary64(
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
		    std::printf(" GTC64Timer::compute_stats: Invalid kurtosis: %.15f\n", kurt);
		    mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	return (true);
}

#else

bool
gms::system::GTC64Timer
::compute_stats(_Out_ double &mean,
				_Out_ double &adev,
				_Out_ double &sdev,
				_Out_ double &skew,
				_Out_ double &kurt) {

	using namespace gms::math::constants;
	if (m_timing_values.size() <= lo_bound) {
		mean = dinf, adev = dinf,
		sdev = dinf, skew = dinf,
		kurt = dinf;
		return (false);
	}
	ULONGLONG s{}, prev{};
	__declspec(align(64)) struct {
		double len{}, sum{}, var{}, t{},
		t2{}, t3{}, t4{}, t5{}, isdev{},
		fracp{}, ct2{}, tmp{};
	} r8_loc;
	r8_loc.len = static_cast<double>(m_delta_values.size());
	// Compute mean (guard against an overflow)
	for (size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		if (s < prev) return (false);
		prev = s;
		s += m_delta_values.operator[](i);
	}
	r8_loc.sum = static_cast<double>(s);
	mean = r8_loc.sum / r8_loc.len;
	//  Compute average deviation and variance
	for (size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		if (zero != m_delta_values[i]) {
			r8_loc.tmp = static_cast<double>(m_delta_values[i]);
			r8_loc.t = std::abs(r8_loc.tmp - mean);
			adev += r8_loc.t; //Potential overflow?
			r8_loc.t2 = (r8_loc.tmp - mean) * (r8_loc.tmp - mean);
			r8_loc.ct2 = r8_loc.t2;
			r8_loc.fracp = r8_loc.ct2 - static_cast<ULONGLONG>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
			   std::printf(" GTC64Timer::compute_stats: Losing a significant digits: %.16f, result may not be an accurate\n", r8_loc.fracp);
			}
			r8_loc.var += r8_loc.t2; // potential overflow?
		}
	}
	adev /= r8_loc.len;
	if (r8_loc.var <= 0.0) {
		std::printf(" GTC64Timer::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		return (false);
	}
	r8_loc.var /= r8_loc.len;
	sdev = std::sqrt(r8_loc.var);
	r8_loc.isdev = 1.0 / sdev;
	r8_loc.fracp = -1.0;
	r8_loc.tmp = 0Ui64;
	for (size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		if (zero != m_delta_values.operator[](i)) {
			r8_loc.tmp = static_cast<double>(m_delta_values[i]);
			r8_loc.t3 = (r8_loc.tmp - mean) * r8_loc.isdev;
			r8_loc.ct2 = r8_loc.t3;
			r8_loc.fracp = r8_loc.t3 - static_cast<ULONGLONG>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" GTC64Timer::compute_stats: Losing a significant digits: %.16f, result may not be an accurrate\n", r8_loc.fracp);
			}
			// potential overflows?
			r8_loc.t4 = r8_loc.t3*r8_loc.t3*r8_loc.t3;
			skew += r8_loc.t4;
			r8_loc.t5 = r8_loc.t4 * r8_loc.t3;
			kurt += r8_loc.t5;
		}
	}
	skew /= r8_loc.len;
	kurt /= r8_loc.len - 3.0;
	if (kurt < 1.0) {
		std::printf(" GTC64Timer::compute_stats: Invalid kurtosis: %.15f\n", kurt);
		return (false);
	}
	return (true);
}

#endif

void
gms::system::
GTC64Timer
::print() const {
	
	std::cout << "	Dumping state of  " << typeid(*this).raw_name() << "\n";
	std::cout << "============================================================\n";
	std::cout << " Number of test runs:      " << m_nruns << "\n"
		      << " Datum size per run:       " << m_datum_size << "\n"
			  << " Tested function address:  " << std::hex << m_ptfunc << "\n"
			  << " Tested wrapper  address:  " << std::hex << m_ptwrap << "\n"
			  << " Tested function name:     " << m_func_name.data() << "\n"
			  << " Located in file name:     " << m_file_name.data() << "\n";
	std::cout << "======================= Printing data set =========================\n";
	for (size_t i = 0Ui64; i != m_nruns; ++i) {
		std::cout << " Function sample  :  Wrapper sample  :  Delta  samples         \n";
		for (size_t j = 0Ui64; j != m_datum_size; ++j) {
			std::cout << std::setw(4)  << m_timing_values[Ix2D(i, m_datum_size, j)] << "\n"
					  << std::setw(8)  << m_overhead_values[Ix2D(i, m_datum_size, j)] << "\n"
				      << std::setw(12) << m_delta_values[Ix2D(i, m_datum_size, j)] << "\n";
		}
		std::cout << "\n";
	}
}

#if !defined (GMS_GTC64_TIMER_CHECK_ARGS_CONFORMANCE)
#define GMS_GTC64_TIMER_CHECK_ARGS_CONFORMANCE(vres,set1,set2)	                    \
		do{														                    \
              if ((vres).size() != (set1).m_delta_values.size() ||                  \
	               (set1).m_delta_values.size() != (set2).m_delta_values.size()) {  \
				   return (false);													\
              }																		\
		} while (0);
#endif

bool
gms::system
::GTC64Timer
::delta_values_eq(_Out_ std::vector<bool> &vres,
				  _In_ const GTC64Timer &set1,
				  _In_ const GTC64Timer &set2) {
	GMS_GTC64_TIMER_CHECK_ARGS_CONFORMANCE(vres, set1, set2)
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (set1.m_delta_values[i] == set2.m_delta_values[i]) {
			vres.operator[](i) = true;
		}
		else {
			vres.operator[](i) = false;
		}
	}
	return (false);
}

bool
gms::system
::GTC64Timer
::delta_values_ineq(_Out_ std::vector<bool> &vres,
					_In_ const GTC64Timer &set1,
					_In_ const GTC64Timer &set2) {
	GMS_GTC64_TIMER_CHECK_ARGS_CONFORMANCE(vres, set1, set2)
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (set1.m_delta_values[i] != set2.m_delta_values[i]) {
			vres.operator[](i) = true;
		}
		else{
			vres.operator[](i) = false;
		}
	}
	return (true);
}

bool
gms::system
::GTC64Timer
::delta_values_gt(_Out_ std::vector<bool> &vres,
				  _In_ const GTC64Timer &set1,
				  _In_ const GTC64Timer &set2) {
	GMS_GTC64_TIMER_CHECK_ARGS_CONFORMANCE(vres, set1, set2)
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (set1.m_delta_values[i] > set2.m_delta_values[i]) {
			vres.operator[](i) = true;
		}
		else {
			vres.operator[](i) = false;
		}
	}
	return (true);
}

bool
gms::system
::GTC64Timer
::delta_values_lt(_Out_ std::vector<bool> &vres,
				  _In_ const GTC64Timer &set1,
				  _In_ const GTC64Timer &set2){
	GMS_GTC64_TIMER_CHECK_ARGS_CONFORMANCE(vres, set1, set2)
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (set1.m_delta_values[i] < set2.m_delta_values[i]) {
			vres.operator[](i) = true;
		}
		else {
			vres.operator[](i) = false;
		}
	}
	return (true);
}
