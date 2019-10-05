
#include <cstdio>
#include <intrin.h>
#include "GMS_win_utils.h"
#include "../GMS_config.h"

uint64_t
gms::system::free_funcs
::cpu_freq_approximated(){
	// Crude approximation of current CPU frequency. On failure return 0ULL frequency.
	// Get system clock time interval
	MMRESULT sys_clock_res = 9999;
	TIMECAPS timecaps;
	ZeroMemory(&timecaps, sizeof(timecaps));
	sys_clock_res = ::timeGetDevCaps(&timecaps, sizeof(timecaps));
#if (GMS_DEBUG_ON) == 1
	std::printf("timeGetDevCaps returned: min clock period: %d, max clock period: %d\n", 
	            timecaps.wPeriodMin, timecaps.wPeriodMax);
#endif
	if (0ULL > sys_clock_res) {
#if (GMS_DEBUG_ON) == 1
		std::printf("timeGetDevCaps failed with an error: %d\n", sys_clock_res);
#endif
		return (0ULL);
	}
	int32_t dummy1[4] = {}, dummy2{};
	volatile uint64_t start{}, stop{};
	__cpuid(&dummy1[0],dummy2);
	start = __rdtsc();
	::Sleep(timecaps.wPeriodMin);
	stop = __rdtsc();
	__cpuid(&dummy1[0],dummy2);
	// Guard against negative/wraparound value.
	if (stop > start)
		return (stop - start);
    else
		return (0ULL);
}

bool
gms::system::free_funcs
::set_cpu_affinity(_In_ const uint32_t mask) {
	
	SYSTEM_INFO system_info;
	::GetSystemInfo(&system_info);
	DWORD ncpus{system_info.dwNumberOfProcessors};
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(mask >= 0 && mask <= ncpus);
#else
	if (mask < 0 || mask > ncpus) {
		return (false);
	}
#endif
	HANDLE cthread = ::GetCurrentThread(); 
	DWORD_PTR dw_result = ::SetThreadAffinityMask(cthread,mask);
	if (0 == dw_result) {
		std::printf(" SetThreadAffinityMask: -- failed with an error: -- 0x%x\n", ::GetLastError());
		return (false);
	}
	return (true);
}

bool
gms::system::free_funcs
::set_thread_ideal_processor() {

	bool bOk{};
	HANDLE thand = ::GetCurrentThread();
	PROCESSOR_NUMBER pnum;
	::GetCurrentProcessorNumberEx(&pnum);
	bOk = ::SetThreadIdealProcessorEx(thand,&pnum,NULL);
	if (!bOk) {
		std::printf(" SetThreadIdealProcessorEx: -- failed with an error: -- 0x%x\n", ::GetLastError());
		return (false);
	}
	return (true);
}

bool
gms::system::free_funcs
::set_thread_ideal_processor(_In_ const HANDLE thand,
						     _In_ const WORD number) {
	if (0 > number) {
		return (false);
	}
	bool bOk{};
	PROCESSOR_NUMBER pnum;
	pnum.Group = 0;
	pnum.Number = number;
	bOk = ::SetThreadIdealProcessorEx(thand,&pnum,NULL);
	if (!bOk) {
		std::printf(" SetThreadIdealProcessorEx: -- failed with an error: -- 0x%x\n", ::GetLastError());
		return (false);
	}
	return (true);
}

bool
gms::system::free_funcs
::set_thread_priority(_In_ const HANDLE thand,
				      _In_ const int32_t priority) {

	bool bOk{};
	bOk = ::SetThreadPriority(thand,priority);
	if (!bOk) {
		std::printf(" SetThreadPriority: -- failed with an error: -- 0x%x\n", ::GetLastError());
		return (false);
	}
	return (true);
}

HANDLE
gms::system::free_funcs::
get_current_process() {
	return (::GetCurrentProcess());
}

DWORD
gms::system::free_funcs::
get_current_processid() {
	return (::GetCurrentProcessId());
}

HANDLE
gms::system::free_funcs
::get_current_thread() {
	return (::GetCurrentThread());
}

DWORD
gms::system::free_funcs
::get_current_threadid() {
	return (::GetCurrentThreadId());
}
