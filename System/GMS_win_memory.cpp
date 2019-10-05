

#include "GMS_win_memory.h"
#include "../GMS_config.h"
#include "../GMS_error_macros.h"
#include <DbgHelp.h>

std::pair<SIZE_T,SIZE_T> 
gms::system::helpers::
get_working_set_size(gms::system::helpers::OpenProcessArgs * pArgs) {
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(NULL != pArgs);
#else
	if(NULL == pArgs) {
		PRINT_ERROR_INFO("get_working_set_size: -- Invalid argument !!!")
	    return (make_pair(0ULL,0ULL));
	}
#endif
	SIZE_T dwMin{},dwMax{};
	HANDLE hProc;
	BOOL bOK = FALSE;
	hProc = OpenProcess(pArgs->dwDesiredAccess,pArgs->bInheritHandle,pArgs->dwProcessId);
	if (NULL == hProc) {
		PRINT_ERROR_VALUE("get_working_set_size: -- OpenProcess failed !!!", GetLastError())
	    return (std::make_pair(0ULL,0ULL));
	}
	bOK = GetProcessWorkingSetSize(hProc,&dwMin,&dwMax);
	if (bOK == FALSE) {
		PRINT_ERROR_VALUE("get_working_set_size: -- GetProcessWorkingSetSize failed !!!", GetLastError())
			return (std::make_pair(0ULL,0ULL));
	}
	return (std::make_pair(dwMin/KiBytes,dwMax/KiBytes));
}

bool
gms::system::helpers::
monitor_working_set_size(_In_ gms::system::helpers::OpenProcessArgs * pArgs,
_Inout_ VecSamples &vsamples, _In_ const UINT32 col_policy ) {
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(NULL != pArgs && vsamples.size() > 10ULL);
#else
	if(NULL == pArgs || vsamples.size() <= 10 ULL) {
		PRINT_ERROR_INFO("monitor_working_set_size: -- Invalid arguments !!!")
		return (false);
	}
#endif
  using idx = VecSamples::size_type;
  DWORD minmsec = 1UL;
  const UINT ms = 1U;
  SIZE_T dwMin{}, dwMax{};
  HANDLE hProc;
  BOOL bOK = FALSE;
  TIMECAPS time;
  MMRESULT mres;
  mres = timeGetDevCaps(&time, sizeof(time));
  
  if (mres != MMSYSERR_NOERROR) {
	  PRINT_ERROR_VALUE("monitor_working_set_size: -- timeGetDevCaps failed -- inaccurate measurements !!!",mres)
	  return (false);
  }
  if (minmsec < time.wPeriodMin || minmsec > time.wPeriodMax )
	  mres = timeBeginPeriod(ms);
  if (mres != TIMERR_NOERROR) {
	  PRINT_ERROR_VALUE("monitor_working_set_size: -- timeBeginPeriod failed !!!", mres)
	  return (false);
  }
  

  hProc = OpenProcess(pArgs->dwDesiredAccess,pArgs->bInheritHandle,pArgs->dwProcessId);
  if (NULL == hProc) {
	  PRINT_ERROR_VALUE("monitor_working_set_size: -- OpenProcess failed !!!", GetLastError())
	  return (false);
  }
 
  switch (col_policy) {

  case 0: { // Stop collecting when samples do not vary in time.

         for (idx i = 0ULL; i != vsamples.size(); ++i) {
	          bOK = GetProcessWorkingSetSize(hProc,&dwMin,&dwMax);
	          if (bOK == FALSE) {
		          PRINT_ERROR_VALUE("monitor_working_set_size: -- GetProcessWorkingSetSize failed !!!", GetLastError())
		          return (false);
	      }
	       vsamples.operator[](i).operator=(std::make_pair(dwMin/KiBytes,dwMax/KiBytes));
	       if (i > 0ULL && i < vsamples.size()) {
		        SIZE_T dmin = vsamples.operator[](i).first - vsamples.operator[](i-1ULL).first;
		        SIZE_T dmax = vsamples.operator[](i).second - vsamples.operator[](i-1ULL).second;
		  if (dmin == 0ULL || dmax == 0ULL) {
			    std::cerr << "Process: " << pArgs->dwProcessId << " Working set size not varying in time... -- terminating measurement!! \n";
			    mres = timeEndPeriod(ms);
			    return (false);
		  }
	  }
	  if (i > 2ULL) {
	      minmsec += static_cast<DWORD>(i);
	      Sleep(minmsec);
	  }
  }

}
  case 1: {
       // Collect samples without checking time variation.
			  for (idx i = 0ULL; i != vsamples.size(); ++i) {
				  bOK = GetProcessWorkingSetSize(hProc,&dwMin,&dwMax);
				  if (bOK == FALSE) {
					  PRINT_ERROR_VALUE("monitor_working_set_size: -- GetProcessWorkingSetSize failed !!!", GetLastError())
					  return (false);
				  }
				  vsamples.operator[](i).operator=(std::make_pair(dwMin/KiBytes,dwMax/KiBytes));
				  if (i > 2ULL) {
					  minmsec += static_cast<DWORD>(i);
					  Sleep(minmsec);
				  }
			  }
  }
  default: {
			   PRINT_ERROR_INFO("monitor_working_set_size: -- Invalid parameter to switch statement!!!")
			   mres = timeEndPeriod(ms);
			   return (false);
  }
}
  mres = timeEndPeriod(ms);
  if (mres != TIMERR_NOERROR) {
	  PRINT_ERROR_VALUE("monitor_working_set_size: -- timeBeginPeriod failed !!!", mres)
	  return (false);
  }
  return (true);
}

SIZE_T
gms::system::helpers::
working_set_size_delta(_In_ OpenProcessArgs * pArgs) {
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(NULL != pArgs);
#else
	if(NULL == pArgs) {
		PRINT_ERROR_INFO("working_set_size_delta: -- Invalid argument!!!")
		return (0ULL);
	}
#endif
	auto delta = get_working_set_size(pArgs);
	if (0ULL != delta.first && 0ULL != delta.second)
		return (delta.second-delta.first);
	else
	    return (0ULL);
}

bool
gms::system::helpers::
working_set_samples_delta(_In_ const VecSamples &vsamples,
						  _Out_ VecSamples &delta) {

	if (vsamples.empty() || (vsamples.size() != delta.size()))
		return (false);
	SIZE_T del_lo{}, del_hi{};
	using idx = VecSamples::size_type;
	for (idx i = 1; i != vsamples.size(); ++i) {
		del_lo = vsamples.operator[](i).first  - vsamples.operator[](i-1ULL).first;
		del_hi = vsamples.operator[](i).second - vsamples.operator[](i-1ULL).second;
		if (0ULL != del_lo && 0ULL != del_hi)
		    delta.operator[](i).operator=(std::make_pair(del_lo,del_hi));
	}
	return (true);
}

