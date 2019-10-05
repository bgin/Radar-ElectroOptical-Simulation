
#ifndef __GMS_WIN_MEMORY_H__
#define __GMS_WIN_MEMORY_H__



#include <utility>
#include <Windows.h>
#include <vector>

namespace file_info {
  const unsigned int gGMS_WIN_MEMORY_MAJOR = 1U;
  const unsigned int gGMS_WIN_MEMORY_MINOR = 0U;
  const unsigned int gGMS_WIN_MEMORY_MICRO = 0U;
  const unsigned int gGMS_WIN_MEMORY_FULLVER =
    1000U*gGMS_WIN_MEMORY_MAJOR+100U*gGMS_WIN_MEMORY_MINOR+10U*gGMS_WIN_MEMORY_MICRO;
  const char * const pgGMS_WIN_MEMORY_CREATE_DATE = "05-10-2019 16:32 +00200 (SAT 05 OCT 2019 GMT+2)";
  const char * const pgGMS_WIN_MEMORY_BUILD_DATE  = "00-00-0000 00:00";
  const char * const pgGMS_WIN_MEMORY_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_WIN_MEMORY_SYNOPSIS    = "System memory related procedures.";
}

namespace lam {
	namespace system {
		namespace helpers {

			static const SIZE_T KiBytes = 1024ULL;

			 struct OpenProcessArgs {
				DWORD dwDesiredAccess;
				BOOL  bInheritHandle;
				DWORD dwProcessId;
			};

			using VecSamples = std::vector<std::pair<SIZE_T,SIZE_T>>;
			//
			// Get Process working set size i.e. minimum and maximum in kibytes.
			//
			std::pair<SIZE_T, SIZE_T> 
			get_working_set_size(_In_ OpenProcessArgs *);

			//
			// Get process working set size i.e. minimum and maximum in kibytes.
			// Attempt to monitor usage of aforementioned metrics.
			//
			bool
			monitor_working_set_size(_In_ OpenProcessArgs *,
			                         _Inout_ VecSamples &, const UINT32 );

		    //
			// Analyze collected samples of working set size.
			//
			bool 
			working_set_samples_delta(_In_ const VecSamples &,
										_Out_ VecSamples &);
			//
			// Return process working set delta (max-min).
			//
			SIZE_T 
			working_set_size_delta(_In_ OpenProcessArgs *);

			//
			// StackWalk64Wrapper
			//
			BOOL
			StackWalk64Wrapper(_In_ const HANDLE, _In_ const HANDLE);
		}
		
	}
}



#endif /*__GMS_WIN_MEMORY_H__*/
