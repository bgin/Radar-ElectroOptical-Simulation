
#ifndef __GMS_WIN_UTILS_H__
#define __GMS_WIN_UTILS_H__



#include <cstdint>
#include <Windows.h>

namespace file_info {
  const unsigned int gGMS_WIN_UTILS_MAJOR = 1U;
  const unsigned int gGMS_WIN_UTILS_MINOR = 0U;
  const unsigned int gGMS_WIN_UTILS_MICRO = 0U;
  const unsigned int gGMS_WIN_UTILS_FULLVER =
    1000U*gGMS_WIN_UTILS_MAJOR + 100U*gGMS_WIN_UTILS_MINOR + 10U*gGMS_WIN_UTILS_MICRO;
  const char * const pgGMS_WIN_UTILS_CREATE_DATE = "05-10-2019 15:28 +00200 (SAT 05 OCT 2019 GMT+2)";
  const char * const pgGMS_WIN_UTILS_BUILD_DATE  = "00-00-0000 00:00";
  const char * const pgGMS_WIN_UTILS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_WIN_UTILS_SYNOPSIS    = "Collection of small Windows utility functions";
} 

namespace gms {
	namespace system {
		namespace free_funcs {


			
			uint64_t   cpu_freq_approximated();

			bool       set_cpu_affinity(_In_ const uint32_t);

			//
			// Process/thread info
			//

			//@function: wset_thread_ideal_processor
			//@Remark:
			//			 Ideal processor is set for
			//           current thread.
			//@Calls:
			//         SetThreadIdealProcessorEx
			bool       set_thread_ideal_processor();

			// Same as above with the only difference
			// a HANDLE parameter passed
			bool       set_thread_ideal_processor(_In_ const HANDLE,
												  _In_ const WORD);

			//@function: wset_thread_priority
			//@Calls:
			//		   SetThreadPriority
			bool       set_thread_priority(_In_ const HANDLE,
										   _In_ const int32_t);

			//@function: wget_current_process
			//@Calls:
			//		   GetCurrentProcess
			HANDLE     get_current_process();

			//@function: wget_current_processid
			//@Calls:
			//		   GetCurrentProcessId
			DWORD      get_current_processid();

			//@function: wget_current_thread
			//@Calls:
			//		     GetCurrentThread
			HANDLE     get_current_thread();

			//@function: wget_current_threadid
			//@Calls:
			//           GetCurrentThreadId
			DWORD      get_current_threadid();


				
		}
	}
}




#endif /*__GMS_WIN_UTILS_H_*/
