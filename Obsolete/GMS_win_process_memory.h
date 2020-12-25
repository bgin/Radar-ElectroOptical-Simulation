
#ifndef __GMS_WIN_PROCESS_MEMORY_H__
#define __GMS_WIN_PROCESS_MEMORY_H__



#include <Windows.h>
#include <list>
#include <Psapi.h>

namespace file_info {
  const unsigned int gGMS_WIN_PROCESS_MEMORY_MAJOR = 1U;
  const unsigned int gGMS_WIN_PROCESS_MEMORY_MINOR = 0U;
  const unsigned int gGMS_WIN_PROCESS_MEMORY_MICRO = 0U;
  const unsigned int gGMS_WIN_PROCESS_MEMORY_FULLVER =
    1000U*gGMS_WIN_PROCESS_MEMORY_MAJOR+100U*gGMS_WIN_PROCESS_MEMORY_MINOR+10U*gGMS_WIN_PROCESS_MEMORY_MICRO;
  const char * const pgGMS_WIN_PROCESS_MEMORY_CREATE_DATE = "05-10-2019 17:54 +00200 (SAT 05 OCT 2019 GMT+2)";
  const char * const pgGMS_WIN_PROCESS_MEMORY_BUILD_DATE  = "00-00-0000 00:00";
  const char * const pgGMS_WIN_PROCESS_MEMORY_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_WIN_PROCESS_MEMORY_SYNOPSIS    = "Process memory monitoring functionality.";
}

namespace gms {
	namespace system {

			//
			// Declare class ProcessMemInfo which will collect
			// varying time snapshots of memory process usage.
			//

#if !defined (PROCESS_MEMORY_INFO_ERROR)
#define PROCESS_MEMORY_INFO_ERROR(msg)  \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "," << "system error value: " << std::hex << "0x" << GetLastError() << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	throw std::runtime_error(msg);
#endif

		class ProcessMemInfo {
		   
			private:
				
				int		m_nstates;

				DWORD   m_procID;
				HANDLE  m_phandle;

				struct proc_mem {
					
					DWORD m_cb;
					DWORD m_PageFaultCount;
					SIZE_T m_PeakWorkingSetSize;
					SIZE_T m_WorkingSetSize;
					SIZE_T m_QuotaPeakPagedPoolUsage;
					SIZE_T m_QuotaPagedPoolUsage;
					SIZE_T m_QuotaPeakNonPagedPoolUsage;
					SIZE_T m_QuotaNonPagedPoolUsage;
					SIZE_T m_PageFileUsage;
					SIZE_T m_PeakPageFileUsage;

				}m_pmem;

				using memlist = std::list<proc_mem>;

				
				memlist m_mlist;

				public:

					ProcessMemInfo() = delete;

				        ProcessMemInfo(_In_ const int, _In_ DWORD, _In_ BOOL);

					ProcessMemInfo(_In_ const ProcessMemInfo &) = delete;

					~ProcessMemInfo() noexcept(true);

					// Class member operators (deleted)

					ProcessMemInfo & operator=(_In_ const ProcessMemInfo &) = delete;

					// Accessors

					int get_nstates() const;

					HANDLE get_phandle() const; // Warning: HANDLE will be closed upon object destruction

					DWORD  get_procID() const;

					const memlist & get_mlist() const; // Warning:  beware of object lifetime here!!

					memlist::const_iterator get_mlist_cbegin() const;

					memlist::const_iterator get_mlist_cend() const;

					// Mutators

					bool add_snapshot();

					void print_state() const;
		};


	}
}

#endif /*__GMS_WIN_PROCESS_MEMORY_H__*/
