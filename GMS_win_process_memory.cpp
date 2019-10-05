
#include <Psapi.h>
#include "../GMS_config.h"
#include "GMS_win_process_memory.h"

//
//	Implementation
//

gms::system::ProcessMemInfo::
ProcessMemInfo(_In_ const int nstates, _In_ DWORD dwDesiredAccess,
_In_ BOOL bInheritHandle)
:
m_nstates{ nstates },
m_procID{ GetCurrentProcessId() },
m_phandle{ OpenProcess(dwDesiredAccess, bInheritHandle, m_procID) } {
	if (NULL == m_phandle) {
		PROCESS_MEMORY_INFO_ERROR("ProcessMemoryInfo::ProcessMemoryInfo(int,DWORD,BOOL) -- !!! Failed to open process !!!")
	}
	
	m_pmem.m_cb = 0;
	m_pmem.m_PageFaultCount = 0;
	m_pmem.m_PeakWorkingSetSize = 0ULL;
	m_pmem.m_WorkingSetSize     = 0ULL;
	m_pmem.m_QuotaPeakPagedPoolUsage = 0ULL;
	m_pmem.m_QuotaPagedPoolUsage = 0ULL;
	m_pmem.m_QuotaPeakNonPagedPoolUsage = 0ULL;
	m_pmem.m_QuotaNonPagedPoolUsage = 0ULL;
	m_pmem.m_PageFileUsage = 0ULL;
	m_pmem.m_PeakPageFileUsage = 0ULL;
	m_mlist = std::list<proc_mem>{};
}

gms::system::ProcessMemInfo::
~ProcessMemInfo() {
	
}

int gms::system::ProcessMemInfo::get_nstates() const {
	return (this->m_nstates);
}

HANDLE gms::system::ProcessMemInfo::get_phandle() const {
	return (this->m_phandle);
}

DWORD gms::system::ProcessMemInfo::get_procID() const {
	return (this->m_procID);
}

const gms::system::ProcessMemInfo::memlist & 
gms::system::ProcessMemInfo::get_mlist() const {
	return (this->m_mlist);
}

gms::system::ProcessMemInfo::memlist::const_iterator
gms::system::ProcessMemInfo::get_mlist_cbegin() const {
	return (this->m_mlist.cbegin());
}

gms::system::ProcessMemInfo::memlist::const_iterator
gms::system::ProcessMemInfo::get_mlist_cend() const {
	return (this->m_mlist.cend());
}

bool
gms::system::ProcessMemInfo::add_snapshot() {
	BOOL bOK = FALSE;
	PROCESS_MEMORY_COUNTERS pmc;
	bOK = GetProcessMemoryInfo(this->m_phandle,&pmc,sizeof(pmc));
	if (NULL == bOK) {
		PRINT_ERROR_VALUE("ProcessMemInfo::add_snapshot: !!! GetProcessMemoryInfo failed !!!",GetLastError())
		CloseHandle(this->m_phandle);
		return (false);
	}
	this->m_pmem.m_cb = pmc.cb;
	this->m_pmem.m_PageFaultCount = pmc.PageFaultCount;
	this->m_pmem.m_PeakWorkingSetSize = pmc.PeakWorkingSetSize;
	this->m_pmem.m_WorkingSetSize = pmc.WorkingSetSize;
	this->m_pmem.m_QuotaPeakPagedPoolUsage = pmc.QuotaPeakPagedPoolUsage;
	this->m_pmem.m_QuotaPagedPoolUsage = pmc.QuotaPagedPoolUsage;
	this->m_pmem.m_QuotaPeakNonPagedPoolUsage = pmc.QuotaPeakNonPagedPoolUsage;
	this->m_pmem.m_QuotaNonPagedPoolUsage = pmc.QuotaNonPagedPoolUsage;
	this->m_pmem.m_PageFileUsage = pmc.PagefileUsage;
	this->m_pmem.m_PeakPageFileUsage = pmc.PeakPagefileUsage;
	this->m_mlist.push_back(this->m_pmem);
	CloseHandle(this->m_phandle);
	return (true);
}

void
gms::system::ProcessMemInfo::print_state() const {
	if (!this->m_mlist.empty()) {
	    int ordinal = 0;
		std::cout << "======================================================\n"
			      << " Data collected for:    " << "\n"
			      << " Process handle:        " << this->m_phandle << "\n"
			      << " Process ID:            " << this->m_procID << "\n"
				  << " Number of snapshots:   " << this->m_mlist.size() << "\n"
			      << "=======================================================\n";
		for (auto it = this->m_mlist.cbegin(); it != this->m_mlist.cend(); ++it) {
		     ++ordinal;
			 std::cout << "==========================================================================================\n"
			           << " Snapshot #: " << ordinal  << "\n"
			           << " Size of collected data:      "  << it->m_cb  << " bytes." << "\n"
				       << " Page fault count:            "  << it->m_PageFaultCount << "\n"
					   << " PeakWorkingSetSize:          "  << it->m_PeakWorkingSetSize << " bytes. \n" 
					   << " Working set size:            "  << it->m_WorkingSetSize << " bytes. \n"
					   << " Quota peak paged pool usage: "  << it->m_QuotaPeakPagedPoolUsage << " bytes. \n"
					   << " Quota paged pool useage:     "  << it->m_QuotaPagedPoolUsage << " bytes. \n"
					   << " Quota peak non paged pool usage: " << it->m_QuotaPeakNonPagedPoolUsage << " bytes. \n"
					   << " Quota non paged pool usage:      " << it->m_QuotaNonPagedPoolUsage << " bytes. \n"
					   << " Pagefile usage:					 " << it->m_PageFileUsage << " bytes. \n"
					   << " Peagefile peak usage:            " << it->m_PeakPageFileUsage << " bytes. \n"
					   << "=============================================================================================\n";

           
		}
		
	}
}
