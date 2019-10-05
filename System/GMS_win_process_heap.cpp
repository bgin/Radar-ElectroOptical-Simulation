

#include <Psapi.h>
#include "../GMS_config.h"
#include "GMS_win_process_heap.h"
//
//	Implementation
//

gms::system::ProcessHeapInfo::
ProcessHeapInfo(_In_ const int nstates, _In_ const WORD flags)
:
m_nstates{ nstates },
m_flags{ flags },
m_heapHandle{},
m_procHandle{},
m_procID{},
m_lock_heap{} {
	m_heap.m_lpData = NULL;
	m_heap.m_cbData = 0;
	m_heap.m_cbOverhead = 0;
	m_heap.m_iRegionIndex = 0;
	m_heap.m_wFlags = m_flags;
	m_heap.m_Block.m_hMem = 0;
	m_heap.m_Block.m_dwReserved[3] = 0;
	m_heap.m_Region.m_dwCommitedSize = 0;
	m_heap.m_Region.m_dwUncommitedSize = 0;
	m_heap.m_Region.m_lpFirstBlock = NULL;
	m_heap.m_Region.m_lpLastBlock = NULL;
	m_heapList = std::list<heap>{};
}

gms::system::ProcessHeapInfo::
~ProcessHeapInfo(){
	if (!m_heapList.empty()) {
	    for (std::list<heap>::iterator it = m_heapList.begin(); it != m_heapList.end(); ++it) {
			if (it->m_Region.m_lpFirstBlock != NULL)
		          it->m_Region.m_lpFirstBlock = NULL;
			if (it->m_Region.m_lpLastBlock != NULL)
			      it->m_Region.m_lpLastBlock = NULL;
	     }
	}
}

int gms::system::ProcessHeapInfo::get_nstates() const {
	return (this->m_nstates);
}

WORD gms::system::ProcessHeapInfo::get_flags() const {
	return (this->m_flags);
}

HANDLE gms::system::ProcessHeapInfo::get_heapHandle() const {
	return (this->m_heapHandle);
}

HANDLE gms::system::ProcessHeapInfo::get_procHandle() const {
	return (this->m_procHandle);
}

DWORD gms::system::ProcessHeapInfo::get_procID() const {
	return (this->m_procID);
}

bool  gms::system::ProcessHeapInfo::get_heap_lock() const {
	return (this->m_lock_heap);
}

// beware of object lifetime here!!
// Reference maybe dangling when object will be destroyed.
const std::list<gms::system::ProcessHeapInfo::heap>  & 
gms::system::ProcessHeapInfo::get_heapList() const {
	return (this->m_heapList);
}

 std::list<gms::system::ProcessHeapInfo::heap>::const_iterator
gms::system::ProcessHeapInfo::get_heapListIterator() const {
	return (this->m_heapList.cbegin());
}

void gms::system::ProcessHeapInfo::print_heap_state() const{
	if (!this->m_heapList.empty()) {
		std::cout << "========================================================================\n"
		          <<    "Process Handle: " << this->m_procHandle << "\n"
			      <<    "Heap handle:    " << std::hex << this->m_heapHandle << "\n"
				  <<    "Process ID:     " << this->m_procID << "\n";
		std::cout << "=========================================================================\n";
		for (auto it = this->m_heapList.cbegin(); it != this->m_heapList.cend(); ++it) {
			std::cout << "Data at address:      " << std::hex << "0x" << it->m_lpData << "\n"
			          << "Size of data:         " << it->m_cbData << " bytes. \n"
					  << "Size of overhead:     " << it->m_cbOverhead << " bytes. \n"
					  << "iRegion index:        " << it->m_iRegionIndex << "\n"
					  << "Region commited:      " << it->m_Region.m_dwCommitedSize << " bytes. \n"
					  << "Region uncommited:    " << it->m_Region.m_dwUncommitedSize << " bytes. \n"
					  << "First block at:       " << std::hex << "0x" << it->m_Region.m_lpFirstBlock << "\n"
					  << "Last block  at:       " << std::hex << "0x" << it->m_Region.m_lpLastBlock << "\n";
		}
		std::cout << "=======================================================================================\n";
	}
}

void  gms::system::ProcessHeapInfo::heap_enumerate() {
	DWORD hwerr = 0;
	BOOL bhl = false;
	BOOL bhw = false;
	PROCESS_HEAP_ENTRY entry;
	this->m_heapHandle = GetProcessHeap();
	if (NULL == this->m_heapHandle) {
		PROCESS_HEAP_INFO_ERROR("ProcessHeapInfo::heap_enumerate -- !!! GetProcessHeap failed !!!")
	}
	this->m_procHandle = GetCurrentProcess();
	this->m_procID     = GetProcessId(this->m_procHandle);
	if (this->m_lock_heap == true) {
		bhl = HeapLock(this->m_heapHandle);
		if (FALSE == bhl) {
			PROCESS_HEAP_INFO_ERROR("ProcessHeapInfo::heap_enumerate -- !!! HeapLock: failed !!! ")
		}
	 bhw  = HeapWalk(this->m_heapHandle, &entry); 
	 if (FALSE == bhw) {
		 hwerr = GetLastError();
		 PRINT_ERROR_VALUE("HeapWalk terminated abnormally!!",hwerr)
		 if (hwerr != ERROR_NO_MORE_ITEMS) {
			 PROCESS_HEAP_INFO_ERROR("ProcessHeapInfo::enumerate: -- !! HeapWalk: failed!!")
			 }
		 }
	 if ((entry.wFlags & PROCESS_HEAP_REGION) != 0) {

				this->m_heap.m_lpData = entry.lpData;
				this->m_heap.m_cbData = entry.cbData;
				this->m_heap.m_cbOverhead = entry.cbOverhead;
				this->m_heap.m_iRegionIndex = entry.iRegionIndex;
				this->m_heap.m_Region.m_dwCommitedSize = entry.Region.dwCommittedSize;
				this->m_heap.m_Region.m_dwUncommitedSize = entry.Region.dwUnCommittedSize;
				this->m_heap.m_Region.m_lpFirstBlock = entry.Region.lpFirstBlock;
				this->m_heap.m_Region.m_lpLastBlock = entry.Region.lpLastBlock;
				this->m_heapList.push_back(this->m_heap);
				
			}
	  bhl = HeapUnlock(this->m_heapHandle);
	  if (FALSE == bhl) {
		    PROCESS_HEAP_INFO_ERROR("ProcessHeapInfo::heap_enumearte -- !!! HeapUnlock: failed !!!")
	        }
	}
	else {
		// Do not lock the heap.
		bhw = HeapWalk(this->m_heapHandle, &entry);
		if (FALSE == bhw) {
			hwerr = GetLastError();
			PRINT_ERROR_VALUE("HeapWalk terminated abnormally!!", hwerr)
			if (hwerr != ERROR_NO_MORE_ITEMS) {
				PROCESS_HEAP_INFO_ERROR("ProcessHeapInfo::enumerate: -- !! HeapWalk: failed!!")
			}
		}
		if ((entry.wFlags & PROCESS_HEAP_REGION) != 0) {

			this->m_heap.m_lpData = entry.lpData;
			this->m_heap.m_cbData = entry.cbData;
			this->m_heap.m_cbOverhead = entry.cbOverhead;
			this->m_heap.m_iRegionIndex = entry.iRegionIndex;
			this->m_heap.m_Region.m_dwCommitedSize = entry.Region.dwCommittedSize;
			this->m_heap.m_Region.m_dwUncommitedSize = entry.Region.dwUnCommittedSize;
			this->m_heap.m_Region.m_lpFirstBlock = entry.Region.lpFirstBlock;
			this->m_heap.m_Region.m_lpLastBlock = entry.Region.lpLastBlock;
			this->m_heapList.push_back(this->m_heap);
		}
	}
}	
	
	
	
	


