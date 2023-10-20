
#ifndef __GMS_MALLOC_H__
#define __GMS_MALLOC_H__



namespace  file_info {

    const unsigned int gGMS_MALLOC_MAJOR = 2U;
    const unsigned int gGMS_MALLOC_MINOR = 1U;
    const unsigned int gGMS_MALLOC_MICRO = 0U;
    const unsigned int gGMS_MALLOC_FULLVER =
      1000U*gGMS_MALLOC_MAJOR+100U*gGMS_MALLOC_MINOR+10U*gGMS_MALLOC_MICRO;
    const char * const pgGMS_MALLOC_CREATE_DATE = "01-10-2019 19:14 +00200 (TUE 01 OCT 2019 GMT+2)";
    const char * const pgGMS_MALLOC_BUILD_DATE  = __DATE__ ":" __TIME__;
    const char * const pgGMS_MALLOC_AUTHOR      =  "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";
    const char * const pgGMS_MALLOC_DESCRIPT    =  "Malloc wrappers.";
}

#include <cstdint>
#include <malloc.h>
#include <sys/mman.h>
#include <linux/mman.h>
#include <assert>
#include <omp.h> // OMP allocators
#include <alloca.h>
#include "GMS_error_macros.h"


//Handy macros for page alignment (mmap)
#if !defined(NORMAL_PAGE_4KiB)
    #define NORMAL_PAGE_4KiB (4* 1024)
#endif

#if !defined(ALIGN_TO_PAGE_4KiB)
     #define ALIGN_TO_PAGE_4KiB(x) \
        (((x) + NORMAL_PAGE_4KiB - 1) / NORMAL_PAGE_4KiB * NORMAL_PAGE_4KiB)
#endif

#if !defined(HUGE_PAGE_2MiB)
    #define HUGE_PAGE_2MiB (2 * 1024 * 1024)
#endif

#if !defined(ALIGN_TO_PAGE_2MiB)
     #define ALIGN_TO_PAGE_2MiB(x) \
        (((x) + HUGE_PAGE_2MiB - 1) / HUGE_PAGE_2MiB * HUGE_PAGE_2MiB)
#endif

#if !defined(HUGE_PAGE_1GiB)
    #define HUGE_PAGE_1GiB (1024 * 1024 * 1024)
#endif

#if !defined(ALIGN_TO_PAGE_1GiB)
    #define ALIGN_TO_PAGE_1GiB(x)
       (((x) + HUGE_PAGE_1GiB - 1) / HUGE_PAGE_1GiB * HUGE_PAGE_1GiB)
#endif


namespace gms {
	namespace common {


	  /*
               Removing Windows support!!
               Bumping version:  major to 2 (complete rewrite).
               Adding mmap allocations aligned to various page size
               boundaries.
               Modified: 11-10-2020 11:14AM +00200
               Modified: 30-08-2021 09:44AM +00200
            */

            __ATTR_ALWAYS_INLINE__
	    __ATTR_ALIGN__(32)
	    std::size_t
	    add_cacheline_pad_float(const std::size_t len) {
                constexpr std::size_t misalign = 0ULL;
		std::size_t pad = 0ULL;
		misalign = len % 16;
		pad = len + (misalign == 0 ? 0 : 16 - misalign);
		return (pad);
	    }


	    __ATTR_ALWAYS_INLINE__
	    __ATTR_ALIGN__(32)
	    std::size_t
	    add_cacheline_pad_double(const std::size_t len) {
                constexpr std::size_t misalign = 0ULL;
		std::size_t pad = 0ULL;
		misalign = len % 8;
		pad = len + (misalign == 0 ? 0 : 8 - misalign);
		return (pad);
   	    }
		
	
	      __ATTR_COLD__
	      __attribute__ ((malloc))
	      __attribute__ ((returns_nonnull))
	      __attribute__ ((assume_aligned(64)))
	      __attribute__ ((alloc_size(1)));
              void * gms_mm_malloc(const std::size_t len,
	                           const std::size_t alignment) {

		     void * __restrict ptr = NULL;
	             ptr = _mm_malloc(len, alignment));
                     if (NULL == ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                 std::cerr << " Not implemented yet!!";
#endif
	                 ABORT_ON_ERROR(__PRETTY_FUNCTION__, MALLOC_FAILED)
	            }
                     return (ptr);
	      }


	      __ATTR_COLD__
	      void gms_mm_free(void * __restrict ptr) {
                   _mm_free(ptr);
	      }


#define ALIGNED_ALLOCA(p, s, a) { \
  char *c_ptr = alloca((s) + (a)); \
  int64 off1 = (a) - 1; \
  int64 i_ptr64 = ((int64)c_ptr) + off1); \
  (p) = (void*)(i_ptr64 - (i_ptr64 & off1)); \
}

             


	        //
		// Using MMAP allocations for 4KiB page, 2MiB page and 1GiB page.
		// Subroutines include an error handling
		// See: https://man7.org/linux/man-pages/man2/mmap.2.html
		//
		// void *mmap(void *addr, size_t length, int prot, int flags,
                //  int fd, off_t offset);


#if !defined (MMAP_FUNC_4KIB_BODY)
#define MMAP_FUNC_4KIB_BODY  \   
            do{       \                        
                  void * ptr = nullptr;  \
                  int32_t n4KiB_pages = 0;     \
                  int32_t len = 0;             \
                  len = static_cast<int32_t>(length);  \
                  n4KiB_pages = len/4096;              \
                  if(len != n4KiB_pages*4096) {        \
                      n4KiB_pages++;                    \
                  }                                    \
                  len = n4KiB_pages*4096;                \
                  ptr = mmap(NULL,static_cast<std::size_t>(len),  \
                  prot,flag,fd,offset);  \
            }while(0)
#endif


#if !defined (MMAP_FUNC_2MIB_BODY)
#define MMAP_FUNC_2MIB_BODY     \
           do{         \
               void * ptr = nullptr;  \
               int32_t n2MiB_pages = 0;  \
               int32_t len = 0;  \
               len = static_cast<int32_t>(length);  \
               n2MiB_pages = len/2097152; \
               if(len != n2MiB_pages*2097152) { \
                  n2MiB_pages++; \
               } \
               len = n2MiB_pages*2097152; \
               ptr = mmap(NULL,static_cast<std::size_t>(len), \
               prot,flag,fd,offset); \
           }while(0)   
#endif


#if !defined (MMAP_FUNC_1GIB_BODY)
#define MMAP_FUNC_1GIB_BODY \
           do{ \
              void  * ptr = nullptr; \
              int32_t n1GiB_pages = 0; \
              int32_t len = 0;    \
              len = static_cast<int32_t>(length); \
              n1GiB_pages = len/1073741824; \
              if(len != n1GiB_pages*1073741824) { \
                 n1GiB_pages++;  \
              } \
              len = n1GiB_pages*1073741824;
              ptr = mmap(NULL,static_cast<std::size_t>(len), \
              prot,flag,fd,offset); \
           }while(0)
#endif

		 __ATTR_COLD__
		 __attribute__ ((malloc))
		 __attribute__ ((alloc_size(1)))
		 __attribute__ ((returns_nonnull))
		 __attribute__ ((assume_aligned(4096)))
		 void * gms_mmap_4KiB(const std::size_t len,
		                      const int32_t prot,
				      const int32_t flags,
				      const int32_t fd,
				      const off_t offset) {

                      MMAP_FUNC_4KIB_BODY(double)
                      if((ptr == (void*)(-1))) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                   std::cerr << "Requested stack-backtrace -- not implemented yet!!"
#endif
		           ABORT_ON_ERROR(__PRETTY_FUNCTION__, MALLOC_FAILED)
                        }
			return (ptr);
		 }
				          

                __ATTR_COLD__
		__attribute__ ((malloc))
		__attribute__ ((alloc_size(1)))
		__attribute__ ((returns_nonnull))
		__attribute__ ((assume_aligned(2097152)))
		void *  gms_mmap_2MiB(const std::size_t len,
		                      const int32_t prot,
				      const int32_t flags,
				      const int32_t fd,
				      const off_t offset) {
				     
				      
                      MMAP_FUNC_2MIB_BODY
		      if((ptr == (void*)(-1))) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                  std::cerr << "Requested stack-backtrace -- not implemented yet!!"
#endif
		          ABORT_ON_ERROR(__PRETTY_FUNCTION__, MALLOC_FAILED)
                       }
                       return (ptr);
		}


		__ATTR_COLD__
		__attribute__ ((malloc))
		__attribute__ ((alloc_size(1)))
		__attribute__ ((returns_nonnull))
		__attribute__ ((assume_aligned(1073741824)))
                void * gms_mmap_1GiB(const std::size_t len,
		                     const int32_t prot,
				     const int32_t flags,
				     const int32_t fd,
				     const off_t offset) {

		      MMAP_FUNC_1GIB_BODY
		      if((ptr == (void*)(-1))) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	                   std::cerr << "Requested stack-backtrace -- not implemented yet!!"
#endif
		           ABORT_ON_ERROR(__PRETTY_FUNCTION__, MALLOC_FAILED)
                       }
		       return (ptr);
		}

		
		__ATTR_COLD__		       
                int gms_ummap(void * __restrict ptr,
		              const std::size_t len) {
                     return (munmap(ptr,size));
		}

	

	}
}


#endif /*__GMS_MALLOC_H__*/
