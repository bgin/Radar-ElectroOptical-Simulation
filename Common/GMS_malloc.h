
#ifndef __GMS_MALLOC_H__
#define __GMS_MALLOC_H__



namespace  file_info {

    const unsigned int gGMS_MALLOC_MAJOR = 1U;
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
#include <complex>


class AVXVec8;
class AVX512Vec16;
class AVXc8f32;


namespace gms {
	namespace common {


	  /*
               Removing Windows support!!
               Bumping version:  minor to 1
               Adding mmap allocations aligned to various page size
               boundaries.
               Modified: 11-10-2020 11:14AM +00200
            */
		
	
	

	


			//
		//	Unaligned malloc wrapper
		//  Returns: double * 
		//  No error handling implemented
		//

	        double * gms_dmallocu( const std::size_t) __ATTR_COLD__
		                                          __attribute__ ((alloc_size(1)))
							  __attribute__ ((malloc))
							  __attribute__ ((returns_nonnull));

		//
		// Unaligned malloc wrapper
		// Returns: float * 
		// No error handling implemented
		//

	        float * gms_fmallocu( const std::size_t)  __ATTR_COLD__
		                                          __attribute__ ((alloc_size(1)))
							  __attribute__ ((malloc))
							  __attribute__ ((returns_nonnull));

		//
		// Unaligned malloc wrapper
		// Returns: int64_t * 
		// No error handling implemented
		//

	        int32_t * gms_imallocu( const std::size_t)  __ATTR_COLD__
		                                          __attribute__ ((alloc_size(1)))
							  __attribute__ ((malloc))
							  __attribute__ ((returns_nonnull));

		

		//
		//	Aligned malloc wrapper
		//  Returns: double * 
		//  No error handling implemented
		//

		double * gms_dmalloca( const std::size_t, const int32_t) __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: float *
		// No error handling implemented
		//

	        float * gms_fmalloca( const std::size_t, const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: int32_t * 
		// No error handling implemented
		//

	        int32_t * gms_imalloca( const std::size_t, const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: std::complex<float> *
		// No error handling implemented
		//

		std::complex<float> * gms_cmplxr4_malloca(const std::size_t, const int32_t)  __ATTR_COLD__
		                                                                             __attribute__ ((malloc))
									                     __attribute__ ((returns_nonnull))
									                     __attribute__ ((assume_aligned(64)))
									                     __attribute__ ((alloc_size(1)));
		//
		// Aligned malloc wrapper
		// Returns: AVXVec8 * 
		// No error handling implemented
		//
	        
		AVXVec8 * gms_avxvec8_malloca(const std::size_t,const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(32)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: AVX512Vec16 * 
		// No error handling implemented
		//

		AVX512Vec16 * gms_avx512vec16_malloca(const std::size_t,const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

                //
		// Aligned malloc wrapper
		// Returns: AVXc8f32 * 
		// No error handling implemented
		//
                
                AVXc8f32 * gms_avxc8f32_malloca(const std::size_t, const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(32)))
									 __attribute__ ((alloc_size(1)));

		//
		// Error handling wrappers
		//

		//
		//	Aligned malloc wrapper
		//  Returns: double *
		//  Error checking and handling (calls std::exit)
		//

	        double * gms_edmalloca(const std::size_t, int32_t)       __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: float *
		// Error checking and handling (calls std::exit)
		//

	        float * gms_efmalloca(const std::size_t, const int32_t)   __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: int32_t *
		// Error checking and handling (calls std::exit)
		//

	        int32_t * gms_eimalloca4(const std::size_t, const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: int64_t *
		// Error checking and handling (calls std::exit)
		//

	        int64_t * gms_eimalloca(const std::size_t, const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: AVXVec8 * 
		// Error checking and handling (calls std::exit)
		//

		std::complex<float> * gms_cmplxr4_emalloca(const std::size_t, const int32_t)  __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(64)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: AVXVec8 * 
		// Error checking and handling (calls std::exit)
		//

		AVXVec8 * gms_avxvec8_emalloca(const std::size_t,
					       const int32_t)             __ATTR_COLD__
		                                                          __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(32)))
									 __attribute__ ((alloc_size(1)));

		//
		// Aligned malloc wrapper
		// Returns: AVX512Vec16 * 
		// Error checking and handling (calls std::exit)
		//
		AVX512Vec16 * gms_avx512vec16_emalloca(const std::size_t,
						       const int32_t)    __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(32)))
									 __attribute__ ((alloc_size(1)));

                //
		// Aligned malloc wrapper
		// Returns: AVXc8f32 * 
		// Error checking and handling (calls std::exit)
		// 

                AVXc8f32 * gms_avxc8f32_emalloca(const std::size_t,
                                                 const int32_t)          __ATTR_COLD__
		                                                         __attribute__ ((malloc))
									 __attribute__ ((returns_nonnull))
									 __attribute__ ((assume_aligned(32)))
									 __attribute__ ((alloc_size(1)));

		//
		// Few special functions for padding of possibly unaligned rows of flat (multidimensional) arrays.
		//

	       
		float * gms_efmalloca_padded2D(const int32_t,
						const int32_t,
						const int32_t,
						int32_t &) __ATTR_COLD__
						           __attribute__ ((malloc))
							   __attribute__ ((alloc_size(1,2,3)))
							   __attribute__ ((returns_nonnull));

	
		double * gms_edmalloca_padded2D(const int32_t,
						const int32_t,
						const int32_t,
					        int32_t &) __ATTR_COLD__
						           __attribute__ ((malloc))
							   __attribute__ ((alloc_size(1,2,3)))
							   __attribute__ ((returns_nonnull));

		
		int32_t * gms_eimalloca4_padded2D(const int32_t,
						  const int32_t,
						  const int32_t,
						  int32_t &) __ATTR_COLD__
						           __attribute__ ((malloc))
							   __attribute__ ((alloc_size(1,2,3)))
							   __attribute__ ((returns_nonnull));


                //
		// Using MMAP allocations for 4KiB page, 2MiB page and 1GiB page.
		// Subroutines include an error handling
		// See: https://man7.org/linux/man-pages/man2/mmap.2.html
		//
		// void *mmap(void *addr, size_t length, int prot, int flags,
                //  int fd, off_t offset);
                
	        // Multiplicity of small page allocation (4KiB)

		double * gms_edmmap_4KiB(const std::size_t,
		                            const int32_t,
					    const int32_t,
					    const int32_t,
					    const off_t) __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(4096)));

		float * gms_efmmap_4KiB(const std::size_t,
		                           const int32_t,
					   const int32_t,
					   const int32_t,
					   const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(4096)));

		int32_t * gms_immap_4KiB(const std::size_t,
		                            const int32_t,
					    const int32_t,
					    const int32_t,
					    const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(4096)));

		std::complex<float> *
		           gms_ec4mmap_4KiB(const std::size_t,
		                               const int32_t,
					       const int32_t,
					       const int32_t,
					       const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(4096)));

		
		AVXc8f32 * gms_avxc8f32_mmap_4KiB(const std::size_t,
		                                     const int32_t,
					             const int32_t,
					             const int32_t,
					             const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(4096)));
#if defined __AVX512F__
		AVX512Vec16 * gms_avx512vec16_mmap_4KiB(const std::size_t,
		                                           const int32_t,
					                   const int32_t,
					                   const int32_t,
					                   const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(4096)));
#endif
		AVXVec8 *  gms_avxvec8_mmap_4KiB(const std::size_t,
		                                    const int32_t,
					            const int32_t,
					            const int32_t,
					            const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(4096)));

		// Multiplicity of medium page allocation (2MiB)

		double * gms_edmmap_2MiB(const std::size_t,
		                            const int32_t,
					    const int32_t,
					    const int32_t,
					    const off_t) __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(2097152)));

		float * gms_efmmap_2MiB(const std::size_t,
		                           const int32_t,
					   const int32_t,
					   const int32_t,
					   const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(2097152)));

		int32_t * gms_immap_2MiB(const std::size_t,
		                            const int32_t,
					    const int32_t,
					    const int32_t,
					    const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(2097152)));

		std::complex<float> *
		           gms_ec4mmap_2MiB(const std::size_t,
		                               const int32_t,
					       const int32_t,
					       const int32_t,
					       const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(2097152)));

		
		AVXc8f32 * gms_avxc8f32_mmap_2MiB(const std::size_t,
		                                     const int32_t,
					             const int32_t,
					             const int32_t,
					             const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(2097152)));
#if defined __AVX512F__
		AVX512Vec16 * gms_avx512vec16_mmap_2MiB(const std::size_t,
		                                           const int32_t,
					                   const int32_t,
					                   const int32_t,
					                   const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(2097152)));
#endif
		AVXVec8 *  gms_avxvec8_mmap_2MiB(const std::size_t,
		                                    const int32_t,
					            const int32_t,
					            const int32_t,
					            const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(2097152)));

	
                // Multiplicity of huge page allocation (1GiB)

		double * gms_edmmap_1GiB(const std::size_t,
		                            const int32_t,
					    const int32_t,
					    const int32_t,
					    const off_t) __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(1073741824)));

		float * gms_efmmap_1GiB(const std::size_t,
		                           const int32_t,
					   const int32_t,
					   const int32_t,
					   const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(1073741824)));

		int32_t * gms_immap_1GiB(const std::size_t,
		                            const int32_t,
					    const int32_t,
					    const int32_t,
					    const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(1073741824)));

		std::complex<float> *
		           gms_ec4mmap_1GiB(const std::size_t,
		                               const int32_t,
					       const int32_t,
					       const int32_t,
					       const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(1073741824)));

		
		AVXc8f32 * gms_avxc8f32_mmap_1GiB(const std::size_t,
		                                     const int32_t,
					             const int32_t,
					             const int32_t,
					             const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(1073741824)));
#if defined __AVX512F__
		AVX512Vec16 * gms_avx512vec16_mmap_1GiB(const std::size_t,
		                                           const int32_t,
					                   const int32_t,
					                   const int32_t,
					                   const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(1073741824)));
#endif
		AVXVec8 *  gms_avxvec8_mmap_1GiB(const std::size_t,
		                                    const int32_t,
					            const int32_t,
					            const int32_t,
					            const off_t)  __ATTR_COLD__
					                 __attribute__ ((malloc))
							 __attribute__ ((alloc_size(1)))
							 __attribute__ ((returns_nonnull))
							 __attribute__ ((assume_aligned(1073741824)));

		

	}
}


#endif /*__GMS_MALLOC_H__*/
