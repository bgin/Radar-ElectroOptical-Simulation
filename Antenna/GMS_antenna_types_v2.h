

#ifndef __GMS_ANTENNA_TYPES_V2_H__
#define __GMS_ANTENNA_TYPES_V2_H__ 201120221117


namespace file_info {

     const unsigned int GMS_ANTENNA_TYPES_V2_MAJOR = 1;
     const unsigned int GMS_ANTENNA_TYPES_V2_MINOR = 1;
     const unsigned int GMS_ANTENNA_TYPES_V2_MICRO = 0;
     const unsigned int GMS_ANTENNA_TYPES_V2_FULLVER =
       1000U*GMS_ANTENNA_TYPES_V2_MAJOR+100U*GMS_ANTENNA_TYPES_V2_MINOR+
       10U*GMS_ANTENNA_TYPES_V2_MICRO;
     const char * const GMS_ANTENNA_TYPES_V2_CREATION_DATE = "20-11-2022 11:17 +00200 (SUN 20 NOV 2022 GMT+2)";
     const char * const GMS_ANTENNA_TYPES_V2_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_ANTENNA_TYPES_V2_SYNOPSIS      = "Antenna common data types."

}


/*
 Purpose:
 !                        Derived data types for 'antenna_sensor' module implementation.
 !                        Various characteristics of different antenna types  
 !                        Based mainly on book titled (rus):          
 !                        Проектирование антенно фидерных устройств. Жук М.С. Молочков Ю.Б
*/


#include <cstdint>
#include <complex>
#include <vector>
#include <cstring> // std::memcpy
#include "GMS_config"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)
#define USE_GMS_ANTENNA_TYPES_V2_NT_STORES 0
#endif

namespace gms {


          namespace  radiolocation {
          

              struct __ATTR_ALIGN__(64) E_c4_t {
                      // Complex electric field.
                      std::complex<float> * __restrict ex
                      std::complex<float> * __restrict ey;
                      std::complex<float> * __restrict ez;
                      std::size_t                      npts;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline E_c4_t() noexcept(true) {
                          
                          this->npts = 0ULL;
                          this->ex  = NULL;
                          this->ey  = NULL;
                          this->ez  = NULL;
                      } 
                          
                     inline E_c4_t(const std::size_t pts) {
                          using namespace gms::common;
                          this->npts = pts;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline E_c4_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->ex = (std::complex<float>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ey = (std::complex<float>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ez = (std::complex<float>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->ex = (std::complex<float>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ey = (std::complex<float>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ez = (std::complex<float>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->ex = (std::complex<float>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ey = (std::complex<float>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ez = (std::complex<float>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline E_c4_t(const std::vector<std::complex<float>> &e_x,    //shall be of the same size (no error checking implemented)
                             const std::vector<std::complex<float>> &e_y,    //shall be of the same size (no error checking implemented)
                             const std::vector<std::complex<float>> &e_z) {  //shall be of the same size (no error checking implemented)
                          using namespace gms::common;
                          
                          this->npts = e_x.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t len = sizeof(std::complex<float>)*this->npts;
                          std::memcpy(this->ex,&e_x[0],len);
                          std::memcpy(this->ey,&e_y[0],len);
                          std::memcpy(this->ez,&e_z[0],len);       
                     }
                             
                             
                      
                     inline E_c4_t(const std::size_t pts,
                             const std::complex<float> * __restrict e_x,
                             const std::complex<float> * __restrict e_y,
                             const std::complex<float> * __restrict e_z) {
                          using namespace gms::common;
                          this->npts = npts;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->ex[0],&e_x[0],this->npts);
	                  avx512_uncached_memmove(&this->ey[0],&e_y[0],this->npts);
	                  avx512_uncached_memmove(&this->ez[0],&e_z[0],this->npts);
#else
	                  avx512_cached_memmove(&this->ex[0],&e_x[0],this->npts);
	                  avx512_cached_memmove(&this->ey[0],&e_y[0],this->npts);
	                  avx512_cached_memmove(&this->ez[0],&e_z[0],this->npts);
#endif
                   }  
                   
                    inline  E_c4_t(E_c4_t && rhs) {
                          
                          this->npts = rhs.npts;
                          this->ex   = &rhs.ex[0];
                          this->ey   = &rhs.ey[0];
                          this->ez   = &rhs.ez[0];
                          rhs.npts   = 0ULL;
                          rhs.ex     = NULL;
                          rhs.ey     = NULL;
                          rhs.ez     = NULL;
                      }
                                 
                      E_c4_t(const E_c4_t &)             = delete;
                      
                     inline ~E_c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap(this->ex,this->npts);
                             gms_unmap(this->ey,this->npts);
                             gms_unmap(this->ez,this->npts);
                          }
                          else {
                              gms_mm_free(this->ex);
                              gms_mm_free(this->ey);
                              gms_mm_free(this->ez);
                          }
                      }
                      
                      E_c4_t & operator=(const E_c4_t &) = delete;
                      
                     inline E_c4_t & operator=(E_c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->ex);
                           gms_mm_free(this->ey);
                           gms_mm_free(this->ez);
                           this->npts = rhs.npts;
                           this->ex   = &rhs.ex[0];
                           this->ey   = &rhs.ey[0];
                           this->ez   = &rhs.ez[0];
                           rhs.npts   = 0ULL;
                           rhs.ex     = NULL;
                           rhs.ey     = NULL;
                           rhs.ez     = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->ex  = (std::complex<float>*)
                                         gms_mm_malloc(this->npts,64ULL);
                          this->ey  = (std::complex<float>*)
                                         gms_mm_malloc(this->npts,64ULL);
                          this->ez  = (std::complex<float>*)
                                         gms_mm_malloc(this->npts,64ULL);
                      }
                    
              };

              
                struct __ATTR_ALIGN__(32) H_c4_t {
                     
                      std::complex<float> * __restrict mx
                      std::complex<float> * __restrict my;
                      std::complex<float> * __restrict mz;
                      std::size_t                      npts;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                      inline H_c4_t() noexcept(true) {
                          
                          this->npts = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                          this->mz  = NULL;
                      } 
                          
                      inline H_c4_t(const std::size_t pts) {
                         
                          this->npts = pts;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                      inline H_c4_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline  H_c4_t(const std::vector<std::complex<float>> &m_x,    //shall be of the same size (no error checking implemented)
                                    const std::vector<std::complex<float>> &m_y,    //shall be of the same size (no error checking implemented)
                                    const std::vector<std::complex<float>> &m_z) {  //shall be of the same size (no error checking implemented)
                                                   
                          this->npts = m_x.size();
                          allocate()
                          this->ismmap = false;
                          const std::size_t len = sizeof(std::complex<float>)*this->npts;
                          std::memcpy(this->mx,&m_x[0],len);
                          std::memcpy(this->my,&m_y[0],len);
                          std::memcpy(this->mz,&m_z[0],len);       
                     }
                      
                     inline H_c4_t(const std::size_t pts,
                                   const std::complex<float> * __restrict m_x,
                                   const std::complex<float> * __restrict m_y,
                                   const std::complex<float> * __restrict m_z) {
                         
                          this->npts = npts;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->npts);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->npts);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->npts);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->npts);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->npts);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->npts);
#endif
                   }  
                   
                    inline  H_c4_t(H_c4_t && rhs) {
                          
                          this->npts = rhs.npts;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          this->mz   = &rhs.mz[0];
                          rhs.npts   = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                          rhs.mz     = NULL;
                      }
                                 
                    inline  H_c4_t(const H_c4_t &)             = delete;
                      
                    inline  ~H_c4_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->mx,this->npts);
                             gms_unmap(this->my,this->npts);
                             gms_unmap(this->mz,this->npts);
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                              gms_mm_free(this->mz);
                          }
                      }
                      
                     inline H_c4_t & operator=(const H_c4_t &) = delete;
                      
                     inline  H_c4_t & operator=(H_c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->mx);
                           gms_mm_free(this->my);
                           gms_mm_free(this->mz);
                           this->npts = rhs.npts;
                           this->mx   = &rhs.mx[0];
                           this->my   = &rhs.my[0];
                           this->mz   = &rhs.mz[0];
                           rhs.npts   = 0ULL;
                           rhs.mx     = NULL;
                           rhs.my     = NULL;
                           rhs.mz     = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (std::complex<float>*)
                                         gms_mm_malloc(this->npts,64ULL);
                          this->my  = (std::complex<float>*)
                                         gms_mm_malloc(this->npts,64ULL);
                          this->mz  = (std::complex<float>*)
                                         gms_mm_malloc(this->npts,64ULL); 
                      }
                    
              };


                struct __ATTR_ALIGN__(32) E_c8_t {
                      // Complex electric field.
                      std::complex<double> * __restrict ex
                      std::complex<double> * __restrict ey;
                      std::complex<double> * __restrict ez;
                      std::size_t                      npts;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline E_c8_t() noexcept(true) {
                          
                          this->npts = 0ULL;
                          this->ex  = NULL;
                          this->ey  = NULL;
                          this->ez  = NULL;
                      } 
                          
                     inline E_c8_t(const std::size_t pts) {
                          this->npts = pts;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline E_c8_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->ex = (std::complex<double>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ey = (std::complex<double>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ez = (std::complex<double>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->ex = (std::complex<double>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ey = (std::complex<double>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ez = (std::complex<double>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->ex = (std::complex<double>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ey = (std::complex<double>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ez = (std::complex<double>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline E_c8_t(const std::vector<std::complex<double>> &e_x,    //shall be of the same size (no error checking implemented)
                                   const std::vector<std::complex<double>> &e_y,    //shall be of the same size (no error checking implemented)
                                   const std::vector<std::complex<double>> &e_z) {  //shall be of the same size (no error checking implemented)
                          
                          this->npts = e_x.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t len = sizeof(std::complex<float>)*this->npts;
                          std::memcpy(this->ex,&e_x[0],len);
                          std::memcpy(this->ey,&e_y[0],len);
                          std::memcpy(this->ez,&e_z[0],len);       
                     }
                           
                      
                     inline E_c8_t(const std::size_t pts,
                                   const std::complex<double> * __restrict e_x,
                                   const std::complex<double> * __restrict e_y,
                                   const std::complex<double> * __restrict e_z) {
                                   
                          this->npts = npts;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->ex[0],&e_x[0],this->npts);
	                  avx512_uncached_memmove(&this->ey[0],&e_y[0],this->npts);
	                  avx512_uncached_memmove(&this->ez[0],&e_z[0],this->npts);
#else
	                  avx512_cached_memmove(&this->ex[0],&e_x[0],this->npts);
	                  avx512_cached_memmove(&this->ey[0],&e_y[0],this->npts);
	                  avx512_cached_memmove(&this->ez[0],&e_z[0],this->npts);
#endif
                   }  
                   
                   inline   E_c8_t(E_c8_t && rhs) {
                          
                          this->npts = rhs.npts;
                          this->ex   = &rhs.ex[0];
                          this->ey   = &rhs.ey[0];
                          this->ez   = &rhs.ez[0];
                          rhs.npts   = 0ULL;
                          rhs.ex     = NULL;
                          rhs.ey     = NULL;
                          rhs.ez     = NULL;
                      }
                                 
                      E_c8_t(const E_c48t &)             = delete;
                      
                    inline ~E_c8_t() {
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->ex,this->npts);
                             gms_unmap(this->ey,this->npts);
                             gms_unmap(this->ez,this->npts);
                          }
                          else {
                             gms_mm_free(this->ex);
                             gms_mm_free(this->ey);
                             gms_mm_free(this->ez);
                          }
                      }
                      
                      E_c8_t & operator=(const E_c8_t &) = delete;
                      
                     inline E_c8_t & operator=(E_c8_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->ex);
                           gms_mm_free(this->ey);
                           gms_mm_free(this->ez);
                           this->npts = rhs.npts;
                           this->ex   = &rhs.ex[0];
                           this->ey   = &rhs.ey[0];
                           this->ez   = &rhs.ez[0];
                           rhs.npts   = 0ULL;
                           rhs.ex     = NULL;
                           rhs.ey     = NULL;
                           rhs.ez     = NULL;
                           return (*this);
                      }
                      
                      void allocate() {
                            using namespace gms::common;
                            this->ex  = (std::complex<double>*)
                                         gms_mm_malloc(this->npts,64ULL);
                            this->ey  = (std::complex<double>*)
                                         gms_mm_malloc(this->npts,64ULL);
                            this->ez  = (std::complex<double>*)
                                         gms_mm_malloc(this->npts,64ULL);
                      }
              };

              
               struct __ATTR_ALIGN__(32) H_c8_t {
                      // Complex magnetic field.
                      std::complex<double> * __restrict mx
                      std::complex<double> * __restrict my;
                      std::complex<double> * __restrict mz;
                      std::size_t                      npts;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline H_c8_t() noexcept(true) {
                          
                          this->npts = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                          this->mz  = NULL;
                      } 
                          
                     inline H_c8_t(const std::size_t pts) {
                     
                          this->npts = pts;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline H_c8_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->mx = (std::complex<double>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->my = (std::complex<double>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mz = (std::complex<double>*)
                                                 gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (std::complex<double>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->my = (std::complex<double>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mz = (std::complex<double>*)
                                                 gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (std::complex<double>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->my = (std::complex<double>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mz = (std::complex<double>*)
                                                 gms_mmap_1GiB(length,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   H_c8_t(const std::vector<std::complex<double>> &m_x,    //shall be of the same size (no error checking implemented)
                                    const std::vector<std::complex<double>> &m_y,    //shall be of the same size (no error checking implemented)
                                    const std::vector<std::complex<double>> &m_z) {  //shall be of the same size (no error checking implemented)
                         
                          this->npts = m_x.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t len = sizeof(std::complex<double>)*this->npts;
                          std::memcpy(this->mx,&m_x[0],len);
                          std::memcpy(this->my,&m_y[0],len);
                          std::memcpy(this->mz,&m_z[0],len);       
                     }
                      
                   inline   H_c8_t(const std::size_t pts,
                                   const std::complex<double> * __restrict m_x,
                                   const std::complex<double> * __restrict m_y,
                                   const std::complex<double> * __restrict m_z) {
                          
                          this->npts = npts;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->npts);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->npts);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->npts);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->npts);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->npts);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->npts);
#endif
                   }  
                   
                  inline  H_c8_t(H_c8_t && rhs) {
                          
                          this->npts = rhs.npts;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          this->mz   = &rhs.mz[0];
                          rhs.npts   = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                          rhs.mz     = NULL;
                      }
                                 
                      H_c8_t(const H_c8_t &)             = delete;
                      
                   inline   ~H_c8_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->mx,this->npts);
                             gms_unmap(this->my,this->npts);
                             gms_unmap(this->mz,this->npts);
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                              gms_mm_free(this->mz);
                          }
                      }
                      
                      H_c8_t & operator=(const H_c8_t &) = delete;
                      
                    inline  H_c8_t & operator=(H_c8_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->mx);
                           gms_mm_free(this->my);
                           gms_mm_free(this->mz);
                           this->npts = rhs.npts;
                           this->mx   = &rhs.mx[0];
                           this->my   = &rhs.my[0];
                           this->mz   = &rhs.mz[0];
                           rhs.npts   = 0ULL;
                           rhs.mx     = NULL;
                           rhs.my     = NULL;
                           rhs.mz     = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() {
                        using namespace gms::common;
                        this->mx  = (std::complex<double>*)
                                         gms_mm_malloc(this->npts,64ULL);
                        this->my  = (std::complex<double>*)
                                         gms_mm_malloc(this->npts,64ULL);
                        this->mz  = (std::complex<double>*)
                                         gms_mm_malloc(this->npts,64ULL);
                    }
                    
              };
   

              struct __ATTR_ALIGN__(64) E_r4_t {
                      //  ! Complex Electric  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                      float * __restrict exr;
                      float * __restrict exi;
                      float * __restrict eyr;
                      float * __restrict eyi;
                      float * __restrict ezr;
                      float * __restrict ezi;
                      std::size_t       npts;
                      bool              ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,7)
#endif
                      inline E_r4_t() {
                         this->npts = 0ULL;
                         this->exr  = NULL;
                         this->exi  = NULL;
                         this->eyr  = NULL;
                         this->eyi  = NULL;
                         this->ezr  = NULL;
                         this->ezi  = NULL;
                      }                    
                      
                      inline E_r4_t(const std::size_t pts) {
                             this->npts = pts;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline E_r4_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->exr = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->exi = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->eyr = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->eyi = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ezr = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ezi = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->exr = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->exi = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->eyr = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->eyi = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ezr = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ezi = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->exr = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->exi = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->eyr = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->eyi = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ezr = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ezi = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline E_r4_t(const std::vector<float> &e_xr,
                                    const std::vector<float> &e_xi,
                                    const std::vector<float> &e_yr,
                                    const std::vector<float> &e_yi,
                                    const std::vector<float> &e_zr,
                                    const std::vector<float> &e_zi) {
                               
                               this->npts = e_xr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t len = sizeof(float)*this->npts;
                               std::memcpy(this->exr,&e_xr[0],len);
                               std::memcpy(this->exi,&e_xi[0],len);
                               std::memcpy(this->eyr,&e_yr[0],len);
                               std::memcpy(this->eyi,&e_yi[0],len);
                               std::memcpy(this->ezr,&e_zr[0],len);
                               std::memcpy(this->ezi,&e_zi[0],len);     
                      }
                      
                      inline E_r4_t(const std::size_t pts,
                                    const float * __restrict e_xr,   
                                    const float * __restrict e_xi,
                                    const float * __restrict e_yr,
                                    const float * __restrict e_yi,
                                    const float * __restrict e_zr,
                                    const float * __restrict e_zi) {
                         
                          this->npts = pts;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->exr[0],&e_xr[0],this->npts);
	                  avx512_uncached_memmove(&this->exi[0],&e_xi[0],this->npts);
	                  avx512_uncached_memmove(&this->eyr[0],&e_yr[0],this->npts);
	                  avx512_uncached_memmove(&this->eyi[0],&e_yi[0],this->npts);
	                  avx512_uncached_memmove(&this->ezr[0],&e_zr[0],this->npts);
	                  avx512_uncached_memmove(&this->ezi[0],&e_zi[0],this->npts);
#else
	                  avx512_cached_memmove(&this->exr[0],&e_xr[0],this->npts);
	                  avx512_cached_memmove(&this->exi[0],&e_xi[0],this->npts);
	                  avx512_cached_memmove(&this->eyr[0],&e_yr[0],this->npts);
	                  avx512_cached_memmove(&this->eyi[0],&e_yi[0],this->npts);
	                  avx512_cached_memmove(&this->ezr[0],&e_zr[0],this->npts);
	                  avx512_cached_memmove(&this->ezi[0],&e_zi[0],this->npts);
#endif          
                      }
                      
                      inline E_r4_t(E_r4_t &&rhs) {
                           
                          this->npts = rhs.npts;
                          this->exr  = &rhs.exr[0];
                          this->exi  = &rhs.exi[0];
                          this->eyr  = &rhs.eyr[0];
                          this->eyi  = &rhs.eyi[0];
                          this->ezr  = &rhs.ezr[0];
                          this->ezi  = &rhs.ezi[0];
                          rhs.npts   = 0ULL;
                          rhs.exr    = NULL;
                          rhs.exi    = NULL;
                          rhs.eyr    = NULL;
                          rhs.eyi    = NULL;
                          rhs.ezr    = NULL;
                          rhs.ezi    = NULL;
                      }   
                        
                      E_r4_t(const E_r4_t &)             = delete;
                      
                      inline ~E_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->exr,this->npts);
                              gms_unmap(this->exi,this->npts);
                              gms_unmap(this->eyr,this->npts); 
                              gms_unmap(this->eyi,this->npts);
                              gms_unmap(this->ezr,this->npts);
                              gms_unmap(this->ezi,this->npts);
                           }
                           else {
                               gms_mm_free(this->exr);
                               gms_mm_free(this->exi);
                               gms_mm_free(this->eyr);
                               gms_mm_free(this->eyi);
                               gms_mm_free(this->ezr);
                               gms_mm_free(this->ezi);
                           }
                      }
                      
                      E_r4_t & operator=(const E_r4_t &) = delete;
                      
                      inline E_r4_t & operator=(E_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->exr);
                            gms_mm_free(this->exi);
                            gms_mm_free(this->eyr);
                            gms_mm_free(this->eyi);
                            gms_mm_free(this->ezr);
                            gms_mm_free(this->ezi);
                            this->npts = rhs.npts;
                            this->exr  = &rhs.exr[0];
                            this->exi  = &rhs.exi[0];
                            this->eyr  = &rhs.eyr[0];
                            this->eyi  = &rhs.eyi[0];
                            this->ezr  = &rhs.ezr[0];
                            this->ezi  = &rhs.ezi[0];
                            rhs.npts   = 0ULL;
                            rhs.exr    = NULL;
                            rhs.exi    = NULL;
                            rhs.eyr    = NULL;
                            rhs.eyi    = NULL;
                            rhs.ezr    = NULL;
                            rhs.ezi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->exr  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->exi  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->eyr  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->eyi  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->ezr  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->ezi  = (float*)gms_mm_malloc(this->npts,64ULL);
                   }    
                   
              };


              struct __ATTR_ALIGN__(64) H_r4_t {
                      //  ! Complex Magnetic  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                      float * __restrict     mxr;
                      float * __restrict     mxi;
                      float * __restrict     myr;
                      float * __restrict     myi;
                      float * __restrict     mzr;
                      float * __restrict     mzi;
                      std::size_t            npts;
                      bool                   ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,7)
#endif
                     
                     inline H_r4_t() {
                         this->npts = 0ULL;
                         this->mxr  = NULL;
                         this->mxi  = NULL;
                         this->myr  = NULL;
                         this->myi  = NULL;
                         this->mzr  = NULL;
                         this->mzi  = NULL;
                      }                    
                      
                      inline H_r4_t(const std::size_t pts) {
                             this->npts = pts;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline H_r4_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->mxr = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mxr = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mxr = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline H_r4_t(const std::vector<float> &mxr,
                                    const std::vector<float> &mxi,
                                    const std::vector<float> &myr,
                                    const std::vector<float> &myi,
                                    const std::vector<float> &mzr,
                                    const std::vector<float> &mzi) {
                               
                               this->npts = mxr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t len = sizeof(float)*this->npts;
                               std::memcpy(this->mxr,&mxr[0],len);
                               std::memcpy(this->mxi,&mxi[0],len);
                               std::memcpy(this->myr,&myr[0],len);
                               std::memcpy(this->myi,&myi[0],len);
                               std::memcpy(this->mzr,&mzr[0],len);
                               std::memcpy(this->mzi,&mzi[0],len);     
                      }
                      
                      inline H_r4_t(const std::size_t pts,
                                    const float * __restrict mxr,   
                                    const float * __restrict mxi,
                                    const float * __restrict myr,
                                    const float * __restrict myi,
                                    const float * __restrict mzr,
                                    const float * __restrict mzi) {
                         
                          this->npts = pts;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&mxr[0],this->npts);
	                  avx512_uncached_memmove(&this->mxi[0],&mxi[0],this->npts);
	                  avx512_uncached_memmove(&this->myr[0],&myr[0],this->npts);
	                  avx512_uncached_memmove(&this->myi[0],&myi[0],this->npts);
	                  avx512_uncached_memmove(&this->mzr[0],&mzr[0],this->npts);
	                  avx512_uncached_memmove(&this->mzi[0],&mzi[0],this->npts);
#else
	                  avx512_cached_memmove(&this->mxr[0],&mxr[0],this->npts);
	                  avx512_cached_memmove(&this->mxi[0],&mxi[0],this->npts);
	                  avx512_cached_memmove(&this->myr[0],&myr[0],this->npts);
	                  avx512_cached_memmove(&this->myi[0],&myi[0],this->npts);
	                  avx512_cached_memmove(&this->mzr[0],&mzr[0],this->npts);
	                  avx512_cached_memmove(&this->mzi[0],&mzi[0],this->npts);
#endif          
                      }
                      
                      inline H_r4_t(H_r4_t &&rhs) {
                           
                          this->npts = rhs.npts;
                          this->mxr  = &rhs.mxr[0];
                          this->mxi  = &rhs.mxi[0];
                          this->myr  = &rhs.myr[0];
                          this->myi  = &rhs.myi[0];
                          this->mzr  = &rhs.mzr[0];
                          this->mzi  = &rhs.mzi[0];
                          rhs.npts   = 0ULL;
                          rhs.mxr    = NULL;
                          rhs.mxi    = NULL;
                          rhs.myr    = NULL;
                          rhs.myi    = NULL;
                          rhs.mzr    = NULL;
                          rhs.mzi    = NULL;
                      }   
                        
                      H_r4_t(const H_r4_t &)    = delete;
                      
                      inline ~H_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->mxr,this->npts);
                              gms_unmap(this->mxi,this->npts);
                              gms_unmap(this->myr,this->npts); 
                              gms_unmap(this->myi,this->npts);
                              gms_unmap(this->mzr,this->npts);
                              gms_unmap(this->mzi,this->npts);
                           }
                           else {
                               gms_mm_free(this->mxr);
                               gms_mm_free(this->mxi);
                               gms_mm_free(this->myr);
                               gms_mm_free(this->myi);
                               gms_mm_free(this->mzr);
                               gms_mm_free(this->mzi);
                           }
                      }
                      
                      H_r4_t & operator=(const H_r4_t &) = delete;
                      
                      inline H_r4_t & operator=(H_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->mxr);
                            gms_mm_free(this->mxi);
                            gms_mm_free(this->myr);
                            gms_mm_free(this->myi);
                            gms_mm_free(this->mzr);
                            gms_mm_free(this->mzi);
                            this->npts = rhs.npts;
                            this->mxr  = &rhs.mxr[0];
                            this->mxi  = &rhs.mxi[0];
                            this->myr  = &rhs.myr[0];
                            this->myi  = &rhs.myi[0];
                            this->mzr  = &rhs.mzr[0];
                            this->mzi  = &rhs.mzi[0];
                            rhs.npts   = 0ULL;
                            rhs.mxr    = NULL;
                            rhs.mxi    = NULL;
                            rhs.myr    = NULL;
                            rhs.myi    = NULL;
                            rhs.mzr    = NULL;
                            rhs.mzi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->mxr  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->mxi  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->myr  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->myi  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->mzr  = (float*)gms_mm_malloc(this->npts,64ULL);
                        this->mzi  = (float*)gms_mm_malloc(this->npts,64ULL);
                   }    
                   
                      
              };


              struct __ATTR_ALIGN__(64) E_r8_t {
                      //  ! Complex Electric  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                     
                      double * __restrict      exr;
                      double * __restrict      exi;
                      double * __restrict      eyr;
                      double * __restrict      eyi;
                      double * __restrict      ezr;
                      double * __restrict      ezi;
                      std::size_t              npts;
                      bool                     ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,7)
#endif

                     inline E_r8_t() {
                         this->npts = 0ULL;
                         this->exr  = NULL;
                         this->exi  = NULL;
                         this->eyr  = NULL;
                         this->eyi  = NULL;
                         this->ezr  = NULL;
                         this->ezi  = NULL;
                      }                    
                      
                     inline E_r8_t(const std::size_t pts) {
                             this->npts = pts;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                     inline E_r8_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->exr = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->exi = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->eyr = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->eyi = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ezr = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ezi = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->exr = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->exi = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->eyr = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->eyi = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ezr = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ezi = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->exr = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->exi = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->eyr = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->eyi = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ezr = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ezi = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline E_r8_t(const std::vector<double> &e_xr,
                                    const std::vector<double> &e_xi,
                                    const std::vector<double> &e_yr,
                                    const std::vector<double> &e_yi,
                                    const std::vector<double> &e_zr,
                                    const std::vector<double> &e_zi) {
                               
                               this->npts = e_xr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t len = sizeof(double)*this->npts;
                               std::memcpy(this->exr,&e_xr[0],len);
                               std::memcpy(this->exi,&e_xi[0],len);
                               std::memcpy(this->eyr,&e_yr[0],len);
                               std::memcpy(this->eyi,&e_yi[0],len);
                               std::memcpy(this->ezr,&e_zr[0],len);
                               std::memcpy(this->ezi,&e_zi[0],len);     
                      }
                      
                      inline E_r8_t(const std::size_t pts,
                                    const double * __restrict e_xr,   
                                    const double * __restrict e_xi,
                                    const double * __restrict e_yr,
                                    const double * __restrict e_yi,
                                    const double * __restrict e_zr,
                                    const double * __restrict e_zi) {
                         
                          this->npts = pts;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->exr[0],&e_xr[0],this->npts);
	                  avx512_uncached_memmove(&this->exi[0],&e_xi[0],this->npts);
	                  avx512_uncached_memmove(&this->eyr[0],&e_yr[0],this->npts);
	                  avx512_uncached_memmove(&this->eyi[0],&e_yi[0],this->npts);
	                  avx512_uncached_memmove(&this->ezr[0],&e_zr[0],this->npts);
	                  avx512_uncached_memmove(&this->ezi[0],&e_zi[0],this->npts);
#else
	                  avx512_cached_memmove(&this->exr[0],&e_xr[0],this->npts);
	                  avx512_cached_memmove(&this->exi[0],&e_xi[0],this->npts);
	                  avx512_cached_memmove(&this->eyr[0],&e_yr[0],this->npts);
	                  avx512_cached_memmove(&this->eyi[0],&e_yi[0],this->npts);
	                  avx512_cached_memmove(&this->ezr[0],&e_zr[0],this->npts);
	                  avx512_cached_memmove(&this->ezi[0],&e_zi[0],this->npts);
#endif          
                      }
                      
                      inline E_r8_t(E_r8_t &&rhs) {
                           
                          this->npts = rhs.npts;
                          this->exr  = &rhs.exr[0];
                          this->exi  = &rhs.exi[0];
                          this->eyr  = &rhs.eyr[0];
                          this->eyi  = &rhs.eyi[0];
                          this->ezr  = &rhs.ezr[0];
                          this->ezi  = &rhs.ezi[0];
                          rhs.npts   = 0ULL;
                          rhs.exr    = NULL;
                          rhs.exi    = NULL;
                          rhs.eyr    = NULL;
                          rhs.eyi    = NULL;
                          rhs.ezr    = NULL;
                          rhs.ezi    = NULL;
                      }   
                        
                      E_r8_t(const E_r8_t &)             = delete;
                      
                      inline ~E_r8_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->exr,this->npts);
                              gms_unmap(this->exi,this->npts);
                              gms_unmap(this->eyr,this->npts); 
                              gms_unmap(this->eyi,this->npts);
                              gms_unmap(this->ezr,this->npts);
                              gms_unmap(this->ezi,this->npts);
                           }
                           else {
                               gms_mm_free(this->exr);
                               gms_mm_free(this->exi);
                               gms_mm_free(this->eyr);
                               gms_mm_free(this->eyi);
                               gms_mm_free(this->ezr);
                               gms_mm_free(this->ezi);
                           }
                      }
                      
                      E_r8_t & operator=(const E_r8_t &) = delete;
                      
                      inline E_r8_t & operator=(E_r8_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->exr);
                            gms_mm_free(this->exi);
                            gms_mm_free(this->eyr);
                            gms_mm_free(this->eyi);
                            gms_mm_free(this->ezr);
                            gms_mm_free(this->ezi);
                            this->npts = rhs.npts;
                            this->exr  = &rhs.exr[0];
                            this->exi  = &rhs.exi[0];
                            this->eyr  = &rhs.eyr[0];
                            this->eyi  = &rhs.eyi[0];
                            this->ezr  = &rhs.ezr[0];
                            this->ezi  = &rhs.ezi[0];
                            rhs.npts   = 0ULL;
                            rhs.exr    = NULL;
                            rhs.exi    = NULL;
                            rhs.eyr    = NULL;
                            rhs.eyi    = NULL;
                            rhs.ezr    = NULL;
                            rhs.ezi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->exr  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->exi  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->eyr  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->eyi  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->ezr  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->ezi  = (double*)gms_mm_malloc(this->npts,64ULL);
                   }    
                   
              };


           struct __ATTR_ALIGN__(64) H_r8_t {
                      //  ! Complex Magnetic  field  decomposed into real and imaginary parts 
                      //  ! To be used mainly by the integrators.
                      double * __restrict mxr;
                      double * __restrict mxi;
                      double * __restrict myr;
                      double * __restrict myi;
                      double * __restrict mzr;
                      double * __restrict mzi;
                      std::size_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,7)
#endif

                        inline H_r8_t() {
                         this->npts = 0ULL;
                         this->mxr  = NULL;
                         this->mxi  = NULL;
                         this->myr  = NULL;
                         this->myi  = NULL;
                         this->mzr  = NULL;
                         this->mzi  = NULL;
                      }                    
                      
                      inline H_r8_t(const std::size_t pts) {
                             this->npts = pts;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline H_r8_t(const std::size_t length,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->npts = length;
                             switch (fsize) {
                                 case:0
                                      this->mxr = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mxi = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->myr = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->myi = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mzr = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->mzi = (double*)
                                                  gms_mmap_4KiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mxr = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mxi = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->myr = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->myi = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mzr = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->mzi = (double*)
                                                  gms_mmap_2MiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mxr = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mxi = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->myr = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->myi = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mzr = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->mzi = (double*)
                                                  gms_mmap_1GiB(length,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline H_r8_t(const std::vector<double> &mxr,
                                    const std::vector<double> &mxi,
                                    const std::vector<double> &myr,
                                    const std::vector<double> &myi,
                                    const std::vector<double> &mzr,
                                    const std::vector<double> &mzi) {
                               
                               this->npts = mxr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t len = sizeof(double)*this->npts;
                               std::memcpy(this->mxr,&mxr[0],len);
                               std::memcpy(this->mxi,&mxi[0],len);
                               std::memcpy(this->myr,&myr[0],len);
                               std::memcpy(this->myi,&myi[0],len);
                               std::memcpy(this->mzr,&mzr[0],len);
                               std::memcpy(this->mzi,&mzi[0],len);     
                      }
                      
                      inline H_r8_t(const std::size_t pts,
                                    const double * __restrict mxr,   
                                    const double * __restrict mxi,
                                    const double * __restrict myr,
                                    const double * __restrict myi,
                                    const double * __restrict mzr,
                                    const double * __restrict mzi) {
                         
                          this->npts = pts;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_TYPES_V2_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&mxr[0],this->npts);
	                  avx512_uncached_memmove(&this->mxi[0],&mxi[0],this->npts);
	                  avx512_uncached_memmove(&this->myr[0],&myr[0],this->npts);
	                  avx512_uncached_memmove(&this->myi[0],&myi[0],this->npts);
	                  avx512_uncached_memmove(&this->mzr[0],&mzr[0],this->npts);
	                  avx512_uncached_memmove(&this->mzi[0],&mzi[0],this->npts);
#else
	                  avx512_cached_memmove(&this->mxr[0],&mxr[0],this->npts);
	                  avx512_cached_memmove(&this->mxi[0],&mxi[0],this->npts);
	                  avx512_cached_memmove(&this->myr[0],&myr[0],this->npts);
	                  avx512_cached_memmove(&this->myi[0],&myi[0],this->npts);
	                  avx512_cached_memmove(&this->mzr[0],&mzr[0],this->npts);
	                  avx512_cached_memmove(&this->mzi[0],&mzi[0],this->npts);
#endif          
                      }
                      
                      inline H_r8_t(H_r8_t &&rhs) {
                           
                          this->npts = rhs.npts;
                          this->mxr  = &rhs.mxr[0];
                          this->mxi  = &rhs.mxi[0];
                          this->myr  = &rhs.myr[0];
                          this->myi  = &rhs.myi[0];
                          this->mzr  = &rhs.mzr[0];
                          this->mzi  = &rhs.mzi[0];
                          rhs.npts   = 0ULL;
                          rhs.mxr    = NULL;
                          rhs.mxi    = NULL;
                          rhs.myr    = NULL;
                          rhs.myi    = NULL;
                          rhs.mzr    = NULL;
                          rhs.mzi    = NULL;
                      }   
                        
                      H_r8_t(const H_r8_t &)    = delete;
                      
                      inline ~H_r8_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->mxr,this->npts);
                              gms_unmap(this->mxi,this->npts);
                              gms_unmap(this->myr,this->npts); 
                              gms_unmap(this->myi,this->npts);
                              gms_unmap(this->mzr,this->npts);
                              gms_unmap(this->mzi,this->npts);
                           }
                           else {
                               gms_mm_free(this->mxr);
                               gms_mm_free(this->mxi);
                               gms_mm_free(this->myr);
                               gms_mm_free(this->myi);
                               gms_mm_free(this->mzr);
                               gms_mm_free(this->mzi);
                           }
                      }
                      
                      H_r8_t & operator=(const H_r8_t &) = delete;
                      
                      inline H_r8_t & operator=(H_r8_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->mxr);
                            gms_mm_free(this->mxi);
                            gms_mm_free(this->myr);
                            gms_mm_free(this->myi);
                            gms_mm_free(this->mzr);
                            gms_mm_free(this->mzi);
                            this->npts = rhs.npts;
                            this->mxr  = &rhs.mxr[0];
                            this->mxi  = &rhs.mxi[0];
                            this->myr  = &rhs.myr[0];
                            this->myi  = &rhs.myi[0];
                            this->mzr  = &rhs.mzr[0];
                            this->mzi  = &rhs.mzi[0];
                            rhs.npts   = 0ULL;
                            rhs.mxr    = NULL;
                            rhs.mxi    = NULL;
                            rhs.myr    = NULL;
                            rhs.myi    = NULL;
                            rhs.mzr    = NULL;
                            rhs.mzi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->mxr  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->mxi  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->myr  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->myi  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->mzr  = (double*)gms_mm_malloc(this->npts,64ULL);
                        this->mzi  = (double*)gms_mm_malloc(this->npts,64ULL);
                   }    
                   
              };


              typedef struct __ATTR_ALIGN__(32) JE_c4_t {
                      // Complex magnetic current
                      JE_c4_t()                            = delete;
                      JE_c4_t(const JE_c4_t &)             = delete;
                      JE_c4_t & operator=(const JE_c4_t &) = delete;
                      std::complex<float> * __restrict je_x;
                      std::complex<float> * __restrict je_y;
                      std::complex<float> * __restrict je_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JE_c4_t;


              typedef struct __ATTR_ALIGN__(32) JM_c4_t {
                      // Complex electric current
                      JM_c4_t()                            = delete;
                      JM_c4_t(const JM_c4_t &)             = delete;
                      JM_c4_t & operator=(const JM_c4_t &) = delete;
                      std::complex<float> * __restrict jm_x;
                      std::complex<float> * __restrict jm_y;
                      std::complex<float> * __restrict jm_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JM_c4_t;


              typedef struct __ATTR_ALIGN__(32) JE_c8_t {
                      // Complex magnetic current
                      JE_c8_t()                            = delete;
                      JE_c8_t(const JE_c8_t &)             = delete;
                      JE_c8_t & operator=(const JE_c8_t &) = delete;
                      std::complex<double> * __restrict je_x;
                      std::complex<double> * __restrict je_y;
                      std::complex<double> * __restrict je_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JE_c8_t;


              typedef struct __ATTR_ALIGN__(32) JM_c8_t {
                      // Complex electric current
                      JM_c8_t()                            = delete;
                      JM_c8_t(const JM_c8_t &)             = delete;
                      JM_c8_t & operator=(const JM_c8_t &) = delete;
                      std::complex<double> * __restrict jm_x;
                      std::complex<double> * __restrict jm_y;
                      std::complex<double> * __restrict jm_z;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif
              } JM_c8_t;


              typedef struct __ATTR_ALIGN__(64) JE_r4_t {
                     // ! Complex Electric Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      JE_r4_t()                            = delete;
                      JE_r4_t(const JE_r4_t &)             = delete;
                      JE_r4_t & operator=(const JE_r4_t &) = delete;
                      float * __restrict je_xr;
                      float * __restrict je_xi;
                      float * __restrict je_yr;
                      float * __restrict je_yi;
                      float * __restrict je_zr;
                      float * __restrict je_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JE_r4_t;


              typedef struct __ATTR_ALIGN__(64) JM_r4_t {
                     // ! Complex Magnetic Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      JM_r4_t()                            = delete;
                      JM_r4_t(const JM_r4_t &)             = delete;
                      JM_r4_t & operator=(const JM_r4_t &) = delete;
                      float * __restrict jm_xr;
                      float * __restrict jm_xi;
                      float * __restrict jm_yr;
                      float * __restrict jm_yi;
                      float * __restrict jm_zr;
                      float * __restrict jm_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JM_r4_t;


              typedef struct __ATTR_ALIGN__(64) JE_r8_t {
                     // ! Complex Electric Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      JE_r8_t()                            = delete;
                      JE_r8_t(const JE_r8_t &)             = delete;
                      JE_r8_t & operator=(const JE_r8_t &) = delete;
                      double * __restrict je_xr;
                      double * __restrict je_xi;
                      double * __restrict je_yr;
                      double * __restrict je_yi;
                      double * __restrict je_zr;
                      double * __restrict je_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JE_r8_t;


              typedef struct __ATTR_ALIGN__(64) JM_r8_t {
                     // ! Complex Magnetic Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      JM_r8_t()                            = delete;
                      JM_r8_t(const JM_r8_t &)             = delete;
                      JM_r8_t & operator=(const JM_r8_t &) = delete;
                      double * __restrict jm_xr;
                      double * __restrict jm_xi;
                      double * __restrict jm_yr;
                      double * __restrict jm_yi;
                      double * __restrict jm_zr;
                      double * __restrict jm_zi;
                      int32_t            npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
              } JM_r8_t;


              typedef struct __ATTR_ALIGN__(32) eikr_c4_t {
                      // Time-Harmonic complex exponential
                      eikr_c4_t()                              = delete;
                      eikr_c4_t(const eikr_c4_t &)             = delete;
                      eikr_c4_t & operator=(const eikr_c4_t &) = delete;
                      float               * __restrict R;
                      std::complex<float> * __restrict ce;
                      float                            k;
                      int32_t                          npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif              
              } eikr_c4_t;


               typedef struct __ATTR_ALIGN__(32) eikr_c8_t {
                      // Time-Harmonic complex exponential
                      eikr_c8_t()                              = delete;
                      eikr_c8_t(const eikr_c8_t &)             = delete;
                      eikr_c8_t & operator=(const eikr_c8_t &) = delete;
                      double               * __restrict R;
                      std::complex<double> * __restrict ce;
                      double                            k;
                      int32_t                           npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif              
              } eikr_c8_t;


              typedef struct __ATTR_ALIGN__(32) eikr_r4_t {
                      // ! Time-Harmonic complex exponential decomposed into 
                      // ! real and imaginary parts
                      eikr_r4_t()                              = delete;
                      eikr_r4_t(const eikr_r4_t &)             = delete;
                      eikr_r4_t & operator=(const eikr_r4_t &) = delete;
                      float                 * __restrict R;
                      float                 * __restrict e_re;
                      float                 * __restrict e_im;
                      float                              k;
                      int32_t                            npts;
              } eikr_r4_t;


              typedef struct __ATTR_ALIGN__(32) eikr_r8_t {
                      // ! Time-Harmonic complex exponential decomposed into 
                      // ! real and imaginary parts
                      eikr_r8_t()                              = delete;
                      eikr_r8_t(const eikr_r8_t &)             = delete;
                      eikr_r8_t & operator=(const eikr_r8_t &) = delete;
                      double                 * __restrict R;
                      double                 * __restrict e_re;
                      double                 * __restrict e_im;
                      double                              k;
                      int32_t                             npts;
              } eikr_r8_t;


              // ! Formula (1-37)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f137_r4_t {
                       f137_r4_t()                             = delete;
                       f137_r4_t(const f137_r4_t &)            = delete;
                       f137_r4_t & operator=(const f137_r4_t &)= delete;
                       float                 * __restrict sinth;
                       float                 * __restrict F;
                       float                              ith;
                       float                              iph;
                       float                              ifac;
                       float                              omega;
                       float                              avsl;
                       int32_t                            nth;
                       int32_t                            nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif  
              } f137_r4_t;


              // ! Formula (1-37)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f137_r8_t {
                       f137_r8_t()                             = delete;
                       f137_r8_t(const f137_r8_t &)            = delete;
                       f137_r8_t & operator=(const f137_r8_t &)= delete;
                       double                 * __restrict sinth;
                       double                 * __restrict F;
                       double                              ith;
                       double                              iph;
                       double                              ifac;
                       double                              omega;
                       double                              avsl;
                       int32_t                             nth;
                       int32_t                             nph;
              } f137_r8_t;


              // ! Formula (1-38)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f138_r4_t {
                       f138_r4_t()                             = delete;
                       f138_r4_t(const f138_r4_t &)            = delete;
                       f138_r4_t & operator=(const f138_r4_t &)= delete;
                       float                 * __restrict sinth;
                       float                 * __restrict Fsqr;
                       float                              ith;
                       float                              iph;
                       float                              ifac;
                       float                              omega;
                       float                              avsl;
                       int32_t                            nth;
                       int32_t                            nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif  
              } f138_r4_t;


              // ! Formula (1-38)
              //! Average level of side lobes
              typedef struct __ATTR_ALIGN__(64) f138_r8_t {
                       f138_r8_t()                             = delete;
                       f138_r8_t(const f138_r8_t &)            = delete;
                       f138_r8_t & operator=(const f138_r8_t &)= delete;
                       double                 * __restrict sinth;
                       double                 * __restrict Fsqr;
                       double                              ith;
                       double                              iph;
                       double                              ifac;
                       double                              omega;
                       double                              avsl;
                       int32_t                             nth;
                       int32_t                             nph;
              } f138_r8_t;


              //! Formula (1-39)
              // ! Dispersion coefficient
              typedef struct __ATTR_ALIGN__(64) f139_r4_t {
                       f139_r4_t()                             = delete;
                       f139_r4_t(const f139_r4_t &)            = delete;
                       f139_r4_t & operator=(const f139_r4_t &)= delete; 
                       float                  * __restrict P;
                       float                  * __restrict sinth;
                       float                               ith;
                       float                               iph;
                       float                               omega;
                       int32_t                             nth;
                       int32_t                             nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,26)
#endif 
              } f139_r4_t;


              //! Formula (1-39)
              // ! Dispersion coefficient
              typedef struct __ATTR_ALIGN__(64) f139_r8_t {
                       f139_r8_t()                             = delete;
                       f139_r8_t(const f139_r8_t &)            = delete;
                       f139_r8_t & operator=(const f139_r8_t &)= delete;
                       double                  * __restrict P;
                       double                  * __restrict sinth;
                       double                               ith;
                       double                               iph;
                       double                               omega;
                       int32_t                              nth;
                       int32_t                              nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,16)
#endif 
              } f139_r8_t; 


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_c4_t {
                       Hve_c4_t()                              = delete();
                       Hve_c4_t(const Hve_c4_t &)              = delete();
                       Hve_c4_t & operator=(const Hve_c4_t &)  = delete();
                       struct JE_c4_t                        jec4;
                       struct eikr_c4_t                      ec4;
                       std::complex<float>                   ifac;
                       std::complex<float>      * __restrict he_x;
                       std::complex<float>      * __restrict he_y;
                       std::complex<float>      * __restrict he_z;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,60)
#endif
             } Hve_c4_t;


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_c8_t {
                       Hve_c8_t()                              = delete();
                       Hve_c8_t(const Hve_c8_t &)              = delete();
                       Hve_c8_t & operator=(const Hve_c8_t &)  = delete();
                       struct JE_c8_t                         jec8;
                       struct eikr_c8_t                       ec8;
                       std::complex<double>                   ifac;
                       std::complex<double>      * __restrict he_x;
                       std::complex<double>      * __restrict he_y;
                       std::complex<double>      * __restrict he_z;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,52)
#endif
             } Hve_c8_t;


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_r4_t {
                       Hve_r4_t()                              = delete();
                       Hve_r4_t(const Hve_r4_t &)              = delete();
                       Hve_r4_t & operator=(const Hve_r4_t &)  = delete();
                       struct JE_r4_t                        jer4;
                       struct eikr_r4_t                      er4;
                       float                                 if_re;
                       float                                 if_im;
                       float                    * __restrict he_xr;
                       float                    * __restrict he_xi;
                       float                    * __restrict he_yr;
                       float                    * __restrict he_yi;
                       float                    * __restrict he_zr;
                       float                    * __restrict he_zi;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Hve_r4_t;


             //  ! Formula (2-13)
             //  ! Hertz vector electric
             typedef struct __ATTR_ALIGN__(64) Hve_r8_t {
                       Hve_r8_t()                              = delete();
                       Hve_r8_t(const Hve_r8_t &)              = delete();
                       Hve_r8_t & operator=(const Hve_r8_t &)  = delete();
                       struct JE_r8_t                         jer8;
                       struct eikr_r8_t                       er8;
                       double                                 if_re;
                       double                                 if_im;
                       double                    * __restrict he_xr;
                       double                    * __restrict he_xi;
                       double                    * __restrict he_yr;
                       double                    * __restrict he_yi;
                       double                    * __restrict he_zr;
                       double                    * __restrict he_zi;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,30)
#endif
             } Hve_r8_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hvm_c4_t {
                       Hvm_c4_t()                              = delete();
                       Hvm_c4_t(const Hvm_c4_t &)              = delete();
                       Hvm_c4_t & operator=(const Hvm_c4_t &)  = delete();
                       struct JM_c4_t                        jmc4;
                       struct eikr_c4_t                      ec4;
                       std::complex<float>                   ifac;
                       std::complex<float>      * __restrict hm_x;
                       std::complex<float>      * __restrict hm_y;
                       std::complex<float>      * __restrict hm_z;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,60)
#endif
             } Hvm_c4_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hvm_c8_t {
                       Hvm_c8_t()                              = delete();
                       Hvm_c8_t(const Hvm_c8_t &)              = delete();
                       Hvm_c8_t & operator=(const Hvm_c8_t &)  = delete();
                       struct JM_c8_t                         jmc8;
                       struct eikr_c8_t                       ec8;
                       std::complex<double>                   ifac;
                       std::complex<double>      * __restrict hm_x;
                       std::complex<double>      * __restrict hm_y;
                       std::complex<double>      * __restrict hm_z;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,52)
#endif
             } Hvm_c8_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hve_r4_t {
                       Hvm_r4_t()                              = delete();
                       Hvm_r4_t(const Hvm_r4_t &)              = delete();
                       Hvm_r4_t & operator=(const Hvm_r4_t &)  = delete(); 
                       struct JM_r4_t                        jmr4;
                       struct eikr_r4_t                      er4;
                       float                                 if_re;
                       float                                 if_im;
                       float                    * __restrict hm_xr;
                       float                    * __restrict hm_xi;
                       float                    * __restrict hm_yr;
                       float                    * __restrict hm_yi;
                       float                    * __restrict hm_zr;
                       float                    * __restrict hm_zi;
                       int32_t                               npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Hvm_r4_t;


             //  ! Formula (2-15)
             //  ! Hertz vector magnetic
             typedef struct __ATTR_ALIGN__(64) Hvm_r8_t {
                       Hvm_r8_t()                              = delete();
                       Hvm_r8_t(const Hvm_r8_t &)              = delete();
                       Hvm_r8_t & operator=(const Hvm_r8_t &)  = delete();
                       struct JM_r8_t                         jmr8;
                       struct eikr_r8_t                       er8;
                       double                                 if_re;
                       double                                 if_im;
                       double                    * __restrict hm_xr;
                       double                    * __restrict hm_xi;
                       double                    * __restrict hm_yr;
                       double                    * __restrict hm_yi;
                       double                    * __restrict hm_zr;
                       double                    * __restrict hm_zi;
                       int32_t                                npts;  
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,30)
#endif
             } Hvm_r8_t;


             //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_r4_t {
                       Nev_r4_t()                              = delete();
                       Nev_r4_t(const Nev_r4_t &)              = delete();
                       Nev_r4_t & operator=(const Nev_r4_t &)  = delete();
                       struct JE_r4_t                         jer4;
                       struct eikr_r4_t                       er4;
                       float                     * __restrict costh;
                       float                     * __restrict ne_xr;
                       float                     * __restrict ne_xi;
                       float                     * __restrict ne_yr;
                       float                     * __restrict ne_yi;
                       float                     * __restrict ne_zr;
                       float                     * __restrict ne_zi;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_r4_t;

              //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_r8_t {
                       Nev_r8_t()                              = delete();
                       Nev_r8_t(const Nev_r8_t &)              = delete();
                       Nev_r8_t & operator=(const Nev_r8_t &)  = delete(); 
                       struct JE_r8_t                          jer8;
                       struct eikr_r8_t                         er8;
                       double                     * __restrict costh;
                       double                     * __restrict ne_xr;
                       double                     * __restrict ne_xi;
                       double                     * __restrict ne_yr;
                       double                     * __restrict ne_yi;
                       double                     * __restrict ne_zr;
                       double                     * __restrict ne_zi;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_r8_t;


             //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_c4_t {
                       Nev_c4_t()                              = delete();
                       Nev_c4_t(const Nev_c4_t &)              = delete();
                       Nev_c4_t & operator=(const Nev_c4_t &)  = delete();
                       struct JE_c4_t                         jec4;
                       struct eikr_c4_t                       ec4;
                       float                     * __restrict costh;
                       std::complex<float>       * __restrict ne_x;
                       std::complex<float>       * __restrict ne_y;
                       std::complex<float>       * __restrict ne_z;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_c4_t;


              //! Formula (2-22,2-23)
             typedef struct __ATTR_ALIGN__(64) Nev_c8_t {
                       Nev_c8_t()                              = delete();
                       Nev_c8_t(const Nev_c8_t &)              = delete();
                       Nev_c8_t & operator=(const Nev_c8_t &)  = delete();
                       struct JE_c8_t                          jcr8;
                       struct eikr_c8_t                         ec8;
                       double                     * __restrict costh;
                       std::complex<double>       * __restrict ne_x;
                       std::complex<double>       * __restrict ne_y;
                       std::complex<double>       * __restrict ne_z;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nev_c8_t;


             //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_r4_t {
                       Nvm_r4_t()                              = delete();
                       Nvm_r4_t(const Nvm_r4_t &)              = delete();
                       Nvm_r4_t & operator=(const Nvm_r4_t &)  = delete();
                       struct JM_r4_t                         jmr4;
                       struct eikr_r4_t                       er4;
                       float                     * __restrict costh;
                       float                     * __restrict nm_xr;
                       float                     * __restrict nm_xi;
                       float                     * __restrict nm_yr;
                       float                     * __restrict nm_yi;
                       float                     * __restrict nm_zr;
                       float                     * __restrict nm_zi;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_r4_t;

              //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_r8_t {
                       Nvm_r8_t()                              = delete();
                       Nvm_r8_t(const Nvm_r8_t &)              = delete();
                       Nvm_r8_t & operator=(const Nvm_r8_t &)  = delete();
                       struct JM_r8_t                          jmr8;
                       struct eikr_r8_t                         er8;
                       double                     * __restrict costh;
                       double                     * __restrict nm_xr;
                       double                     * __restrict nm_xi;
                       double                     * __restrict nm_yr;
                       double                     * __restrict nm_yi;
                       double                     * __restrict nm_zr;
                       double                     * __restrict nm_zi;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_r8_t;


             //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_c4_t {
                       Nvm_c4_t()                              = delete();
                       Nvm_c4_t(const Nvm_c4_t &)              = delete();
                       Nvm_c4_t & operator=(const Nvm_c4_t &)  = delete();
                       struct JM_c4_t                         jmc4;
                       struct eikr_c4_t                       ec4;
                       float                     * __restrict costh;
                       std::complex<float>       * __restrict nm_x;
                       std::complex<float>       * __restrict nm_y;
                       std::complex<float>       * __restrict nm_z;
                       int32_t                                npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_c4_t;


              //! Formula (2-24,2-25)
             typedef struct __ATTR_ALIGN__(64) Nmv_c8_t {
                       Nvm_c8_t()                              = delete();
                       Nvm_c8_t(const Nvm_c8_t &)              = delete();
                       Nvm_c8_t & operator=(const Nvm_c8_t &)  = delete();
                       struct JM_c8_t                          jmc8;
                       struct eikr_c8_t                         ec8;
                       double                     * __restrict costh;
                       std::complex<double>       * __restrict nm_x;
                       std::complex<double>       * __restrict nm_y;
                       std::complex<double>       * __restrict nm_z;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,38)
#endif
             } Nmv_c8_t;



            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_r4_t {
                       FFhev_r4_t()                            = delete;
                       FFhev_r4_t(const FFHev_r4_t &)          = delete;
                       FFhev_r4_t & operator=(const FFhev_r4_t &)= delete;
                       struct Nev_r4_t                         ner4;
                       struct eikr_r4_t                        er4;   // 96-bytes
                       float                      * __restrict hv_xr;
                       float                      * __restrict hv_xi;
                       float                      * __restrict hv_yr;
                       float                      * __restrict hv_yi;
                       float                      * __restrict hv_zr;
                       float                      * __restrict hv_zi;
                       float                                   if_re;
                       float                                   if_im; // 54-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,42)
#endif           
            } FFhev_r4_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_r8_t {
                       FFhev_r8_t()                            = delete;
                       FFhev_r8_t(const FFHev_r8_t &)          = delete;
                       FFhev_r8_t & operator=(const FFhev_r8_t &)= delete;
                       struct Nev_r8_t                          ner8;
                       struct eikr_r8_t                         er8;   // 96-bytes
                       double                      * __restrict hv_xr;
                       double                      * __restrict hv_xi;
                       double                      * __restrict hv_yr;
                       double                      * __restrict hv_yi;
                       double                      * __restrict hv_zr;
                       double                      * __restrict hv_zi;
                       double                                   if_re;
                       double                                   if_im; // 54-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,34)
#endif           
            } FFhev_r8_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_c4_t {
                       FFhev_c4_t()                            = delete;
                       FFhev_c4_t(const FFHev_c4_t &)          = delete;
                       FFhev_c4_t & operator=(const FFhev_c4_t &)= delete;
                       struct Nev_c4_t                         nec4;
                       struct eikr_c4_t                        ec4;   // 96-bytes
                       std::complex<float>        * __restrict hv_x;
                       std::complex<float>        * __restrict hv_y;
                       std::complex<float>        * __restrict hv_z;
                       std::complex<float>                     ifac; // 32-bytes
         
            } FFhev_r4_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhev_c8_t {
                       FFhev_c8_t()                            = delete;
                       FFhev_c8_t(const FFHev_c8_t &)          = delete;
                       FFhev_c8_t & operator=(const FFhev_c8_t &)= delete;
                       struct Nev_c8_t                          nec8;
                       struct eikr_c8_t                         ec8;   // 96-bytes
                       std::complex<double>        * __restrict hv_x;
                       std::complex<double>        * __restrict hv_y;
                       std::complex<double>        * __restrict hv_z;
                       std::complex<double>                     ifac;
                       
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,56)
#endif           
            } FFhev_c8_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector magnetic
            typedef struct __ATTR_ALIGN__(64) FFhmv_r4_t {
                       FFhmv_r4_t()                            = delete;
                       FFhmv_r4_t(const FFhmv_r4_t &)          = delete;
                       FFhmv_r4_t & operator=(const FFhmv_r4_t &)= delete;
                       struct Nmv_r4_t                         nmr4;
                       struct eikr_r4_t                        er4;   // 96-bytes
                       float                      * __restrict hm_xr;
                       float                      * __restrict hm_xi;
                       float                      * __restrict hm_yr;
                       float                      * __restrict hm_yi;
                       float                      * __restrict hm_zr;
                       float                      * __restrict hm_zi;
                       float                                   if_re;
                       float                                   if_im; // 54-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,42)
#endif           
            } FFhmv_r4_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhmv_r8_t {
                       FFhmv_r8_t()                            = delete;
                       FFhmv_r8_t(const FFHmv_r8_t &)          = delete;
                       FFhmv_r8_t & operator=(const FFhmv_r8_t &)= delete;
                       struct Nmv_r8_t                          nmr8;
                       struct eikr_r8_t                         er8;   // 96-bytes
                       double                      * __restrict hm_xr;
                       double                      * __restrict hm_xi;
                       double                      * __restrict hm_yr;
                       double                      * __restrict hm_yi;
                       double                      * __restrict hm_zr;
                       double                      * __restrict hm_zi;
                       double                                   if_re;
                       double                                   if_im; // 54-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,34)
#endif           
            } FFhmv_r8_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhmv_c4_t {
                       FFhmv_c4_t()                            = delete;
                       FFhmv_c4_t(const FFhmv_c4_t &)          = delete;
                       FFhmv_c4_t & operator=(const FFhmv_r4_t &)= delete;
                       struct Nmv_c4_t                         nec4;
                       struct eikr_c4_t                        ec4;   // 96-bytes
                       std::complex<float>        * __restrict hm_x;
                       std::complex<float>        * __restrict hm_y;
                       std::complex<float>        * __restrict hm_z;
                       std::complex<float>                     ifac; // 32-bytes
         
            } FFhmv_r4_t;


            //! Formula (2-22)
            //! "Far-field" Hertz vector electric
            typedef struct __ATTR_ALIGN__(64) FFhmv_c8_t {
                       FFhmv_c8_t()                            = delete;
                       FFhmv_c8_t(const FFhmv_c8_t &)          = delete;
                       FFhmv_c8_t & operator=(const FFhmv_c8_t &)= delete;
                       struct Nmv_c8_t                          nmc8;
                       struct eikr_c8_t                         ec8;   // 96-bytes
                       std::complex<double>        * __restrict hm_x;
                       std::complex<double>        * __restrict hm_y;
                       std::complex<double>        * __restrict hm_z;
                       std::complex<double>                     ifac;
                       
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,56)
#endif           
            } FFhmv_c8_t;


           //! Formula (2-26)
           //! Cosine of integration angle.
           typedef struct __ATTR_ALIGN__(64) f226_r4_t {
                       f226_r4_t()                             = delete;
                       f226_r4_t(const f226_r4_t &)            = delete;
                       f226_r4_t & operator=(const f226_r4_t) &  = delete;
                       float                      * __restrict cth;
                       float                      * __restrict cthi;
                       float                      * __restrict sth;
                       float                      * __restrict sthi;
                       float                      * __restrict cphphi;
                       float                      * __restrict ciang;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif   
           } f226_r4_t;


           //! Formula (2-26)
           //! Cosine of integration angle.
           typedef struct __ATTR_ALIGN__(64) f226_r8_t {
                       f226_r8_t()                             = delete;
                       f226_r8_t(const f226_r8_t &)            = delete;
                       f226_r8_t & operator=(const f226_r8_t) &  = delete;
                       double                      * __restrict cth;
                       double                      * __restrict cthi;
                       double                      * __restrict sth;
                       double                      * __restrict sthi;
                       double                      * __restrict cphphi;
                       double                      * __restrict ciang;
                       int32_t                                 npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif   
           } f226_r8_t;


           /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nesph_r4_t {
                       Nesph_r4_t()                               = delete;
                       Nesph_r4_t(const Nesph_r4_t &)             = delete;
                       Nesph_r4_t & operator=(const Nesph_r4_t &) = delete;
                       struct Nev_r4_t                          nev; //64-bytes
                       float                      * __restrict  cphi;
                       float                      * __restrict  cth;
                       float                      * __restrict  sphi;
                       float                      * __restrict  sth;
                       float                      * __restrict  nth_re;
                       float                      * __restrict  nth_im;
                       float                      * __restrict  nph_re;
                       float                      * __restrict  nph_im;
                       int32_t                                  npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif 
            } Nesph_r4_t;


              /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nesph_r8_t {
                       Nesph_r8_t()                               = delete;
                       Nesph_r8_t(const Nesph_r8_t &)             = delete;
                       Nesph_r8_t & operator=(const Nesph_r8_t &) = delete;
                       struct Nev_r8_t                           nev; //64-bytes
                       double                      * __restrict  cphi;
                       double                      * __restrict  cth;
                       double                      * __restrict  sphi;
                       double                      * __restrict  sth;
                       double                      * __restrict  nth_re;
                       double                      * __restrict  nth_im;
                       double                      * __restrict  nph_re;
                       double                      * __restrict  nph_im;
                       int32_t                                   npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif 
            } Nesph_r8_t;


              /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nesph_c4_t {
                       Nesph_c4_t()                               = delete;
                       Nesph_c4_t(const Nesph_c4_t &)             = delete;
                       Nesph_c4_t & operator=(const Nesph_c4_t &) = delete;
                       struct Nev_c4_t                          nev; //64-bytes
                       float                      * __restrict  cphi;
                       float                      * __restrict  cth;
                       float                      * __restrict  sphi;
                       float                      * __restrict  sth;
                       std::complex<float>        * __restrict  nth;
                       std::complex<float>        * __restrict  nph;
                       int32_t                                  npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } Nesph_c4_t;


              /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nesph_c8_t {
                       Nesph_c8_t()                               = delete;
                       Nesph_c8_t(const Nesph_c8_t &)             = delete;
                       Nesph_c8_t & operator=(const Nesph_c8_t &) = delete;
                       struct Nev_c8_t                           nev; //64-bytes
                       double                      * __restrict  cphi;
                       double                      * __restrict  cth;
                       double                      * __restrict  sphi;
                       double                      * __restrict  sth;
                       std::complex<double>        * __restrict  nth;
                       std::complex<double>        * __restrict  nph;
                       int32_t                                   npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } Nesph_c8_t;


             /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nmsph_r4_t {
                       Nmsph_r4_t()                               = delete;
                       Nmsph_r4_t(const Nmsph_r4_t &)             = delete;
                       Nmsph_r4_t & operator=(const Nmsph_r4_t &) = delete;
                       struct Nmv_r4_t                          nmv; //64-bytes
                       float                      * __restrict  cphi;
                       float                      * __restrict  cth;
                       float                      * __restrict  sphi;
                       float                      * __restrict  sth;
                       float                      * __restrict  nth_re;
                       float                      * __restrict  nth_im;
                       float                      * __restrict  nph_re;
                       float                      * __restrict  nph_im;
                       int32_t                                  npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif 
            } Nmsph_r4_t;


              /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nmsph_r8_t {
                       Nmsph_r8_t()                               = delete;
                       Nmsph_r8_t(const Nmsph_r8_t &)             = delete;
                       Nmsph_r8_t & operator=(const Nmsph_r8_t &) = delete;
                       struct Nmv_r8_t                           nmv; //64-bytes
                       double                      * __restrict  cphi;
                       double                      * __restrict  cth;
                       double                      * __restrict  sphi;
                       double                      * __restrict  sth;
                       double                      * __restrict  nth_re;
                       double                      * __restrict  nth_im;
                       double                      * __restrict  nph_re;
                       double                      * __restrict  nph_im;
                       int32_t                                   npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif 
            } Nmsph_r8_t;


              /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nmsph_c4_t {
                       Nmsph_c4_t()                               = delete;
                       Nmsph_c4_t(const Nmsph_c4_t &)             = delete;
                       Nmsph_c4_t & operator=(const Nmsph_c4_t &) = delete; 
                       struct Nmv_c4_t                          nmv; //64-bytes
                       float                      * __restrict  cphi;
                       float                      * __restrict  cth;
                       float                      * __restrict  sphi;
                       float                      * __restrict  sth;
                       std::complex<float>        * __restrict  nth;
                       std::complex<float>        * __restrict  nph;
                       int32_t                                  npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } Nmsph_c4_t;


              /*
               ! Formula (2-36)
               ! 'N'-electric and 'N'-magnetic vector functions of fields voltage
               ! in spherical coordinate system.
            */
            typedef struct __ATTR_ALIGN__(64) Nmsph_c8_t {
                       Nmsph_c8_t()                               = delete;
                       Nmsph_c8_t(const Nmsph_c8_t &)             = delete;
                       Nmsph_c8_t & operator=(const Nmsph_c8_t &) = delete;
                       struct Nmv_c8_t                           nmv; //64-bytes
                       double                      * __restrict  cphi;
                       double                      * __restrict  cth;
                       double                      * __restrict  sphi;
                       double                      * __restrict  sth;
                       std::complex<double>        * __restrict  nth;
                       std::complex<double>        * __restrict  nph;
                       int32_t                                   npts; //60-bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } Nmsph_c8_t;


            //! Formula (2-27)
            //! Far field (zone) 'electric current' for R>=2*D^2/gamma
            typedef struct __ATTR_ALIGN__(64) FFec_r4_t {
                        FFec_r4_t()                                = delete;
                        FFec_r4_t(const FFec_r4_t &)               = delete;
                        FFec_r4_t & operator=(const FFec_r4_t &)   = delete;
                        struct eikr_r4_t                          eikr;
                        struct Nesph_r4_t           * __restrict  nesp;
                        struct Nev_r4_t             * __restrict  nev;
                        float                       * __restrict  e_r;
                        float                       * __restrict  e_th;
                        float                       * __restrict  e_ph;
                        float                       * __restrict  ffec_xr;
                        float                       * __restrict  ffec_xi;
                        float                       * __restrict  ffec_yr;
                        float                       * __restrict  ffec_yi;
                        float                       * __restrict  ffec_zr;
                        float                       * __restrict  ffec_zi;
                        float                                     if_re;
                        float                                     if_im;
                        int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,60)
#endif 
            } FFec_r4_t;


            //! Formula (2-27)
            //! Far field (zone) 'electric current' for R>=2*D^2/gamma
            typedef struct __ATTR_ALIGN__(64) FFec_r8_t {
                        FFec_r8_t()                                = delete;
                        FFec_r8_t(const FFec_r8_t &)               = delete;
                        FFec_r8_t & operator=(const FFec_r8_t &)   = delete;
                        struct eikr_r8_t                           eikr;
                        struct Nesph_r8_t            * __restrict  nesp;
                        struct Nev_r8_t              * __restrict  nev;
                        double                       * __restrict  e_r;
                        double                       * __restrict  e_th;
                        double                       * __restrict  e_ph;
                        double                       * __restrict  ffec_xr;
                        double                       * __restrict  ffec_xi;
                        double                       * __restrict  ffec_yr;
                        double                       * __restrict  ffec_yi;
                        double                       * __restrict  ffec_zr;
                        double                       * __restrict  ffec_zi;
                        double                                     if_re;
                        double                                     if_im;
                        int32_t                                    npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,52)
#endif 
            } FFec_r8_t;


            //! Formula (2-27)
            //! Far field (zone) 'electric current' for R>=2*D^2/gamma
            typedef struct __ATTR_ALIGN__(64) FFec_c4_t {
                        FFec_c4_t()                                = delete;
                        FFec_c4_t(const FFec_c4_t &)               = delete;
                        FFec_c4_t & operator=(const FFec_c4_t &)   = delete;
                        struct eikr_c4_t                          eikr;
                        struct Nesph_c4_t           * __restrict  nesp;
                        struct Nev_c4_t             * __restrict  nev;
                        float                       * __restrict  e_r;
                        float                       * __restrict  e_th;
                        float                       * __restrict  e_ph;
                        std::complex<float>         * __restrict  ffec_x;
                        std::complex<float>         * __restrict  ffec_y;
                        std::complex<float>         * __restrict  ffec_z;
                        std::complex<float>                       ifac;
                        int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif 
            } FFec_c4_t;


            //! Formula (2-27)
            //! Far field (zone) 'electric current' for R>=2*D^2/gamma
            typedef struct __ATTR_ALIGN__(64) FFec_c8_t {
                        FFec_c8_t()                                = delete;
                        FFec_c8_t(const FFec_c8_t &)               = delete;
                        FFec_c8_t & operator=(const FFec_c8_t &)   = delete;
                        struct eikr_c8_t                           eikr;
                        struct Nesph_c8_t            * __restrict  nesp;
                        struct Nev_c8_t              * __restrict  nev;
                        double                       * __restrict  e_r;
                        double                       * __restrict  e_th;
                        double                       * __restrict  e_ph;
                        std::complex<double>         * __restrict  ffec_x;
                        std::complex<double>         * __restrict  ffec_y;
                        std::complex<double>         * __restrict  ffec_z;
                        std::complex<double>                       ifac;
                        int32_t                                    npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } FFec_c8_t;


            //    ! Formula (2-27)
            // ! Far field (zone) 'magnetic current' for R>=2*D^2/gamma
              typedef struct __ATTR_ALIGN__(64) FFmc_r4_t {
                        FFmc_r4_t()                                = delete;
                        FFmc_r4_t(const FFmc_r4_t &)               = delete;
                        FFmc_r4_t & operator=(const FFmc_r4_t &)   = delete;
                        struct eikr_r4_t                          eikr;
                        struct Nmsph_r4_t           * __restrict  nmsp;
                        struct Nmv_r4_t             * __restrict  nmv;
                        float                       * __restrict  e_r;
                        float                       * __restrict  e_th;
                        float                       * __restrict  e_ph;
                        float                       * __restrict  ffmc_xr;
                        float                       * __restrict  ffmc_xi;
                        float                       * __restrict  ffmc_yr;
                        float                       * __restrict  ffmc_yi;
                        float                       * __restrict  ffmc_zr;
                        float                       * __restrict  ffmc_zi;
                        float                                     if_re;
                        float                                     if_im;
                        int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,60)
#endif 
            } FFmc_r4_t;


            //! Formula (2-27)
            //! Far field (zone) 'magnetic current' for R>=2*D^2/gamma
            typedef struct __ATTR_ALIGN__(64) FFmc_r8_t {
                        FFmc_r8_t()                                = delete;
                        FFmc_r8_t(const FFmc_r8_t &)               = delete;
                        FFmc_r8_t & operator=(const FFmc_r8_t &)   = delete;
                        struct eikr_r8_t                           eikr;
                        struct Nmsph_r8_t            * __restrict  nmsp;
                        struct Nmv_r8_t              * __restrict  nmv;
                        double                       * __restrict  e_r;
                        double                       * __restrict  e_th;
                        double                       * __restrict  e_ph;
                        double                       * __restrict  ffmc_xr;
                        double                       * __restrict  ffmc_xi;
                        double                       * __restrict  ffmc_yr;
                        double                       * __restrict  ffmc_yi;
                        double                       * __restrict  ffmc_zr;
                        double                       * __restrict  ffmc_zi;
                        double                                     if_re;
                        double                                     if_im;
                        int32_t                                    npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,52)
#endif 
            } FFmc_r8_t;


            //! Formula (2-27)
            //! Far field (zone) 'electric current' for R>=2*D^2/gamma
            typedef struct __ATTR_ALIGN__(64) FFmc_c4_t {
                        FFmc_c4_t()                                = delete;
                        FFmc_c4_t(const FFmc_c4_t &)               = delete;
                        FFmc_c4_t & operator=(const FFmc_c4_t &)   = delete;
                        struct eikr_c4_t                          eikr;
                        struct Nmsph_c4_t           * __restrict  nmsp;
                        struct Nmv_c4_t             * __restrict  nmv;
                        float                       * __restrict  e_r;
                        float                       * __restrict  e_th;
                        float                       * __restrict  e_ph;
                        std::complex<float>         * __restrict  ffmc_x;
                        std::complex<float>         * __restrict  ffmc_y;
                        std::complex<float>         * __restrict  ffmc_z;
                        std::complex<float>                       ifac;
                        int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif 
            } FFmc_c4_t;


            //! Formula (2-27)
            //! Far field (zone) 'electric current' for R>=2*D^2/gamma
            typedef struct __ATTR_ALIGN__(64) FFmc_c8_t {
                        FFmc_c8_t()                                = delete;
                        FFmc_c8_t(const FFmc_c8_t &)               = delete;
                        FFmc_c8_t & operator=(const FFmc_c8_t &)   = delete;
                        struct eikr_c8_t                           eikr;
                        struct Nmsph_c8_t            * __restrict  nmsp;
                        struct Nmv_c8_t              * __restrict  nmv;
                        double                       * __restrict  e_r;
                        double                       * __restrict  e_th;
                        double                       * __restrict  e_ph;
                        std::complex<double>         * __restrict  ffmc_x;
                        std::complex<double>         * __restrict  ffmc_y;
                        std::complex<double>         * __restrict  ffmc_z;
                        std::complex<double>                       ifac;
                        int32_t                                    npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } FFmc_c8_t;


           //  ! Formula (2-38)
           typedef struct __ATTR_ALIGN__(64) poynting_avg_r4_t {
                        poynting_avg_r4_t()                             = delete;
                        poynting_avg_r4_t(const poynting_avg_r4_t &)    = delete;
                        poynting_avg_r4_t & operator=(const poynting_avg_r4_t &) = delete;
                        struct Nmsph_r4_t             * __restrict nmsp;
                        struct Nesph_r4_t             * __restrict nesp;
                        float                         * __restrict R;
                        float                         * __restrict e_r;
                        float                         * __restrict S_re;
                        float                         * __restrict S_im;
                        float                                      k2;
                        float                                      eps_r;
                        float                                      eps_i;
                        float                                      mu_r;
                        float                                      mu_i;
                        int32_t                                    nth;
                        int32_t                                    nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,52)
#endif 
           } poynting_avg_r4_t;


             //  ! Formula (2-38)
           typedef struct __ATTR_ALIGN__(64) poynting_avg_r8_t {
                        poynting_avg_r8_t()                             = delete;
                        poynting_avg_r8_t(const poynting_avg_r8_t &)    = delete;
                        poynting_avg_r8_t & operator=(const poynting_avg_r8_t &) = delete;
                        struct Nmsph_r8_t             * __restrict nmsp;
                        struct Nesph_r8_t             * __restrict nesp;
                        double                        * __restrict R;
                        double                        * __restrict e_r;
                        double                        * __restrict S_re;
                        double                        * __restrict S_im;
                        double                                     k2;
                        double                                     eps_r;
                        double                                     eps_i;
                        double                                     mu_r;
                        double                                     mu_i;
                        int32_t                                    nth;
                        int32_t                                    nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,32)
#endif 
           } poynting_avg_r8_t;


           
           //  ! Formula (2-38)
           typedef struct __ATTR_ALIGN__(64) poynting_avg_c4_t {
                        poynting_avg_c4_t()                             = delete;
                        poynting_avg_c4_t(const poynting_avg_c4_t &)    = delete;
                        poynting_avg_c4_t & operator=(const poynting_avg_c4_t &) = delete;
                        struct Nmsph_c4_t             * __restrict nmsp;
                        struct Nesph_c4_t             * __restrict nesp;
                        float                         * __restrict R;
                        float                         * __restrict e_r;
                        std::complex<float>           * __restrict S;
                        std::complex<float>                        eps;
                        std::complex<float>                        mu;
                        float                                      k2;
                        int32_t                                    nth;
                        int32_t                                    nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,60)
#endif 
           } poynting_avg_c4_t;


           //  ! Formula (2-38)
           typedef struct __ATTR_ALIGN__(64) poynting_avg_c8_t {
                        poynting_avg_c8_t()                             = delete;
                        poynting_avg_c8_t(const poynting_avg_c8_t &)    = delete;
                        poynting_avg_c8_t & operator=(const poynting_avg_c8_t &) = delete;
                        struct Nmsph_c8_t             * __restrict nmsp;
                        struct Nesph_c8_t             * __restrict nesp;
                        double                        * __restrict R;
                        double                        * __restrict e_r;
                        std::complex<double>          * __restrict S;
                        std::complex<double>                       eps;
                        std::complex<double>                       mu;
                        double                                     k2;
                        int32_t                                    nth;
                        int32_t                                    nph;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,40)
#endif 
           } poynting_avg_c8_t;


           //  ! Formula (2-37)
           //! Intensity of radiation (steradian angle) in directions: phi,theta
           typedef struct __ATTR_ALIGN__(64) f237_r4_t {
                        f237_r4_t()                                = delete;
                        f237_r4_t(const f237_r4_t &)               = delete;
                        f237_r4_t & operator=(const f237_r4_t &)   = delete;
                        struct poynting_avg_r4_t                   pavg;
                        float                         * __restrict P_re;
                        float                         * __restrict P_im;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,54)
#endif 
           } f237_r4_t;


           //  ! Formula (2-37)
           //! Intensity of radiation (steradian angle) in directions: phi,theta
           typedef struct __ATTR_ALIGN__(64) f237_r8_t {
                        f237_r8_t()                                = delete;
                        f237_r8_t(const f237_r8_t &)               = delete;
                        f237_r8_t & operator=(const f237_r8_t &)   = delete;
                        struct poynting_avg_r8_t                   pavg;
                        double                         * __restrict P_re;
                        double                         * __restrict P_im;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,54)
#endif 
           } f237_r8_t;


           //  ! Formula (2-37)
           //! Intensity of radiation (steradian angle) in directions: phi,theta
           typedef struct __ATTR_ALIGN__(64) f237_c4_t {
                        f237_c4_t()                                = delete;
                        f237_c4_t(const f237_c4_t &)               = delete;
                        f237_c4_t & operator=(const f237_c4_t &)   = delete;
                        struct poynting_avg_c4_t                   pavg;
                        std::complex<float>           * __restrict P;
                        
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,62)
#endif 
           } f237_c4_t;


           //  ! Formula (2-37)
           //! Intensity of radiation (steradian angle) in directions: phi,theta
           typedef struct __ATTR_ALIGN__(64) f237_c8_t {
                        f237_c8_t()                                = delete;
                        f237_c8_t(const f237_c8_t &)               = delete;
                        f237_c8_t & operator=(const f237_c8_t &)   = delete; 
                        struct poynting_avg_c8_t                   pavg;
                        std::complex<double>          * __restrict P;
                        
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,62)
#endif 
           } f237_c8_t;


          /*
                 ! Formula (2-39)
                 ! Normalized radiation pattern expressed by the ratio
                 ! of poynting vectors.
           */
           typedef struct __ATTR_ALIGN__(64) f239_r4_t {
                         f239_r4_t()                                = delete;
                         f239_r4_t(const f239_r4_t &)               = delete;
                         f239_r4_t & operator=(const f239_r4_t &)   = delete;
                         struct poynting_avg_r4_t                  pv;
                         struct poynitng_avg_r4_t                  pvmax;
                         float                        * __restrict psi_re;
                         float                        * __restrict psi_im; //272 bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,48)
#endif 
           } f239_r4_t;


            /*
                 ! Formula (2-39)
                 ! Normalized radiation pattern expressed by the ratio
                 ! of poynting vectors.
           */
           typedef struct __ATTR_ALIGN__(64) f239_r8_t {
                         f239_r8_t()                                = delete;
                         f239_r8_t(const f239_r8_t &)               = delete;
                         f239_r8_t & operator=(const f239_r8_t &)   = delete;
                         struct poynting_avg_r8_t                   pv;
                         struct poynitng_avg_r8_t                   pvmax;
                         double                        * __restrict psi_re;
                         double                        * __restrict psi_im; //272 bytes
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,48)
#endif 
           } f239_r8_t;


              /*
                 ! Formula (2-39)
                 ! Normalized radiation pattern expressed by the ratio
                 ! of poynting vectors.
           */
           typedef struct __ATTR_ALIGN__(64) f239_c4_t {
                         f239_c4_t()                                = delete;
                         f239_c4_t(const f239_c4_t &)               = delete;
                         f239_c4_t & operator=(const f239_c4_t &)   = delete;
                         struct poynting_avg_c4_t                  pv;
                         struct poynitng_avg_c4_t                  pvmax;
                         std::complex<float>          * __restrict psi;
                         
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,56)
#endif 
           } f239_c4_t;


            /*
                 ! Formula (2-39)
                 ! Normalized radiation pattern expressed by the ratio
                 ! of poynting vectors.
           */
           typedef struct __ATTR_ALIGN__(64) f239_c8_t {
                         f239_c8_t()                                = delete;
                         f239_c8_t(const f239_r4_t &)               = delete;
                         f239_c8_t & operator=(const f239_c8_t &)   = delete; 
                         struct poynting_avg_c8_t                   pv;
                         struct poynitng_avg_c8_t                   pvmax;
                         std::complex<double>          * __restrict psi;
                         
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,56)
#endif 
           } f239_c8_t;

           /*
                ! Formula (2-43)
                ! The current distribution along the 'z' coordiante
            */
            typedef struct __ATTR_ALIGN__(32) f243_r4_t {
                          f243_r4_t()                                = delete;
                          f243_r4_t(const f243_r4_t &)               = delete;
                          f243_r4_t & operator=(const f243_r4_t &)   = delete;
                          float                        * __restrict Iz;
                          float                                     I0;
                          float                                     k;
                          int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } f243_r4_t;


             /*
                ! Formula (2-43)
                ! The current distribution along the 'z' coordiante
            */
            typedef struct __ATTR_ALIGN__(32) f243_r8_t {
                          f243_r8_t()                                = delete;
                          f243_r8_t(const f243_r8_t &)               = delete;
                          f243_r8_t & operator=(const f243_r8_t &)   = delete;
                          double                        * __restrict Iz;
                          double                                     I0;
                          double                                     k;
                          int32_t                                    npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif 
            } f243_r8_t;


            /*
                ! Formula (2-44)
                ! Radiation pattern of thin wire (odd m)
             */
            typedef struct __ATTR_ALIGN__(32) f244_r4_t {
                           f244_r4_t()                                = delete;
                           f244_r4_t(const f244_r4_t &)               = delete;
                           f244_r4_t & operator=(const f244_r4_t &)   = delete;
                           float                        * __restrict Fth;
                           float                                     A;
                           int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,16)
#endif 
            } f244_r4_t;


               /*
                ! Formula (2-44)
                ! Radiation pattern of thin wire (odd m)
             */
            typedef struct __ATTR_ALIGN__(32) f244_r8_t {
                           f244_r8_t()                                = delete;
                           f244_r8_t(const f244_r8_t &)               = delete;
                           f244_r8_t & operator=(const f244_r8_t &)   = delete;
                           double                        * __restrict Fth;
                           double                                     A;
                           int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif 
            } f244_r8_t;


           /*
               ! Formula (2-45)
               ! Radiation pattern of thin wire (even m)
            */ 
            typedef struct __ATTR_ALIGN__(32) f245_r4_t {
                           f245_r4_t()                                = delete;
                           f245_r4_t(const f245_r4_t &)               = delete;
                           f245_r4_t & operator=(const f245_r4_t &)   = delete;
                           float                        * __restrict Fth;
                           float                                     A;
                           int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,16)
#endif  
           } f245_r4_t;


           /*
               ! Formula (2-45)
               ! Radiation pattern of thin wire (even m)
            */ 
            typedef struct __ATTR_ALIGN__(32) f245_r8_t {
                           f245_r8_t()                                = delete;
                           f245_r8_t(const f245_r8_t &)               = delete;
                           f245_r8_t & operator=(const f245_r8_t &)   = delete; 
                           double                        * __restrict Fth;
                           double                                     A;
                           int32_t                                   npts;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif  
           } f245_r8_t;


           /*
               ! Formula (2-48)
               ! For a wire with current modulated by phase speed
           */
           typedef struct __ATTR_ALIGN__(32) f248_r4_t {
                            f248_r4_t()                                = delete;
                            f248_r4_t(const f248_r4_t &)               = delete;
                            f248_r4_t & operator=(const f248_r4_t &)   = delete; 
                            float                       * __restrict Iz_re;
                            float                       * __restrict Iz_im;
                            float                                    I0_re;
                            float                                    I0_im;
                            float                                    beta
                            int32_t                                  nz;
 
           } f248_r4_t;


          /*
               ! Formula (2-48)
               ! For a wire with current modulated by phase speed
           */
           typedef struct __ATTR_ALIGN__(64) f248_r8_t {
                            f248_r8_t()                                = delete;
                            f248_r8_t(const f248_r8_t &)               = delete;
                            f248_r8_t & operator=(const f248_r8_t &)   = delete; 
                            double                       * __restrict Iz_re;
                            double                       * __restrict Iz_im;
                            double                                    I0_re;
                            double                                    I0_im;
                            double                                    beta
                            int32_t                                   nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif  
           } f248_r8_t;


           /*
               ! Formula (2-48)
               ! For a wire with current modulated by phase speed
           */
           typedef struct __ATTR_ALIGN__(32) f248_c4_t {
                            f248_c4_t()                                = delete;
                            f248_c4_t(const f248_c4_t &)               = delete;
                            f248_c4_t & operator=(const f248_c4_t &)   = delete; 
                            std::complex<float>          * __restrict Iz;
                            std::complex<float>                       I0;
                            float                                     beta
                            int32_t                                   nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif   
           } f248_c4_t;


          /*
               ! Formula (2-48)
               ! For a wire with current modulated by phase speed
           */
           typedef struct __ATTR_ALIGN__(64) f248_c8_t {
                            f248_c8_t()                                = delete;
                            f248_c8_t(const f248_c8_t &)               = delete;
                            f248_c8_t & operator=(const f248_c8_t &)   = delete; 
                            std::complex<double>          * __restrict Iz;
                            std::complex<double>                       I0;
                            double                                     beta
                            int32_t                                    nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,28)
#endif  
           } f248_c8_t;


           /*
               ! Formula (2-46)
               ! Symmetric sinusoidal current distribution.
            */
           typedef struct __ATTR_ALIGN__(64) f246_r4_t {
                            f246_r4_t()                                = delete;
                            f246_r4_t(const f246_r4_t &)               = delete;
                            f246_r4_t & operator=(const f246_r4_t &)   = delete; 
                            float                         * __restrict Iz_re;
                            float                         * __restrict Iz_im;
                            float                                      I0_re;
                            float                                      I0_im;
                            float                                      L;
                            float                                      k;
                            int32_t                                    nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,28)
#endif  
          } f246_r4_t;


           /*
               ! Formula (2-46)
               ! Symmetric sinusoidal current distribution.
            */
           typedef struct __ATTR_ALIGN__(64) f246_r8_t {
                            f246_r8_t()                                = delete;
                            f246_r8_t(const f246_r8_t &)               = delete;
                            f246_r8_t & operator=(const f246_r8_t &)   = delete; 
                            double                         * __restrict Iz_re;
                            double                         * __restrict Iz_im;
                            double                                      I0_re;
                            double                                      I0_im;
                            double                                      L;
                            double                                      k;
                            int32_t                                    nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,12)
#endif  
          } f246_r8_t;


             /*
               ! Formula (2-46)
               ! Symmetric sinusoidal current distribution.
            */
           typedef struct __ATTR_ALIGN__(32) f246_c4_t {
                            f246_c4_t()                                = delete;
                            f246_c4_t(const f246_c4_t &)               = delete;
                            f246_c4_t & operator=(const f246_c4_t &)   = delete; 
                            std::complex<float>            * __restrict Iz;
                            std::complex<float>                         I0;
                            float                                       L;
                            float                                       k;
                            int32_t                                     nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif  
          } f246_c4_t;


           /*
               ! Formula (2-46)
               ! Symmetric sinusoidal current distribution.
            */
           typedef struct __ATTR_ALIGN__(64) f246_c8_t {
                            f246_c8_t()                                = delete;
                            f246_c8_t(const f246_c8_t &)               = delete;
                            f246_c8_t & operator=(const f246_c8_t &)   = delete; 
                            std::complex<double>           * __restrict Iz;
                            std::complex<double>                        I0;
                            double                                      L;
                            double                                      k;
                            int32_t                                    nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif  
          } f246_c8_t;


          /*
              ! Formula (2-47)
              ! Wire symmetric current radiation pattern
          */
           typedef struct __ATTR_ALIGN__(32) f247_r4_t {
                            f247_r4_t()                                = delete;
                            f247_r4_t(const f247_r4_t &)               = delete;
                            f247_r4_t & operator=(const f247_r4_t &)   = delete; 
                            float                          * __restrict Fth;
                            float                                       k;
                            float                                       L;
                            float                                       A;
                            int32_t                                     nth;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,8)
#endif  
           } f247_f4_t;


           /*
              ! Formula (2-47)
              ! Wire symmetric current radiation pattern
          */
           typedef struct __ATTR_ALIGN__(64) f247_r8_t {
                            f247_r8_t()                                = delete;
                            f247_r8_t(const f247_r8_t &)               = delete;
                            f247_r8_t & operator=(const f247_r8_t &)   = delete; 
                            double                          * __restrict Fth;
                            double                                       k;
                            double                                       L;
                            double                                       A;
                            int32_t                                     nth;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,28)
#endif  
           } f247_f8_t;


           /*
               ! Formula (2-49)
               ! Wire (running-wave) radiation pattern
           */
          typedef struct __ATTR_ALIGN__(32) f249_r4_t {
                            f249_r4_t()                                = delete;
                            f249_r4_t(const f249_r4_t &)               = delete;
                            f249_r4_t & operator=(const f249_r4_t &)   = delete; 
                            float                          * __restrict Fth;
                            float                                       k;
                            float                                       L;
                            float                                       A;
                            float                                       beta;
                            int32_t                                     nth;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,4)
#endif   
          } f249_r4_t;


           /*
               ! Formula (2-49)
               ! Wire (running-wave) radiation pattern
           */
          typedef struct __ATTR_ALIGN__(64) f249_r8_t {
                            f249_r8_t()                                = delete;
                            f249_r8_t(const f249_r8_t &)               = delete;
                            f249_r8_t & operator=(const f249_r8_t &)   = delete;  
                            double                          * __restrict Fth;
                            double                                       k;
                            double                                       L;
                            double                                       A;
                            double                                       beta;
                            int32_t                                      nth;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,20)
#endif  
          } f249_r8_t;


          /*
              
          */


       } // radiolocation

} // gms














#endif /*__GMS_ANTENNA_TYPES_V2_H__*/
