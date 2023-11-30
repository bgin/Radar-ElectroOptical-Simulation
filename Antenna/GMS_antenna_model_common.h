

#ifndef __GMS_ANTENNA_MODEL_COMMON_H__
#define __GMS_ANTENNA_MODEL_COMMON_H__ 201120221117


namespace file_info {

     const unsigned int GMS_ANTENNA_MODEL_COMMON_MAJOR = 1;
     const unsigned int GMS_ANTENNA_MODEL_COMMON_MINOR = 1;
     const unsigned int GMS_ANTENNA_MODEL_COMMON_MICRO = 0;
     const unsigned int GMS_ANTENNA_MODEL_COMMON_FULLVER =
       1000U*GMS_ANTENNA_MODEL_COMMON_MAJOR+100U*GMS_ANTENNA_MODEL_COMMON_MINOR+
       10U*GMS_ANTENNA_MODEL_COMMON_MICRO;
     const char * const GMS_ANTENNA_MODEL_COMMON_CREATION_DATE = "20-11-2022 11:17 +00200 (SUN 20 NOV 2022 GMT+2)";
     const char * const GMS_ANTENNA_MODEL_COMMON_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_ANTENNA_MODEL_COMMON_SYNOPSIS      = "Antenna model dynamically allocated abstract data types."

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
#ifdef __INTEL_COMPILER
#include <.../.../intel/oneapi/compiler/2023.1.0/linux/compiler/perf_headers/c++/valarray.h>
#else
#include <valarray>
#endif
#include <cstring> // std::memcpy
#include "GMS_config"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)
#define USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES 0
#endif

namespace gms {


          namespace  radiolocation {
          

              struct __ATTR_ALIGN__(64) E_c4_t {
                      // Complex electric field.
                      std::complex<float> * __restrict ex
                      std::complex<float> * __restrict ey;
                      std::complex<float> * __restrict ez;
                      std::size_t                      nx;
                      std::size_t                      ny;
                      std::size_t                      nz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline E_c4_t() noexcept(true) {
                          
                          this->nx  = 0ULL;
                          this->ny  = 0ULL;
                          this->nz  = 0ULL;
                          this->ex  = NULL;
                          this->ey  = NULL;
                          this->ez  = NULL;
                      } 
                          
                     inline E_c4_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz) {
                          using namespace gms::common;
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline E_c4_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->ex = (std::complex<float>*)
                                                 gms_mmap_4KiB( sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->ey = (std::complex<float>*)
                                                 gms_mmap_4KiB( sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->ez = (std::complex<float>*)
                                                 gms_mmap_4KiB( sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->ex = (std::complex<float>*)
                                                 gms_mmap_2MiB( sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->ey = (std::complex<float>*)
                                                 gms_mmap_2MiB( sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->ez = (std::complex<float>*)
                                                 gms_mmap_2MiB( sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->ex = (std::complex<float>*)
                                                 gms_mmap_1GiB( sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->ey = (std::complex<float>*)
                                                 gms_mmap_1GiB( sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->ez = (std::complex<float>*)
                                                 gms_mmap_1GiB( sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset); 
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
                          
                          this->nx = e_x.size();
                          this->ny = e_y.size();
                          this->nz = e_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->ex,&e_x[0],lenx);
                          std::memcpy(this->ey,&e_y[0],leny);
                          std::memcpy(this->ez,&e_z[0],lenz);       
                     }
                     
                     inline E_c4_t(const std::valarray<std::complex<float>> &e_x,
                                   const std::valarray<std::complex<float>> &e_y,
                                   const std::valarray<std::complex<float>> &e_z) {
                          using namespace gms::common;
                          
                          this->nx = e_x.size();
                          this->ny = e_y.size();
                          this->nz = e_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->ex,&e_x[0],lenx);
                          std::memcpy(this->ey,&e_y[0],leny);
                          std::memcpy(this->ez,&e_z[0],lenz);           
                     }
                             
                             
                      
                     inline E_c4_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz,
                                   const std::complex<float> * __restrict e_x,
                                   const std::complex<float> * __restrict e_y,
                                   const std::complex<float> * __restrict e_z) {
                          using namespace gms::common;
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->ex[0],&e_x[0],this->nx);
	                  avx512_uncached_memmove(&this->ey[0],&e_y[0],this->ny);
	                  avx512_uncached_memmove(&this->ez[0],&e_z[0],this->nz);
#else
	                  avx512_cached_memmove(&this->ex[0],&e_x[0],this->nx);
	                  avx512_cached_memmove(&this->ey[0],&e_y[0],this->ny);
	                  avx512_cached_memmove(&this->ez[0],&e_z[0],this->nz);
#endif
                   }  
                   
                    inline  E_c4_t(E_c4_t && rhs) {
                          
                          this->nx = rhs.nx;
                          this->ny = rhs.ny;
                          this->nz = rhs.nz;
                          this->ex   = &rhs.ex[0];
                          this->ey   = &rhs.ey[0];
                          this->ez   = &rhs.ez[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
                          rhs.ex     = NULL;
                          rhs.ey     = NULL;
                          rhs.ez     = NULL;
                      }
                                 
                     E_c4_t(const E_c4_t &)             = delete;
                      
                     inline ~E_c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap(this->ex,this->nx);
                             gms_unmap(this->ey,this->ny);
                             gms_unmap(this->ez,this->nz);
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
                           this->nx   = rhs.nx;
                           this->ny   = rhs.ny;
                           this->nz   = rhs.nz;
                           this->ex   = &rhs.ex[0];
                           this->ey   = &rhs.ey[0];
                           this->ez   = &rhs.ez[0];
                           rhs.nx     = 0ULL;
                           rhs.ny     = 0ULL;
                           rhs.nz     = 0ULL;
                           rhs.ex     = NULL;
                           rhs.ey     = NULL;
                           rhs.ez     = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->ex  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->nx,64ULL);
                          this->ey  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->ny,64ULL);
                          this->ez  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->ny,64ULL);
                      }
                    
              };

              
                struct __ATTR_ALIGN__(32) H_c4_t {
                     
                      std::complex<float> * __restrict mx
                      std::complex<float> * __restrict my;
                      std::complex<float> * __restrict mz;
                      std::size_t                      nx;
                      std::size_t                      ny;
                      std::size_t                      nz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                      inline H_c4_t() noexcept(true) {
                          
                          this->nx  = 0ULL;
                          this->ny  = 0ULL;
                          this->nz  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                          this->mz  = NULL;
                      } 
                          
                      inline H_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                      inline H_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_4KiB( sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_2MiB( sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_2MiB( sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_2MiB( sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_1GiB( sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_1GiB( sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_1GiB( sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset); 
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
                                                   
                          this->nx = m_x.size();
                          this->ny = m_y.size();
                          this->nz = m_z.size();
                          allocate()
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->mx,&m_x[0],lenx);
                          std::memcpy(this->my,&m_y[0],leny);
                          std::memcpy(this->mz,&m_z[0],lenz);       
                     }
                     
                     inline H_c4_t(const std::valarray<std::complex<float>> &m_x,
                                   const std::valarray<std::complex<float>> &m_y,
                                   const std::valarray<std::complex<float>> &m_z) {
                          using namespace gms::common;
                          
                          this->nx = m_x.size();
                          this->ny = m_y.size();
                          this->nz = m_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->mx,&m_x[0],lenx);
                          std::memcpy(this->my,&m_y[0],leny);
                          std::memcpy(this->mz,&m_z[0],lenz);           
                     }
                      
                     inline H_c4_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz,
                                   const std::complex<float> * __restrict m_x,
                                   const std::complex<float> * __restrict m_y,
                                   const std::complex<float> * __restrict m_z) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->nx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->ny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->nz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->nx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->ny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->nz);
#endif
                   }  
                   
                    inline  H_c4_t(H_c4_t && rhs) {
                          
                          this->nx   = rhs.nx;
                          this->ny   = rhs.ny;
                          this->nz   = rhs.nz;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          this->mz   = &rhs.mz[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                          rhs.mz     = NULL;
                      }
                                 
                    inline  H_c4_t(const H_c4_t &)             = delete;
                      
                    inline  ~H_c4_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->mx,this->nx);
                             gms_unmap(this->my,this->nx);
                             gms_unmap(this->mz,this->nx);
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
                           this->nx   = rhs.nx;
                           this->ny   = rhs.ny;
                           this->nz   = rhs.nz;
                           this->mx   = &rhs.mx[0];
                           this->my   = &rhs.my[0];
                           this->mz   = &rhs.mz[0];
                           rhs.nx     = 0ULL;
                           rhs.ny     = 0ULL;
                           rhs.nz     = 0ULL;
                           rhs.mx     = NULL;
                           rhs.my     = NULL;
                           rhs.mz     = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->nx,64ULL);
                          this->my  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->ny,64ULL);
                          this->mz  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->nz,64ULL); 
                      }
                    
              };


                struct __ATTR_ALIGN__(32) E_c8_t {
                      // Complex electric field.
                      std::complex<double> * __restrict ex
                      std::complex<double> * __restrict ey;
                      std::complex<double> * __restrict ez;
                      std::size_t                       nx;
                      std::size_t                       ny;
                      std::size_t                       nz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline E_c8_t() noexcept(true) {
                          
                          this->nx  = 0ULL;
                          this->ny  = 0ULL;
                          this->nz  = 0ULL;
                          this->ex  = NULL;
                          this->ey  = NULL;
                          this->ez  = NULL;
                      } 
                          
                     inline E_c8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz) {
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline E_c8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->ex = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(td::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->ey = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->ez = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->ex = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->ey = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->ez = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->ex = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->ey = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->ez = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset); 
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
                          
                          this->nx = e_x.size();
                          this->ny = e_y.size();
                          this->nz = e_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->ex,&e_x[0],lenx);
                          std::memcpy(this->ey,&e_y[0],leny);
                          std::memcpy(this->ez,&e_z[0],lenz);       
                     }
                     
                     
                     inline E_c8_t(const std::valarray<std::complex<double>> &e_x,    //shall be of the same size (no error checking implemented)
                                   const std::valarray<std::complex<double>> &e_y,    //shall be of the same size (no error checking implemented)
                                   const std::valarray<std::complex<double>> &e_z) {  //shall be of the same size (no error checking implemented)
                          
                          this->nx = e_x.size();
                          this->ny = e_y.size();
                          this->nz = e_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->ex,&e_x[0],lenx);
                          std::memcpy(this->ey,&e_y[0],leny);
                          std::memcpy(this->ez,&e_z[0],lenz);       
                     }      
                      
                     inline E_c8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz,
                                   const std::complex<double> * __restrict e_x,
                                   const std::complex<double> * __restrict e_y,
                                   const std::complex<double> * __restrict e_z) {
                                   
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->ex[0],&e_x[0],this->nx);
	                  avx512_uncached_memmove(&this->ey[0],&e_y[0],this->ny);
	                  avx512_uncached_memmove(&this->ez[0],&e_z[0],this->nz);
#else
	                  avx512_cached_memmove(&this->ex[0],&e_x[0],this->nx);
	                  avx512_cached_memmove(&this->ey[0],&e_y[0],this->ny);
	                  avx512_cached_memmove(&this->ez[0],&e_z[0],this->nz);
#endif
                   }  
                   
                   inline   E_c8_t(E_c8_t && rhs) {
                          
                          this->nx   = rhs.nx;
                          this->ny   = rhs.ny;
                          this->nz   = rhs.nz;
                          this->ex   = &rhs.ex[0];
                          this->ey   = &rhs.ey[0];
                          this->ez   = &rhs.ez[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
                          rhs.ex     = NULL;
                          rhs.ey     = NULL;
                          rhs.ez     = NULL;
                      }
                                 
                      E_c8_t(const E_c48t &)             = delete;
                      
                    inline ~E_c8_t() {
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->ex,this->nx);
                             gms_unmap(this->ey,this->ny);
                             gms_unmap(this->ez,this->nz);
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
                           this->nx   = rhs.nx;
                           this->ny   = rhs.ny;
                           this->nz   = rhs.nz;
                           this->ex   = &rhs.ex[0];
                           this->ey   = &rhs.ey[0];
                           this->ez   = &rhs.ez[0];
                           rhs.nx     = 0ULL;
                           rhs.ny     = 0ULL;
                           rhs.nz     = 0ULL;
                           rhs.ex     = NULL;
                           rhs.ey     = NULL;
                           rhs.ez     = NULL;
                           return (*this);
                      }
                      
                      void allocate() {
                            using namespace gms::common;
                            this->ex  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nx,64ULL);
                            this->ey  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->ny,64ULL);
                            this->ez  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nz,64ULL);
                      }
              };

              
               struct __ATTR_ALIGN__(32) H_c8_t {
                      // Complex magnetic field.
                      std::complex<double> * __restrict mx
                      std::complex<double> * __restrict my;
                      std::complex<double> * __restrict mz;
                      std::size_t                       nx;
                      std::size_t                       ny;
                      std::size_t                       nz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline H_c8_t() noexcept(true) {
                          
                          this->nx  = 0ULL;
                          this->ny  = 0ULL;
                          this->nz  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                          this->mz  = NULL;
                      } 
                          
                     inline H_c8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz) {
                     
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline H_c8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->nx  = _nx;
                             this->ny  = _ny;
                             this->nz  = _nz;
                             switch (fsize) {
                                 case:0
                                      this->mx = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->my = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->mz = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->my = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->mz = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->my = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->mz = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset); 
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
                         
                          this->nx = m_x.size();
                          this->ny = m_y.size();
                          this->nz = m_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->mx,&m_x[0],lenx);
                          std::memcpy(this->my,&m_y[0],leny);
                          std::memcpy(this->mz,&m_z[0],lenz);       
                     }
                     
                     
                     inline   H_c8_t(const std::valarray<std::complex<double>> &m_x,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<double>> &m_y,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<double>> &m_z) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = m_x.size();
                          this->ny = m_y.size();
                          this->nz = m_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->mx,&m_x[0],lenx);
                          std::memcpy(this->my,&m_y[0],leny);
                          std::memcpy(this->mz,&m_z[0],lenz);       
                     } 
                     
                      
                   inline   H_c8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz
                                   const std::complex<double> * __restrict m_x,
                                   const std::complex<double> * __restrict m_y,
                                   const std::complex<double> * __restrict m_z) {
                          
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->nx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->ny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->nz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->nx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->ny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->nz);
#endif
                   }  
                   
                  inline  H_c8_t(H_c8_t && rhs) {
                          
                          this->nx   = rhs.nx;
                          this->ny   = rhs.ny;
                          this->nz   = rhs.nz;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          this->mz   = &rhs.mz[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                          rhs.mz     = NULL;
                      }
                                 
                      H_c8_t(const H_c8_t &)             = delete;
                      
                   inline   ~H_c8_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->mx,this->nx);
                             gms_unmap(this->my,this->ny);
                             gms_unmap(this->mz,this->nz);
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
                           this->nx   = rhs.nx;
                           this->ny   = rhs.ny;
                           this->nz   = rhs.nz;
                           this->mx   = &rhs.mx[0];
                           this->my   = &rhs.my[0];
                           this->mz   = &rhs.mz[0];
                           rhs.nx     = 0ULL;
                           rhs.ny     = 0ULL;
                           rhs.nz     = 0ULL;
                           rhs.mx     = NULL;
                           rhs.my     = NULL;
                           rhs.mz     = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() {
                        using namespace gms::common;
                        this->mx  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nx,64ULL);
                        this->my  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->ny,64ULL);
                        this->mz  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nz,64ULL);
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
                      std::size_t        nx;
                      std::size_t        ny;
                      std::size_t        nz;
                      bool              ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif
                      inline E_r4_t() {
                         this->nx   = 0ULL;
                         this->ny   = 0ULL;
                         this->nz   = 0ULL;
                         this->exr  = NULL;
                         this->exi  = NULL;
                         this->eyr  = NULL;
                         this->eyi  = NULL;
                         this->ezr  = NULL;
                         this->ezi  = NULL;
                      }                    
                      
                      inline E_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline E_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->exr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->exi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->eyr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->eyi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->ezr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ezi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->exr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->exi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->eyr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->eyi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->ezr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ezi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->exr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->exi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->eyr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->eyi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->ezr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ezi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
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
                               
                               this->nx = e_xr.size(); 
                               this->ny = e_yr.size();
                               this->nz = e_zr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenx = sizeof(float)*this->nz;
                               std::memcpy(this->exr,&e_xr[0],lenx);
                               std::memcpy(this->exi,&e_xi[0],lenx);
                               std::memcpy(this->eyr,&e_yr[0],leny);
                               std::memcpy(this->eyi,&e_yi[0],leny);
                               std::memcpy(this->ezr,&e_zr[0],lenz);
                               std::memcpy(this->ezi,&e_zi[0],lenz);     
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline E_r4_t(const std::valarray<float> &e_xr,
                                    const std::valarray<float> &e_xi,
                                    const std::valarray<float> &e_yr,
                                    const std::valarray<float> &e_yi,
                                    const std::valarray<float> &e_zr,
                                    const std::valarray<float> &e_zi) {
                               
                               this->nx = e_xr.size(); 
                               this->ny = e_yr.size();
                               this->nz = e_zr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenx = sizeof(float)*this->nz;
                               std::memcpy(this->exr,&e_xr[0],lenx);
                               std::memcpy(this->exi,&e_xi[0],lenx);
                               std::memcpy(this->eyr,&e_yr[0],leny);
                               std::memcpy(this->eyi,&e_yi[0],leny);
                               std::memcpy(this->ezr,&e_zr[0],lenz);
                               std::memcpy(this->ezi,&e_zi[0],lenz);     
                      }
                      
                      inline E_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const float * __restrict e_xr,   
                                    const float * __restrict e_xi,
                                    const float * __restrict e_yr,
                                    const float * __restrict e_yi,
                                    const float * __restrict e_zr,
                                    const float * __restrict e_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->exr[0],&e_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->exi[0],&e_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->eyr[0],&e_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->eyi[0],&e_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->ezr[0],&e_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->ezi[0],&e_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->exr[0],&e_xr[0],this->nx);
	                  avx512_cached_memmove(&this->exi[0],&e_xi[0],this->nx);
	                  avx512_cached_memmove(&this->eyr[0],&e_yr[0],this->ny);
	                  avx512_cached_memmove(&this->eyi[0],&e_yi[0],this->ny);
	                  avx512_cached_memmove(&this->ezr[0],&e_zr[0],this->nz);
	                  avx512_cached_memmove(&this->ezi[0],&e_zi[0],this->nz);
#endif          
                      }
                      
                      inline E_r4_t(E_r4_t &&rhs) {
                           
                          this->nx   = rhs.nx;
                          this->ny   = rhs.ny;
                          this->nz   = rhs.nz;
                          this->exr  = &rhs.exr[0];
                          this->exi  = &rhs.exi[0];
                          this->eyr  = &rhs.eyr[0];
                          this->eyi  = &rhs.eyi[0];
                          this->ezr  = &rhs.ezr[0];
                          this->ezi  = &rhs.ezi[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
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
                              gms_unmap(this->exr,this->nx);
                              gms_unmap(this->exi,this->nx);
                              gms_unmap(this->eyr,this->ny); 
                              gms_unmap(this->eyi,this->ny);
                              gms_unmap(this->ezr,this->nz);
                              gms_unmap(this->ezi,this->nz);
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
                            this->nx   = rhs.nx;
                            this->ny   = rhs.ny;
                            this->nz   = rhs.nz;
                            this->exr  = &rhs.exr[0];
                            this->exi  = &rhs.exi[0];
                            this->eyr  = &rhs.eyr[0];
                            this->eyi  = &rhs.eyi[0];
                            this->ezr  = &rhs.ezr[0];
                            this->ezi  = &rhs.ezi[0];
                            rhs.nx     = 0ULL;
                            rhs.ny     = 0ULL;
                            rhs.nz     = 0ULL;
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
                        this->exr  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->exi  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->eyr  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->eyi  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->ezr  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
                        this->ezi  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
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
                      std::size_t            nx;
                      std::size_t            ny;
                      std::size_t            nz;
                      bool                   ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif
                     
                     inline H_r4_t() {
                         this->nx   = 0ULL;
                         this->ny   = 0ULL;
                         this->nz   = 0ULL;
                         this->mxr  = NULL;
                         this->mxi  = NULL;
                         this->myr  = NULL;
                         this->myi  = NULL;
                         this->mzr  = NULL;
                         this->mzi  = NULL;
                      }                    
                      
                      inline H_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline H_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->mxr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mxr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mxr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline H_r4_t(const std::vector<float> &m_xr,
                                    const std::vector<float> &m_xi,
                                    const std::vector<float> &m_yr,
                                    const std::vector<float> &m_yi,
                                    const std::vector<float> &m_zr,
                                    const std::vector<float> &m_zi) {
                               
                               this->nx = m_xr.size();
                               this->ny = m_yr.size();
                               this->nz = m_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenz = sizeof(float)*this->nz;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                               std::memcpy(this->mzr,&m_zr[0],lenz);
                               std::memcpy(this->mzi,&m_zi[0],lenz);     
                      }
                      
                     inline H_r4_t( const std::valarray<float> &m_xr,
                                    const std::valarray<float> &m_xi,
                                    const std::valarray<float> &m_yr,
                                    const std::valarray<float> &m_yi,
                                    const std::valarray<float> &m_zr,
                                    const std::valarray<float> &m_zi) {
                               
                               this->nx = m_xr.size();
                               this->ny = m_yr.size();
                               this->nz = m_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenz = sizeof(float)*this->nz;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                               std::memcpy(this->mzr,&m_zr[0],lenz);
                               std::memcpy(this->mzi,&m_zi[0],lenz);     
                      }
                      
                      inline H_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const float * __restrict m_xr,   
                                    const float * __restrict m_xi,
                                    const float * __restrict m_yr,
                                    const float * __restrict m_yi,
                                    const float * __restrict m_zr,
                                    const float * __restrict m_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&m_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->mxi[0],&m_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->myr[0],&m_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->myi[0],&m_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->mzr[0],&m_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->mzi[0],&m_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->mxr[0],&m_xr[0],this->nx);
	                  avx512_cached_memmove(&this->mxi[0],&m_xi[0],this->nx);
	                  avx512_cached_memmove(&this->myr[0],&m_yr[0],this->ny);
	                  avx512_cached_memmove(&this->myi[0],&m_yi[0],this->ny);
	                  avx512_cached_memmove(&this->mzr[0],&m_zr[0],this->nz);
	                  avx512_cached_memmove(&this->mzi[0],&m_zi[0],this->nz);
#endif          
                      }
                      
                      inline H_r4_t(H_r4_t &&rhs) {
                           
                          this->nx   = rhs.nx;
                          this->ny   = rhs.ny;
                          this->nz   = rhs.nz;
                          this->mxr  = &rhs.mxr[0];
                          this->mxi  = &rhs.mxi[0];
                          this->myr  = &rhs.myr[0];
                          this->myi  = &rhs.myi[0];
                          this->mzr  = &rhs.mzr[0];
                          this->mzi  = &rhs.mzi[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
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
                              gms_unmap(this->mxr,this->nx);
                              gms_unmap(this->mxi,this->nx);
                              gms_unmap(this->myr,this->ny); 
                              gms_unmap(this->myi,this->ny);
                              gms_unmap(this->mzr,this->nz);
                              gms_unmap(this->mzi,this->nz);
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
                            this->nx   = rhs.nx;
                            this->ny   = rhs.ny;
                            this->nz   = rhs.nz;
                            this->mxr  = &rhs.mxr[0];
                            this->mxi  = &rhs.mxi[0];
                            this->myr  = &rhs.myr[0];
                            this->myi  = &rhs.myi[0];
                            this->mzr  = &rhs.mzr[0];
                            this->mzi  = &rhs.mzi[0];
                            rhs.nx     = 0ULL;
                            rhs.ny     = 0ULL;
                            rhs.nz     = 0ULL;
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
                        this->mxr  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->mxi  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->myr  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->myi  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->mzr  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
                        this->mzi  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
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
                      std::size_t              nx;
                      std::size_t              ny;
                      std::size_t              nz;
                      bool                     ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif

                     inline E_r8_t() {
                         this->nx   = 0ULL;
                         this->ny   = 0ULL;
                         this->nz   = 0ULL;
                         this->exr  = NULL;
                         this->exi  = NULL;
                         this->eyr  = NULL;
                         this->eyi  = NULL;
                         this->ezr  = NULL;
                         this->ezi  = NULL;
                      }                    
                      
                     inline E_r8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                     inline E_r8_t(const std::size_t _nx,
                                   const std::size_t _ny,
                                   const std::size_t _nz,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->exr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->exi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->eyr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->eyi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->ezr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ezi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->exr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->exi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->eyr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->eyi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->ezr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ezi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->exr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->exi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->eyr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this-ny,prot,flags,fd,offset);
                                      this->eyi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->ezr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ezi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
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
                               
                               this->nx = e_xr.size();
                               this->ny = e_yr.size();
                               this->nz = e_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->exr,&e_xr[0],lenx);
                               std::memcpy(this->exi,&e_xi[0],lenx);
                               std::memcpy(this->eyr,&e_yr[0],leny);
                               std::memcpy(this->eyi,&e_yi[0],leny);
                               std::memcpy(this->ezr,&e_zr[0],lenz);
                               std::memcpy(this->ezi,&e_zi[0],lenz);     
                      }
                      
                      inline E_r8_t(const std::valarray<double> &e_xr,
                                    const std::valarray<double> &e_xi,
                                    const std::valarray<double> &e_yr,
                                    const std::valarray<double> &e_yi,
                                    const std::valarray<double> &e_zr,
                                    const std::valarray<double> &e_zi) {
                               
                               this->nx = e_xr.size();
                               this->ny = e_yr.size();
                               this->nz = e_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->exr,&e_xr[0],lenx);
                               std::memcpy(this->exi,&e_xi[0],lenx);
                               std::memcpy(this->eyr,&e_yr[0],leny);
                               std::memcpy(this->eyi,&e_yi[0],leny);
                               std::memcpy(this->ezr,&e_zr[0],lenz);
                               std::memcpy(this->ezi,&e_zi[0],lenz);     
                      }
                      
                      inline E_r8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const double * __restrict e_xr,   
                                    const double * __restrict e_xi,
                                    const double * __restrict e_yr,
                                    const double * __restrict e_yi,
                                    const double * __restrict e_zr,
                                    const double * __restrict e_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->exr[0],&e_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->exi[0],&e_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->eyr[0],&e_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->eyi[0],&e_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->ezr[0],&e_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->ezi[0],&e_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->exr[0],&e_xr[0],this->nx);
	                  avx512_cached_memmove(&this->exi[0],&e_xi[0],this->nx);
	                  avx512_cached_memmove(&this->eyr[0],&e_yr[0],this->ny);
	                  avx512_cached_memmove(&this->eyi[0],&e_yi[0],this->ny);
	                  avx512_cached_memmove(&this->ezr[0],&e_zr[0],this->nz);
	                  avx512_cached_memmove(&this->ezi[0],&e_zi[0],this->nz);
#endif          
                      }
                      
                      inline E_r8_t(E_r8_t &&rhs) {
                           
                          this->nx   = rhs.nx;
                          this->ny   = rhs.ny;
                          this->nz   = rhs.nz;
                          this->exr  = &rhs.exr[0];
                          this->exi  = &rhs.exi[0];
                          this->eyr  = &rhs.eyr[0];
                          this->eyi  = &rhs.eyi[0];
                          this->ezr  = &rhs.ezr[0];
                          this->ezi  = &rhs.ezi[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
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
                              gms_unmap(this->exr,this->nx);
                              gms_unmap(this->exi,this->nx);
                              gms_unmap(this->eyr,this->ny); 
                              gms_unmap(this->eyi,this->ny);
                              gms_unmap(this->ezr,this->nz);
                              gms_unmap(this->ezi,this->nz);
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
                            this->nx   = rhs.nx;
                            this->ny   = rhs.ny;
                            this->nz   = rhs.nz;
                            this->exr  = &rhs.exr[0];
                            this->exi  = &rhs.exi[0];
                            this->eyr  = &rhs.eyr[0];
                            this->eyi  = &rhs.eyi[0];
                            this->ezr  = &rhs.ezr[0];
                            this->ezi  = &rhs.ezi[0];
                            rhs.nx     = 0ULL;
                            rhs.ny     = 0ULL;
                            rhs.nz     = 0ULL;
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
                        this->exr  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->exi  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->eyr  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->eyi  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->ezr  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
                        this->ezi  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
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
                      std::size_t            nx;
                      std::size_t            ny;
                      std::size_t            nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif

                        inline H_r8_t() {
                         this->nx   = 0ULL;
                         this->ny   = 0ULL;
                         this->nz   = 0ULL;
                         this->mxr  = NULL;
                         this->mxi  = NULL;
                         this->myr  = NULL;
                         this->myi  = NULL;
                         this->mzr  = NULL;
                         this->mzi  = NULL;
                      }                    
                      
                      inline H_r8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline H_r8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->mxr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->mxi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->myr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->myi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->mzr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->mzi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mxr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->mxi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->myr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->myi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->mzr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->mzi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mxr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->mxi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->myr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->myi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->mzr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->mzi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline H_r8_t(const std::vector<double> &m_xr,
                                    const std::vector<double> &m_xi,
                                    const std::vector<double> &m_yr,
                                    const std::vector<double> &m_yi,
                                    const std::vector<double> &m_zr,
                                    const std::vector<double> &m_zi) {
                               
                               this->nx = m_xr.size();
                               this->ny = m_yr.size();
                               this->nz = m_zr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                               std::memcpy(this->mzr,&m_zr[0],lenz);
                               std::memcpy(this->mzi,&m_zi[0],lenz);     
                      }
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline H_r8_t(const std::valarray<double> &m_xr,
                                    const std::valarray<double> &m_xi,
                                    const std::valarray<double> &m_yr,
                                    const std::valarray<double> &m_yi,
                                    const std::valarray<double> &m_zr,
                                    const std::valarray<double> &m_zi) {
                               
                               this->nx = m_xr.size();
                               this->ny = m_yr.size();
                               this->nz = m_zr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                               std::memcpy(this->mzr,&m_zr[0],lenz);
                               std::memcpy(this->mzi,&m_zi[0],lenz);     
                      }
                      
                      
                      
                      inline H_r8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const double * __restrict m_xr,   
                                    const double * __restrict m_xi,
                                    const double * __restrict m_yr,
                                    const double * __restrict m_yi,
                                    const double * __restrict m_zr,
                                    const double * __restrict m_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&m_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->mxi[0],&m_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->myr[0],&m_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->myi[0],&m_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->mzr[0],&m_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->mzi[0],&m_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->mxr[0],&m_xr[0],this->nx);
	                  avx512_cached_memmove(&this->mxi[0],&m_xi[0],this->nx);
	                  avx512_cached_memmove(&this->myr[0],&m_yr[0],this->ny);
	                  avx512_cached_memmove(&this->myi[0],&m_yi[0],this->ny);
	                  avx512_cached_memmove(&this->mzr[0],&m_zr[0],this->nz);
	                  avx512_cached_memmove(&this->mzi[0],&m_zi[0],this->nz);
#endif          
                      }
                      
                      inline H_r8_t(H_r8_t &&rhs) {
                           
                          this->nx   = rhs.nx;
                          this->ny   = rhs.ny;
                          this->nz   = rhs.nz;
                          this->mxr  = &rhs.mxr[0];
                          this->mxi  = &rhs.mxi[0];
                          this->myr  = &rhs.myr[0];
                          this->myi  = &rhs.myi[0];
                          this->mzr  = &rhs.mzr[0];
                          this->mzi  = &rhs.mzi[0];
                          rhs.nx     = 0ULL;
                          rhs.ny     = 0ULL;
                          rhs.nz     = 0ULL;
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
                              gms_unmap(this->mxr,this->nx);
                              gms_unmap(this->mxi,this->nx);
                              gms_unmap(this->myr,this->ny); 
                              gms_unmap(this->myi,this->ny);
                              gms_unmap(this->mzr,this->nz);
                              gms_unmap(this->mzi,this->nz);
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
                            this->nx   = rhs.nx;
                            this->ny   = rhs.ny;
                            this->nz   = rhs.nz;
                            this->mxr  = &rhs.mxr[0];
                            this->mxi  = &rhs.mxi[0];
                            this->myr  = &rhs.myr[0];
                            this->myi  = &rhs.myi[0];
                            this->mzr  = &rhs.mzr[0];
                            this->mzi  = &rhs.mzi[0];
                            rhs.nx     = 0ULL;
                            rhs.ny     = 0ULL;
                            rhs.nz     = 0ULL;
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
                        this->mxr  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->mxi  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->myr  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->myi  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->mzr  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
                        this->mzi  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
                   }    
                   
              };


               struct __ATTR_ALIGN__(64) Je_c4_t {
                      // Complex electric current
                      
                      std::complex<float> * __restrict jex;
                      std::complex<float> * __restrict jey;
                      std::complex<float> * __restrict jez;
                      std::size_t                      nx;
                      std::size_t                      ny;
                      std::size_t                      nz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
   
                      inline Je_c4_t() noexcept(true) {
                          
                          this->nx   = 0ULL;
                          this->ny   = 0ULL;
                          this->jex  = NULL;
                          this->jey  = NULL;
                          this->jez  = NULL;
                      } 
                          
                     inline Je_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                     
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline Je_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jex = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->jey = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->jez = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jex = (std::complex<float>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->jey = (std::complex<float>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->jez = (std::complex<float>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jex = (std::complex<float>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->jey = (std::complex<float>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->jez = (std::complex<float>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   Je_c4_t(const std::vector<std::complex<float>> &j_ex,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<float>> &j_ey,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<float>> &j_ez) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_ex.size();
                          this->ny = j_ey.size();
                          this->nz = j_ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->jex,&j_ex[0],lenx);
                          std::memcpy(this->jey,&j_ey[0],leny);
                          std::memcpy(this->jez,&j_ez[0],lenz);       
                     }
                     
                   inline   Je_c4_t( const std::valarray<std::complex<float>> &j_ex,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<float>> &j_ey,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<float>> &j_ez) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_ex.size();
                          this->ny = j_ey.size();
                          this->nz = j_ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->jex,&j_ex[0],lenx);
                          std::memcpy(this->jey,&j_ey[0],leny);
                          std::memcpy(this->jez,&j_ez[0],lenz);       
                     }
                      
                   inline   Je_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const std::complex<float> * __restrict j_ex,
                                    const std::complex<float> * __restrict j_ey,
                                    const std::complex<float> * __restrict j_ez) {
                          
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jex[0],&j_ex[0],this->nx);
	                  avx512_uncached_memmove(&this->jey[0],&j_ey[0],this->ny);
	                  avx512_uncached_memmove(&this->jez[0],&j_ez[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jex[0],&j_ex[0],this->nx);
	                  avx512_cached_memmove(&this->jey[0],&j_ey[0],this->ny);
	                  avx512_cached_memmove(&this->jez[0],&j_ez[0],this->nz);
#endif
                   }  
                   
                  inline  Je_c4_t(Je_c4_t && rhs) {
                          
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jex   = &rhs.jex[0];
                          this->jey   = &rhs.jey[0];
                          this->jez   = &rhs.jez[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jex     = NULL;
                          rhs.jey     = NULL;
                          rhs.jez     = NULL;
                      }
                                 
                   Je_c4_t(const Je_c4_t &)             = delete;
                      
                   inline   ~Je_c4_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->jex,this->nx);
                             gms_unmap(this->jey,this->ny);
                             gms_unmap(this->jez,this->nz);
                          }
                          else {
                              gms_mm_free(this->jex);
                              gms_mm_free(this->jey);
                              gms_mm_free(this->jez);
                          }
                      }
                      
                    Je_c4_t & operator=(const Je_c4_t &) = delete;
                      
                    inline  Je_c4_t & operator=(Je_c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->jex);
                           gms_mm_free(this->jey);
                           gms_mm_free(this->jez);
                           this->nx    = rhs.nx;
                           this->ny    = rhs.ny;
                           this->nz    = rhs.nz;
                           this->jex   = &rhs.jex[0];
                           this->jey   = &rhs.jey[0];
                           this->jez   = &rhs.jez[0];
                           rhs.nx      = 0ULL;
                           rhs.ny      = 0ULL;
                           rhs.nz      = 0ULL;
                           rhs.jex     = NULL;
                           rhs.jey     = NULL;
                           rhs.jez     = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() {
                        using namespace gms::common;
                        this->jex  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<float>)*this->nx,64ULL);
                        this->jey  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<float>)*this->ny,64ULL);
                        this->jez  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<float>)*this->nz,64ULL);
                    }
                    
              };


               struct __ATTR_ALIGN__(64) Jm_c4_t {
                      // Complex magnetic current
                      std::complex<float> * __restrict jmx;
                      std::complex<float> * __restrict jmy;
                      std::complex<float> * __restrict jmz;
                      std::size_t                      nx;
                      std::size_t                      ny;
                      std::size_t                      nz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline Jm_c4_t() noexcept(true) {
                          
                          this->nx   = 0ULL;
                          this->ny   = 0ULL;
                          this->nz   = 0ULL;
                          this->jmx  = NULL;
                          this->jmy  = NULL;
                          this->jmz  = NULL;
                      } 
                          
                     inline Jm_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                     
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline Jm_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jmx = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->jmy = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->jmz = (std::complex<float>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jmx = (std::complex<float>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->jmy = (std::complex<float>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->jmz = (std::complex<float>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jmx = (std::complex<float>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<float>)*this->nx,prot,flags,fd,offset);
                                      this->jmy = (std::complex<float>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<float>)*this->ny,prot,flags,fd,offset);
                                      this->jmz = (std::complex<float>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<float>)*this->nz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   Jm_c4_t(const std::vector<std::complex<float>> &j_mx,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<float>> &j_my,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<float>> &j_mz) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_mx.size();
                          this->ny = j_my.size();
                          this->nz = j_mz.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->jmx,&j_mx[0],lenx);
                          std::memcpy(this->jmy,&j_my[0],leny);
                          std::memcpy(this->jmz,&j_mz[0],lenz);       
                     }
                     
                    inline   Jm_c4_t(const std::valarray<std::complex<float>> &j_mx,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<float>> &j_my,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<float>> &j_mz) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_mx.size();
                          this->ny = j_my.size();
                          this->nz = j_mz.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<float>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<float>)*this->nz;
                          std::memcpy(this->jmx,&j_mx[0],lenx);
                          std::memcpy(this->jmy,&j_my[0],leny);
                          std::memcpy(this->jmz,&j_mz[0],lenz);       
                     }
                      
                   inline   Jm_c4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const std::complex<float> * __restrict j_mx,
                                    const std::complex<float> * __restrict j_my,
                                    const std::complex<float> * __restrict j_mz) {
                          
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jmx[0],&j_mx[0],this->nx);
	                  avx512_uncached_memmove(&this->jmy[0],&j_my[0],this->ny);
	                  avx512_uncached_memmove(&this->jmz[0],&j_mz[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jmx[0],&j_mx[0],this->nx);
	                  avx512_cached_memmove(&this->jmy[0],&j_my[0],this->ny);
	                  avx512_cached_memmove(&this->jmz[0],&j_mz[0],this->nz);
#endif
                   }  
                   
                  inline  Jm_c4_t(Jm_c4_t && rhs) {
                          
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jmx   = &rhs.jmx[0];
                          this->jmy   = &rhs.jmy[0];
                          this->jmz   = &rhs.jmz[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jmx     = NULL;
                          rhs.jmy     = NULL;
                          rhs.jmz     = NULL;
                      }
                                 
                   Jm_c4_t(const Jm_c4_t &)             = delete;
                      
                   inline   ~Jm_c4_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->jmx,this->nx);
                             gms_unmap(this->jmy,this->ny);
                             gms_unmap(this->jmz,this->nz);
                          }
                          else {
                              gms_mm_free(this->jmx);
                              gms_mm_free(this->jmy);
                              gms_mm_free(this->jmz);
                          }
                      }
                      
                    Jm_c4_t & operator=(const Jm_c4_t &) = delete;
                      
                    inline  Jm_c4_t & operator=(Jm_c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->jmx);
                           gms_mm_free(this->jmy);
                           gms_mm_free(this->jmz);
                           this->nx    = rhs.nx;
                           this->ny    = rhs.ny;
                           this->nz    = rhs.nz;
                           this->jmx   = &rhs.jmx[0];
                           this->jmy   = &rhs.jmy[0];
                           this->jmz   = &rhs.jmz[0];
                           rhs.nx      = 0ULL;
                           rhs.ny      = 0ULL;
                           rhs.nz      = 0ULL;
                           rhs.jmx     = NULL;
                           rhs.jmy     = NULL;
                           rhs.jmz     = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() {
                        using namespace gms::common;
                        this->jmx  = (std::complex<float>*)
                                         gms_mm_malloc(sizeof(std::complex<float>)*this->nx,64ULL);
                        this->jmy  = (std::complex<float>*)
                                         gms_mm_malloc(sizeof(std::complex<float>)*this->ny,64ULL);
                        this->jmz  = (std::complex<float>*)
                                         gms_mm_malloc(sizeof(std::complex<float>)*this->nz,64ULL);
                    }
                    
              };


              typedef struct __ATTR_ALIGN__(64) Je_c8_t {
                      // Complex electric current
                     
                      std::complex<double> * __restrict jex;
                      std::complex<double> * __restrict jey;
                      std::complex<double> * __restrict jez;
                      std::size_t                       nx;
                      std::size_t                       ny;
                      std::size_t                       nz;
                      bool                              ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                      inline Je_c8_t() noexcept(true) {
                          
                          this->nx   = 0ULL;
                          this->ny   = 0ULL;
                          this->nz   = 0ULL;
                          this->jex  = NULL;
                          this->jey  = NULL;
                          this->jez  = NULL;
                      } 
                          
                     inline Je_c8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                     
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline Je_c8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jex = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->jey = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->jez = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jex = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->jey = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->jez = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jex = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->jey = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->jez = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   Je_c8_t(const std::vector<std::complex<double>> &j_ex,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<double>> &j_ey,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<double>> &j_ez) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_ex.size();
                          this->ny = j_ey.size();
                          this->nz = j_ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->jex,&j_ex[0],lenx);
                          std::memcpy(this->jey,&j_ey[0],leny);
                          std::memcpy(this->jez,&j_ez[0],lenz);       
                     }
                     
                   inline   Je_c8_t( const std::valarray<std::complex<double>> &j_ex,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<double>> &j_ey,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<double>> &j_ez) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_ex.size();
                          this->ny = j_ey.size();
                          this->nz = j_ez.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->jex,&j_ex[0],lenx);
                          std::memcpy(this->jey,&j_ey[0],leny);
                          std::memcpy(this->jez,&j_ez[0],lenz);       
                     }
                      
                   inline   Je_c8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const std::complex<double> * __restrict j_ex,
                                    const std::complex<double> * __restrict j_ey,
                                    const std::complex<double> * __restrict j_ez) {
                          
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jex[0],&j_ex[0],this->nx);
	                  avx512_uncached_memmove(&this->jey[0],&j_ey[0],this->ny);
	                  avx512_uncached_memmove(&this->jez[0],&j_ez[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jex[0],&j_ex[0],this->nx);
	                  avx512_cached_memmove(&this->jey[0],&j_ey[0],this->ny);
	                  avx512_cached_memmove(&this->jez[0],&j_ez[0],this->nz);
#endif
                   }  
                   
                  inline  Je_c8_t(Je_c8_t && rhs) {
                          
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jex   = &rhs.jex[0];
                          this->jey   = &rhs.jey[0];
                          this->jez   = &rhs.jez[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jex     = NULL;
                          rhs.jey     = NULL;
                          rhs.jez     = NULL;
                      }
                                 
                   Je_c8_t(const Je_c8_t &)             = delete;
                      
                   inline   ~Je_c8_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->jex,this->nx);
                             gms_unmap(this->jey,this->ny);
                             gms_unmap(this->jez,this->nz);
                          }
                          else {
                              gms_mm_free(this->jex);
                              gms_mm_free(this->jey);
                              gms_mm_free(this->jez);
                          }
                      }
                      
                    Je_c8_t & operator=(const Je_c8_t &) = delete;
                      
                    inline  Je_c8_t & operator=(Je_c8_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->jex);
                           gms_mm_free(this->jey);
                           gms_mm_free(this->jez);
                           this->nx    = rhs.nx;
                           this->ny    = rhs.ny;
                           this->nz    = rhs.nz;
                           this->jex   = &rhs.jex[0];
                           this->jey   = &rhs.jey[0];
                           this->jez   = &rhs.jez[0];
                           rhs.nx      = 0ULL;
                           rhs.ny      = 0ULL;
                           rhs.nz      = 0ULL;
                           rhs.jex     = NULL;
                           rhs.jey     = NULL;
                           rhs.jez     = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() {
                        using namespace gms::common;
                        this->jex  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nx,64ULL);
                        this->jey  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->ny,64ULL);
                        this->jez  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nz,64ULL);
                    }
                    
              };


              struct __ATTR_ALIGN__(64) Jm_c8_t {
                      // Complex electric current
                     
                      std::complex<double> * __restrict jmx;
                      std::complex<double> * __restrict jmy;
                      std::complex<double> * __restrict jmz;
                      std::size_t                       nx;
                      std::size_t                       ny;
                      std::size_t                       nz;
                      bool                              ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline Jm_c8_t() noexcept(true) {
                          
                          this->nx   = 0ULL;
                          this->ny   = 0ULL;
                          this->nz   = 0ULL;
                          this->jmx  = NULL;
                          this->jmy  = NULL;
                          this->jmz  = NULL;
                      } 
                          
                     inline Jm_c8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz) {
                     
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline Jm_c8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jmx = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->jmy = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->jmz = (std::complex<double>*)
                                                 gms_mmap_4KiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jmx = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->jmy = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->jmz = (std::complex<double>*)
                                                 gms_mmap_2MiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jmx = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nx,prot,flags,fd,offset);
                                      this->jmy = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->ny,prot,flags,fd,offset);
                                      this->jmz = (std::complex<double>*)
                                                 gms_mmap_1GiB(sizeof(std::complex<double>)*this->nz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                      
                    inline   Jm_c8_t(const std::vector<std::complex<double>> &j_mx,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<double>> &j_my,    //shall be of the same size (no error checking implemented)
                                     const std::vector<std::complex<double>> &j_mz) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_mx.size();
                          this->ny = j_my.size();
                          this->nz = j_mz.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->jmx,&j_mx[0],lex);
                          std::memcpy(this->jmy,&j_my[0],leny);
                          std::memcpy(this->jmz,&j_mz[0],lenz);       
                     }
                     
                    inline   Jm_c8_t(const std::valarray<std::complex<double>> &j_mx,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<double>> &j_my,    //shall be of the same size (no error checking implemented)
                                     const std::valarray<std::complex<double>> &j_mz) {  //shall be of the same size (no error checking implemented)
                         
                          this->nx = j_mx.size();
                          this->ny = j_my.size();
                          this->nz = j_mz.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->nx;
                          const std::size_t leny = sizeof(std::complex<double>)*this->ny;
                          const std::size_t lenz = sizeof(std::complex<double>)*this->nz;
                          std::memcpy(this->jmx,&j_mx[0],lex);
                          std::memcpy(this->jmy,&j_my[0],leny);
                          std::memcpy(this->jmz,&j_mz[0],lenz);       
                     }
                      
                   inline   Jm_c8_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const std::complex<double> * __restrict j_mx,
                                    const std::complex<double> * __restrict j_my,
                                    const std::complex<double> * __restrict j_mz) {
                          
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jmx[0],&j_mx[0],this->nx);
	                  avx512_uncached_memmove(&this->jmy[0],&j_my[0],this->ny);
	                  avx512_uncached_memmove(&this->jmz[0],&j_mz[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jmx[0],&j_mx[0],this->nx);
	                  avx512_cached_memmove(&this->jmy[0],&j_my[0],this->ny);
	                  avx512_cached_memmove(&this->jmz[0],&j_mz[0],this->nz);
#endif
                   }  
                   
                  inline  Jm_c8_t(Jm_c8_t && rhs) {
                          
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jmx   = &rhs.jmx[0];
                          this->jmy   = &rhs.jmy[0];
                          this->jmz   = &rhs.jmz[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jmx     = NULL;
                          rhs.jmy     = NULL;
                          rhs.jmz     = NULL;
                      }
                                 
                      Jm_c8_t(const Jm_c8_t &)             = delete;
                      
                   inline   ~Jm_c8_t() {
                      
                          using namespace gms::common;
                          if(this->ismmap) {
                             gms_unmap(this->jmx,this->nx);
                             gms_unmap(this->jmy,this->ny);
                             gms_unmap(this->jmz,this->nz);
                          }
                          else {
                              gms_mm_free(this->jmx);
                              gms_mm_free(this->jmy);
                              gms_mm_free(this->jmz);
                          }
                      }
                      
                      Jm_c8_t & operator=(const Jm_c8_t &) = delete;
                      
                    inline  Jm_c8_t & operator=(Jm_c8_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           gms_mm_free(this->jmx);
                           gms_mm_free(this->jmy);
                           gms_mm_free(this->jmz);
                           this->nx    = rhs.nx;
                           this->ny    = rhs.ny;
                           this->nz    = rhs.nz;
                           this->jmx   = &rhs.jmx[0];
                           this->jmy   = &rhs.jmy[0];
                           this->jmz   = &rhs.jmz[0];
                           rhs.nx      = 0ULL;
                           rhs.ny      = 0ULL;
                           rhs.nz      = 0ULL;
                           rhs.jmx     = NULL;
                           rhs.jmy     = NULL;
                           rhs.jmz     = NULL;
                           return (*this);
                      }
                      
                    inline void allocate() {
                        using namespace gms::common;
                        this->jmx  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nx,64ULL);
                        this->jmy  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->ny,64ULL);
                        this->jmz  = (std::complex<double>*)
                                         gms_mm_malloc(sizeof(std::complex<double>)*this->nz,64ULL);
                    }
                    
              };


               struct __ATTR_ALIGN__(64) Je_r4_t {
                     // ! Complex Electric Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                     
                      float * __restrict     jexr;
                      float * __restrict     jexi;
                      float * __restrict     jeyr;
                      float * __restrict     jeyi;
                      float * __restrict     jezr;
                      float * __restrict     jezi;
                      std::size_t            nx;
                      std::size_t            ny;
                      std::size_t            nz;
                      bool                   ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif 
                     inline Je_r4_t() {
                         this->nx    = 0ULL;
                         this->ny    = 0ULL;
                         this->nz    = 0ULL;
                         this->jexr  = NULL;
                         this->jexi  = NULL;
                         this->jeyr  = NULL;
                         this->jeyi  = NULL;
                         this->jezr  = NULL;
                         this->jezi  = NULL;
                      }                    
                      
                      inline Je_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline Je_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const int32_t prot,
                                     const int32_t flags,
                                     const int32_t fd,
                                     const int32_t offset,
                                     const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jexr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jexi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jeyr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jeyi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jezr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->jezi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jexr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jexi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jeyr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jeyi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jezr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->jezi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jexr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jexi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jeyr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jeyi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jezr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->jezi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline Je_r4_t(const std::vector<float> &je_xr,
                                     const std::vector<float> &je_xi,
                                     const std::vector<float> &je_yr,
                                     const std::vector<float> &je_yi,
                                     const std::vector<float> &je_zr,
                                     const std::vector<float> &je_zi) {
                               
                               this->nx = je_xr.size();
                               this->ny = je_yr.size();
                               this->nz = je_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenz = sizeof(float)*this->nz;
                               std::memcpy(this->jexr,&je_xr[0],lenx);
                               std::memcpy(this->jexi,&je_xi[0],lenx);
                               std::memcpy(this->jeyr,&je_yr[0],leny);
                               std::memcpy(this->jeyi,&je_yi[0],leny);
                               std::memcpy(this->jezr,&je_zr[0],lenz);
                               std::memcpy(this->jezi,&je_zi[0],lenz);     
                      }
                      
                      inline Je_r4_t(const std::valarray<float> &je_xr,
                                     const std::valarray<float> &je_xi,
                                     const std::valarray<float> &je_yr,
                                     const std::valarray<float> &je_yi,
                                     const std::valarray<float> &je_zr,
                                     const std::valarray<float> &je_zi) {
                               
                               this->nx = je_xr.size();
                               this->ny = je_yr.size();
                               this->nz = je_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenz = sizeof(float)*this->nz;
                               std::memcpy(this->jexr,&je_xr[0],lenx);
                               std::memcpy(this->jexi,&je_xi[0],lenx);
                               std::memcpy(this->jeyr,&je_yr[0],leny);
                               std::memcpy(this->jeyi,&je_yi[0],leny);
                               std::memcpy(this->jezr,&je_zr[0],lenz);
                               std::memcpy(this->jezi,&je_zi[0],lenz);     
                      }
                      
                      inline Je_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const float * __restrict je_xr,   
                                     const float * __restrict je_xi,
                                     const float * __restrict je_yr,
                                     const float * __restrict je_yi,
                                     const float * __restrict je_zr,
                                     const float * __restrict je_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jexr[0],&je_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->jexi[0],&je_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->jeyr[0],&je_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->jeyi[0],&je_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->jezr[0],&je_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->jezi[0],&je_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jexr[0],&je_xr[0],this->nx);
	                  avx512_cached_memmove(&this->jexi[0],&je_xi[0],this->nx);
	                  avx512_cached_memmove(&this->jeyr[0],&je_yr[0],this->ny);
	                  avx512_cached_memmove(&this->jeyi[0],&je_yi[0],this->ny);
	                  avx512_cached_memmove(&this->jezr[0],&je_zr[0],this->nz);
	                  avx512_cached_memmove(&this->jezi[0],&je_zi[0],this->nz);
#endif          
                      }
                      
                      inline Je_r4_t(Je_r4_t &&rhs) {
                           
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jexr  = &rhs.jexr[0];
                          this->jexi  = &rhs.jexi[0];
                          this->jeyr  = &rhs.jeyr[0];
                          this->jeyi  = &rhs.jeyi[0];
                          this->jezr  = &rhs.jezr[0];
                          this->jezi  = &rhs.jezi[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jexr    = NULL;
                          rhs.jexi    = NULL;
                          rhs.jeyr    = NULL;
                          rhs.jeyi    = NULL;
                          rhs.jezr    = NULL;
                          rhs.jezi    = NULL;
                      }   
                        
                      Je_r4_t(const Je_r4_t &)             = delete;
                      
                      inline ~Je_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->jexr,this->nx);
                              gms_unmap(this->jexi,this->nx);
                              gms_unmap(this->jeyr,this->ny); 
                              gms_unmap(this->jeyi,this->ny);
                              gms_unmap(this->jezr,this->nz);
                              gms_unmap(this->jezi,this->nz);
                           }
                           else {
                               gms_mm_free(this->jexr);
                               gms_mm_free(this->jexi);
                               gms_mm_free(this->jeyr);
                               gms_mm_free(this->jeyi);
                               gms_mm_free(this->jezr);
                               gms_mm_free(this->jezi);
                           }
                      }
                      
                      Je_r4_t & operator=(const Je_r4_t &) = delete;
                      
                      inline Je_r4_t & operator=(Je_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->jexr);
                            gms_mm_free(this->jexi);
                            gms_mm_free(this->jeyr);
                            gms_mm_free(this->jeyi);
                            gms_mm_free(this->jezr);
                            gms_mm_free(this->jezi);
                            this->nx    = rhs.nx;
                            this->ny    = rhs.ny;
                            this->nz    = rhs.nz;
                            this->jexr  = &rhs.jexr[0];
                            this->jexi  = &rhs.jexi[0];
                            this->jeyr  = &rhs.jeyr[0];
                            this->jeyi  = &rhs.jeyi[0];
                            this->jezr  = &rhs.jezr[0];
                            this->jezi  = &rhs.jezi[0];
                            rhs.nx      = 0ULL;
                            rhs.ny      = 0ULL;
                            rhs.nz      = 0ULL;
                            rhs.jexr    = NULL;
                            rhs.jexi    = NULL;
                            rhs.jeyr    = NULL;
                            rhs.jeyi    = NULL;
                            rhs.jezr    = NULL;
                            rhs.jezi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->jexr  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->jexi  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->jeyr  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->jeyi  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->jezr  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
                        this->jezi  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
                   }    
                   
              };


            struct __ATTR_ALIGN__(64) Jm_r4_t {
                     // ! Complex Magnetic Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      
                      float * __restrict jmxr;
                      float * __restrict jmxi;
                      float * __restrict jmyr;
                      float * __restrict jmyi;
                      float * __restrict jmzr;
                      float * __restrict jmzi;
                      std::size_t        nx;
                      std::size_t        ny;
                      std::size_t        nz;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif 
                     inline Jm_r4_t() {
                         this->nx    = 0ULL;
                         this->ny    = 0ULL;
                         this->nz    = 0ULL;
                         this->jmxr  = NULL;
                         this->jmxi  = NULL;
                         this->jmyr  = NULL;
                         this->jmyi  = NULL;
                         this->jmzr  = NULL;
                         this->jmzi  = NULL;
                      }                    
                      
                      inline Jm_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline Jm_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const int32_t prot,
                                     const int32_t flags,
                                     const int32_t fd,
                                     const int32_t offset,
                                     const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jmxr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jmxi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jmyr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jmyi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jmzr = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->jmzi = (float*)
                                                  gms_mmap_4KiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jmxr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jmxi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jmyr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jmyi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jmzr = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->jmzi = (float*)
                                                  gms_mmap_2MiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jmxr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jmxi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nx,prot,flags,fd,offset);
                                      this->jmyr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jmyi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->ny,prot,flags,fd,offset);
                                      this->jmzr = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->jmzi = (float*)
                                                  gms_mmap_1GiB(sizeof(float)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline Jm_r4_t(const std::vector<float> &jm_xr,
                                     const std::vector<float> &jm_xi,
                                     const std::vector<float> &jm_yr,
                                     const std::vector<float> &jm_yi,
                                     const std::vector<float> &jm_zr,
                                     const std::vector<float> &jm_zi) {
                               
                               this->nx = jm_xr.size();
                               this->ny = jm_yr.size();
                               this->nz = jm_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenz = sizeof(float)*this->nz;
                               std::memcpy(this->jmxr,&je_xr[0],lenx);
                               std::memcpy(this->jmxi,&je_xi[0],lenx);
                               std::memcpy(this->jmyr,&je_yr[0],leny);
                               std::memcpy(this->jmyi,&je_yi[0],leny);
                               std::memcpy(this->jmzr,&je_zr[0],lenz);
                               std::memcpy(this->jmzi,&je_zi[0],lenz);     
                      }
                      
                      inline Jm_r4_t(const std::valarray<float> &jm_xr,
                                     const std::valarray<float> &jm_xi,
                                     const std::valarray<float> &jm_yr,
                                     const std::valarray<float> &jm_yi,
                                     const std::valarray<float> &jm_zr,
                                     const std::valarray<float> &jm_zi) {
                               
                               this->nx = jm_xr.size();
                               this->ny = jm_yr.size();
                               this->nz = jm_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->nx;
                               const std::size_t leny = sizeof(float)*this->ny;
                               const std::size_t lenz = sizeof(float)*this->nz;
                               std::memcpy(this->jmxr,&jm_xr[0],lenx);
                               std::memcpy(this->jmxi,&jm_xi[0],lenx);
                               std::memcpy(this->jmyr,&jm_yr[0],leny);
                               std::memcpy(this->jmyi,&jm_yi[0],leny);
                               std::memcpy(this->jmzr,&jm_zr[0],lenz);
                               std::memcpy(this->jmzi,&jm_zi[0],lenz);     
                      }
                      
                      inline Jm_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const float * __restrict jm_xr,   
                                     const float * __restrict jm_xi,
                                     const float * __restrict jm_yr,
                                     const float * __restrict jm_yi,
                                     const float * __restrict jm_zr,
                                     const float * __restrict jm_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jmxr[0],&jm_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->jmxi[0],&jm_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->jmyr[0],&jm_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->jmyi[0],&jm_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->jmzr[0],&jm_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->jmzi[0],&jm_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jmxr[0],&jm_xr[0],this->nx);
	                  avx512_cached_memmove(&this->jmxi[0],&jm_xi[0],this->nx);
	                  avx512_cached_memmove(&this->jmyr[0],&jm_yr[0],this->ny);
	                  avx512_cached_memmove(&this->jmyi[0],&jm_yi[0],this->ny);
	                  avx512_cached_memmove(&this->jmzr[0],&jm_zr[0],this->nz);
	                  avx512_cached_memmove(&this->jmzi[0],&jm_zi[0],this->nz);
#endif          
                      }
                      
                      inline Jm_r4_t(Jm_r4_t &&rhs) {
                           
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jmxr  = &rhs.jmxr[0];
                          this->jmxi  = &rhs.jmxi[0];
                          this->jmyr  = &rhs.jmyr[0];
                          this->jmyi  = &rhs.jmyi[0];
                          this->jmzr  = &rhs.jmzr[0];
                          this->jmzi  = &rhs.jmzi[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jmxr    = NULL;
                          rhs.jmxi    = NULL;
                          rhs.jmyr    = NULL;
                          rhs.jmyi    = NULL;
                          rhs.jmzr    = NULL;
                          rhs.jmzi    = NULL;
                      }   
                        
                      Jm_r4_t(const Jm_r4_t &)             = delete;
                      
                      inline ~Jm_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->jmxr,this->nx);
                              gms_unmap(this->jmxi,this->nx);
                              gms_unmap(this->jmyr,this->ny); 
                              gms_unmap(this->jmyi,this->ny);
                              gms_unmap(this->jmzr,this->nz);
                              gms_unmap(this->jmzi,this->nz);
                           }
                           else {
                               gms_mm_free(this->jmxr);
                               gms_mm_free(this->jmxi);
                               gms_mm_free(this->jmyr);
                               gms_mm_free(this->jmyi);
                               gms_mm_free(this->jmzr);
                               gms_mm_free(this->jmzi);
                           }
                      }
                      
                      Jm_r4_t & operator=(const Jm_r4_t &) = delete;
                      
                      inline Jm_r4_t & operator=(Jm_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->jmxr);
                            gms_mm_free(this->jmxi);
                            gms_mm_free(this->jmyr);
                            gms_mm_free(this->jmyi);
                            gms_mm_free(this->jmzr);
                            gms_mm_free(this->jmzi);
                            this->nx    = rhs.nx;
                            this->ny    = rhs.ny;
                            this->nz    = rhs.nz;
                            this->jmxr  = &rhs.jmxr[0];
                            this->jmxi  = &rhs.jmxi[0];
                            this->jmyr  = &rhs.jmyr[0];
                            this->jmyi  = &rhs.jmyi[0];
                            this->jmzr  = &rhs.jmzr[0];
                            this->jmzi  = &rhs.jmzi[0];
                            rhs.nx      = 0ULL;
                            rhs.ny      = 0ULL;
                            rhs.nz      = 0ULL;
                            rhs.jmxr    = NULL;
                            rhs.jmxi    = NULL;
                            rhs.jmyr    = NULL;
                            rhs.jmyi    = NULL;
                            rhs.jmzr    = NULL;
                            rhs.jmzi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->jmxr  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->jmxi  = (float*)gms_mm_malloc(sizeof(float)*this->nx,64ULL);
                        this->jmyr  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->jmyi  = (float*)gms_mm_malloc(sizeof(float)*this->ny,64ULL);
                        this->jmzr  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
                        this->jmzi  = (float*)gms_mm_malloc(sizeof(float)*this->nz,64ULL);
                   }    
                    
              };


              typedef struct __ATTR_ALIGN__(64) Je_r8_t {
                     // ! Complex Electric Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      
                      double * __restrict jexr;
                      double * __restrict jexi;
                      double * __restrict jeyr;
                      double * __restrict jeyi;
                      double * __restrict jezr;
                      double * __restrict jezi;
                      std::size_t         nx;
                      std::size_t         ny;
                      std::size_t         nz;
                      bool                ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif 
                     inline Je_r8_t() {
                         this->nx    = 0ULL;
                         this->ny    = 0ULL;
                         this->nz    = 0ULL;
                         this->jexr  = NULL;
                         this->jexi  = NULL;
                         this->jeyr  = NULL;
                         this->jeyi  = NULL;
                         this->jezr  = NULL;
                         this->jezi  = NULL;
                      }                    
                      
                      inline Je_r8_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline J8_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const int32_t prot,
                                     const int32_t flags,
                                     const int32_t fd,
                                     const int32_t offset,
                                     const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jexr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jexi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jeyr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jeyi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jezr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->jezi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jexr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jexi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jeyr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jeyi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jezr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->jezi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jexr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jexi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jeyr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jeyi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jezr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->jezi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline Je_r4_t(const std::vector<double> &je_xr,
                                     const std::vector<double> &je_xi,
                                     const std::vector<double> &je_yr,
                                     const std::vector<double> &je_yi,
                                     const std::vector<double> &je_zr,
                                     const std::vector<double> &je_zi) {
                               
                               this->nx = je_xr.size();
                               this->ny = je_yr.size();
                               this->nz = je_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->jexr,&je_xr[0],lenx);
                               std::memcpy(this->jexi,&je_xi[0],lenx);
                               std::memcpy(this->jeyr,&je_yr[0],leny);
                               std::memcpy(this->jeyi,&je_yi[0],leny);
                               std::memcpy(this->jezr,&je_zr[0],lenz);
                               std::memcpy(this->jezi,&je_zi[0],lenz);     
                      }
                      
                      inline Je_r4_t(const std::valarray<double> &je_xr,
                                     const std::valarray<double> &je_xi,
                                     const std::valarray<double> &je_yr,
                                     const std::valarray<double> &je_yi,
                                     const std::valarray<double> &je_zr,
                                     const std::valarray<double> &je_zi) {
                               
                               this->nx = je_xr.size();
                               this->ny = je_yr.size();
                               this->nz = je_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->jexr,&je_xr[0],lenx);
                               std::memcpy(this->jexi,&je_xi[0],lenx);
                               std::memcpy(this->jeyr,&je_yr[0],leny);
                               std::memcpy(this->jeyi,&je_yi[0],leny);
                               std::memcpy(this->jezr,&je_zr[0],lenz);
                               std::memcpy(this->jezi,&je_zi[0],lenz);     
                      }
                      
                      inline Je_r4_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const double * __restrict je_xr,   
                                     const double * __restrict je_xi,
                                     const double * __restrict je_yr,
                                     const double * __restrict je_yi,
                                     const double * __restrict je_zr,
                                     const double * __restrict je_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jexr[0],&je_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->jexi[0],&je_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->jeyr[0],&je_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->jeyi[0],&je_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->jezr[0],&je_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->jezi[0],&je_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jexr[0],&je_xr[0],this->nx);
	                  avx512_cached_memmove(&this->jexi[0],&je_xi[0],this->nx);
	                  avx512_cached_memmove(&this->jeyr[0],&je_yr[0],this->ny);
	                  avx512_cached_memmove(&this->jeyi[0],&je_yi[0],this->ny);
	                  avx512_cached_memmove(&this->jezr[0],&je_zr[0],this->nz);
	                  avx512_cached_memmove(&this->jezi[0],&je_zi[0],this->nz);
#endif          
                      }
                      
                      inline Je_r4_t(Je_r4_t &&rhs) {
                           
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jexr  = &rhs.jexr[0];
                          this->jexi  = &rhs.jexi[0];
                          this->jeyr  = &rhs.jeyr[0];
                          this->jeyi  = &rhs.jeyi[0];
                          this->jezr  = &rhs.jezr[0];
                          this->jezi  = &rhs.jezi[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jexr    = NULL;
                          rhs.jexi    = NULL;
                          rhs.jeyr    = NULL;
                          rhs.jeyi    = NULL;
                          rhs.jezr    = NULL;
                          rhs.jezi    = NULL;
                      }   
                        
                      Je_r4_t(const Je_r4_t &)             = delete;
                      
                      inline ~Je_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->jexr,this->nx);
                              gms_unmap(this->jexi,this->nx);
                              gms_unmap(this->jeyr,this->ny); 
                              gms_unmap(this->jeyi,this->ny);
                              gms_unmap(this->jezr,this->nz);
                              gms_unmap(this->jezi,this->nz);
                           }
                           else {
                               gms_mm_free(this->jexr);
                               gms_mm_free(this->jexi);
                               gms_mm_free(this->jeyr);
                               gms_mm_free(this->jeyi);
                               gms_mm_free(this->jezr);
                               gms_mm_free(this->jezi);
                           }
                      }
                      
                      Je_r4_t & operator=(const Je_r4_t &) = delete;
                      
                      inline Je_r4_t & operator=(Je_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->jexr);
                            gms_mm_free(this->jexi);
                            gms_mm_free(this->jeyr);
                            gms_mm_free(this->jeyi);
                            gms_mm_free(this->jezr);
                            gms_mm_free(this->jezi);
                            this->nx    = rhs.nx;
                            this->ny    = rhs.ny;
                            this->nz    = rhs.nz;
                            this->jexr  = &rhs.jexr[0];
                            this->jexi  = &rhs.jexi[0];
                            this->jeyr  = &rhs.jeyr[0];
                            this->jeyi  = &rhs.jeyi[0];
                            this->jezr  = &rhs.jezr[0];
                            this->jezi  = &rhs.jezi[0];
                            rhs.nx      = 0ULL;
                            rhs.ny      = 0ULL;
                            rhs.nz      = 0ULL;
                            rhs.jexr    = NULL;
                            rhs.jexi    = NULL;
                            rhs.jeyr    = NULL;
                            rhs.jeyi    = NULL;
                            rhs.jezr    = NULL;
                            rhs.jezi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->jexr  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->jexi  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->jeyr  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->jeyi  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->jezr  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
                        this->jezi  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
                   }    
              };


              typedef struct __ATTR_ALIGN__(64) Jm_r8_t {
                     // ! Complex Magnetic Current  decomposed into real and imaginary parts 
                     // ! To be used mainly by the integrators.
                      
                      double * __restrict jm_xr;
                      double * __restrict jm_xi;
                      double * __restrict jm_yr;
                      double * __restrict jm_yi;
                      double * __restrict jm_zr;
                      double * __restrict jm_zi;
                      std::size_t         nx;
                      std::size_t         ny;
                      std::size_t         nz;
                      bool                ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif 
                     inline Jm_r8_t() {
                         this->nx    = 0ULL;
                         this->ny    = 0ULL;
                         this->nz    = 0ULL;
                         this->jmxr  = NULL;
                         this->jmxi  = NULL;
                         this->jmyr  = NULL;
                         this->jmyi  = NULL;
                         this->jmzr  = NULL;
                         this->jmzi  = NULL;
                      }                    
                      
                      inline Jm_r8_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz) {
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline Jm_r8_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const int32_t prot,
                                     const int32_t flags,
                                     const int32_t fd,
                                     const int32_t offset,
                                     const int32_t fsize) {
                             using namespace gms::common;
                             this->nx = _nx;
                             this->ny = _ny;
                             this->nz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->jmxr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jmxi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jmyr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jmyi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jmzr = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->jmzi = (double*)
                                                  gms_mmap_4KiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->jmxr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jmxi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jmyr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jmyi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jmzr = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->jmzi = (double*)
                                                  gms_mmap_2MiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->jmxr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jmxi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nx,prot,flags,fd,offset);
                                      this->jmyr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jmyi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->ny,prot,flags,fd,offset);
                                      this->jmzr = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->jmzi = (double*)
                                                  gms_mmap_1GiB(sizeof(double)*this->nz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline Jm_r8_t(const std::vector<double> &jm_xr,
                                     const std::vector<double> &jm_xi,
                                     const std::vector<double> &jm_yr,
                                     const std::vector<double> &jm_yi,
                                     const std::vector<double> &jm_zr,
                                     const std::vector<double> &jm_zi) {
                               
                               this->nx = jm_xr.size();
                               this->ny = jm_yr.size();
                               this->nz = jm_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->jmxr,&je_xr[0],lenx);
                               std::memcpy(this->jmxi,&je_xi[0],lenx);
                               std::memcpy(this->jmyr,&je_yr[0],leny);
                               std::memcpy(this->jmyi,&je_yi[0],leny);
                               std::memcpy(this->jmzr,&je_zr[0],lenz);
                               std::memcpy(this->jmzi,&je_zi[0],lenz);     
                      }
                      
                      inline Jm_r8_t(const std::valarray<double> &jm_xr,
                                     const std::valarray<double> &jm_xi,
                                     const std::valarray<double> &jm_yr,
                                     const std::valarray<double> &jm_yi,
                                     const std::valarray<double> &jm_zr,
                                     const std::valarray<double> &jm_zi) {
                               
                               this->nx = jm_xr.size();
                               this->ny = jm_yr.size();
                               this->nz = jm_zr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(double)*this->nx;
                               const std::size_t leny = sizeof(double)*this->ny;
                               const std::size_t lenz = sizeof(double)*this->nz;
                               std::memcpy(this->jmxr,&jm_xr[0],lenx);
                               std::memcpy(this->jmxi,&jm_xi[0],lenx);
                               std::memcpy(this->jmyr,&jm_yr[0],leny);
                               std::memcpy(this->jmyi,&jm_yi[0],leny);
                               std::memcpy(this->jmzr,&jm_zr[0],lenz);
                               std::memcpy(this->jmzi,&jm_zi[0],lenz);     
                      }
                      
                      inline Jm_r8_t(const std::size_t _nx,
                                     const std::size_t _ny,
                                     const std::size_t _nz,
                                     const double * __restrict jm_xr,   
                                     const double * __restrict jm_xi,
                                     const double * __restrict jm_yr,
                                     const double * __restrict jm_yi,
                                     const double * __restrict jm_zr,
                                     const double * __restrict jm_zi) {
                         
                          this->nx = _nx;
                          this->ny = _ny;
                          this->nz = _nz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->jmxr[0],&jm_xr[0],this->nx);
	                  avx512_uncached_memmove(&this->jmxi[0],&jm_xi[0],this->nx);
	                  avx512_uncached_memmove(&this->jmyr[0],&jm_yr[0],this->ny);
	                  avx512_uncached_memmove(&this->jmyi[0],&jm_yi[0],this->ny);
	                  avx512_uncached_memmove(&this->jmzr[0],&jm_zr[0],this->nz);
	                  avx512_uncached_memmove(&this->jmzi[0],&jm_zi[0],this->nz);
#else
	                  avx512_cached_memmove(&this->jmxr[0],&jm_xr[0],this->nx);
	                  avx512_cached_memmove(&this->jmxi[0],&jm_xi[0],this->nx);
	                  avx512_cached_memmove(&this->jmyr[0],&jm_yr[0],this->ny);
	                  avx512_cached_memmove(&this->jmyi[0],&jm_yi[0],this->ny);
	                  avx512_cached_memmove(&this->jmzr[0],&jm_zr[0],this->nz);
	                  avx512_cached_memmove(&this->jmzi[0],&jm_zi[0],this->nz);
#endif          
                      }
                      
                      inline Jm_r8_t(Jm_r8_t &&rhs) {
                           
                          this->nx    = rhs.nx;
                          this->ny    = rhs.ny;
                          this->nz    = rhs.nz;
                          this->jmxr  = &rhs.jmxr[0];
                          this->jmxi  = &rhs.jmxi[0];
                          this->jmyr  = &rhs.jmyr[0];
                          this->jmyi  = &rhs.jmyi[0];
                          this->jmzr  = &rhs.jmzr[0];
                          this->jmzi  = &rhs.jmzi[0];
                          rhs.nx      = 0ULL;
                          rhs.ny      = 0ULL;
                          rhs.nz      = 0ULL;
                          rhs.jmxr    = NULL;
                          rhs.jmxi    = NULL;
                          rhs.jmyr    = NULL;
                          rhs.jmyi    = NULL;
                          rhs.jmzr    = NULL;
                          rhs.jmzi    = NULL;
                      }   
                        
                      Jm_r8_t(const Jm_r8_t &)             = delete;
                      
                      inline ~Jm_r8_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap(this->jmxr,this->nx);
                              gms_unmap(this->jmxi,this->nx);
                              gms_unmap(this->jmyr,this->ny); 
                              gms_unmap(this->jmyi,this->ny);
                              gms_unmap(this->jmzr,this->nz);
                              gms_unmap(this->jmzi,this->nz);
                           }
                           else {
                               gms_mm_free(this->jmxr);
                               gms_mm_free(this->jmxi);
                               gms_mm_free(this->jmyr);
                               gms_mm_free(this->jmyi);
                               gms_mm_free(this->jmzr);
                               gms_mm_free(this->jmzi);
                           }
                      }
                      
                      Jm_r8_t & operator=(const Jm_r8_t &) = delete;
                      
                      inline Jm_r8_t & operator=(Jm_r8_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            gms_mm_free(this->jmxr);
                            gms_mm_free(this->jmxi);
                            gms_mm_free(this->jmyr);
                            gms_mm_free(this->jmyi);
                            gms_mm_free(this->jmzr);
                            gms_mm_free(this->jmzi);
                            this->nx    = rhs.nx;
                            this->ny    = rhs.ny;
                            this->nz    = rhs.nz;
                            this->jmxr  = &rhs.jmxr[0];
                            this->jmxi  = &rhs.jmxi[0];
                            this->jmyr  = &rhs.jmyr[0];
                            this->jmyi  = &rhs.jmyi[0];
                            this->jmzr  = &rhs.jmzr[0];
                            this->jmzi  = &rhs.jmzi[0];
                            rhs.nx      = 0ULL;
                            rhs.ny      = 0ULL;
                            rhs.nz      = 0ULL;
                            rhs.jmxr    = NULL;
                            rhs.jmxi    = NULL;
                            rhs.jmyr    = NULL;
                            rhs.jmyi    = NULL;
                            rhs.jmzr    = NULL;
                            rhs.jmzi    = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->jmxr  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->jmxi  = (double*)gms_mm_malloc(sizeof(double)*this->nx,64ULL);
                        this->jmyr  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->jmyi  = (double*)gms_mm_malloc(sizeof(double)*this->ny,64ULL);
                        this->jmzr  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
                        this->jmzi  = (double*)gms_mm_malloc(sizeof(double)*this->nz,64ULL);
                   }    
              };


            







       } // radiolocation

} // gms














#endif /*__GMS_ANTENNA_MODEL_COMMON_H__*/
