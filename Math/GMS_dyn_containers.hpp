
#ifndef __GMS_DYN_CONTAINERS_HPP__
#define __GMS_DYN_CONTAINERS_HPP__ 051220230457


namespace file_info {

     const unsigned int GMS_DYN_CONTAINERS_MAJOR = 1;
     const unsigned int GMS_DYN_CONTAINERS_MINOR = 1;
     const unsigned int GMS_DYN_CONTAINERS_MICRO = 0;
     const unsigned int GMS_DYN_CONTAINERS_FULLVER =
       1000U*GMS_DYN_CONTAINERS_MAJOR+100U*GMS_DYN_CONTAINERS_MINOR+
       10U*GMS_DYN_CONTAINERS_MICRO;
     const char * const GMS_DYN_CONTAINERS_CREATION_DATE = "05-12-2023 04:57 +00200 (TUE 05 DEC 2023 GMT+2)";
     const char * const GMS_DYN_CONTAINERS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_DYN_CONTAINERS_SYNOPSIS      = "Dynamically allocated, 64-byte aligned dense containers."

}





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


           struct __ATTR_ALIGN__(64) DC3D_c4_t {
                      // Complmx electric field.
                      std::complex<float> * __restrict mx
                      std::complex<float> * __restrict my;
                      std::complex<float> * __restrict mz;
                      std::size_t                      mmnx;
                      std::size_t                      mny;
                      std::size_t                      mnz;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                     inline DC3D_c4_t() nomxcept(true) {
                          
                          this->mnx  = 0ULL;
                          this->mny  = 0ULL;
                          this->mnz  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                          this->mz  = NULL;
                      } 
                          
                     inline DC3D_c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::size_t _mnz) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC3D_c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::size_t _mnz,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             this->mny = _mny;
                             this->mnz = _mnz;
                             switch (fsize) {
                                 case:0
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->mz = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnz,prot,flags,fd,offset); 
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline DC3D_c4_t(const std::vector<std::complex<float>> &m_x,    //shall be of the same size (no error checking implemented)
                                   const std::vector<std::complex<float>> &m_y,    //shall be of the same size (no error checking implemented)
                                   const std::vector<std::complex<float>> &m_z) {  //shall be of the same size (no error checking implemented)
                          using namespace gms::common;
                          
                          this->mnx = m_x.size();
                          this->mny = m_y.size();
                          this->mnz = m_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          const std::size_t lemnz = sizeof(std::complex<float>)*this->mnz;
                          std::memcpy(this->mx,&m_x[0],lemnx);
                          std::memcpy(this->my,&m_y[0],lemny);
                          std::memcpy(this->mz,&m_z[0],lemnz);       
                     }
                     
                     inline DC3D_c4_t(const std::valarray<std::complex<float>> &m_x,
                                   const std::valarray<std::complex<float>> &m_y,
                                   const std::valarray<std::complex<float>> &m_z) {
                          using namespace gms::common;
                          
                          this->mnx = m_x.size();
                          this->mny = m_y.size();
                          this->mnz = m_z.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          const std::size_t lemnz = sizeof(std::complex<float>)*this->mnz;
                          std::memcpy(this->mx,&e_x[0],lemnx);
                          std::memcpy(this->my,&e_y[0],lemny);
                          std::memcpy(this->mz,&e_z[0],lemnz);           
                     }
                             
                             
                      
                     inline DC3D_c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::size_t _mnz,
                                   const std::complex<float> * __restrict m_x,
                                   const std::complex<float> * __restrict m_y,
                                   const std::complex<float> * __restrict m_z) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->mnz);
#endif
                   }  
                   
                    inline  DC3D_c4_t(DC3D_c4_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mnz = rhs.mnz;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          this->mz   = &rhs.mz[0];
                          rhs.mnx    = 0ULL;
                          rhs.mny    = 0ULL;
                          rhs.mnz    = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                          rhs.mz     = NULL;
                      }
                                 
                     DC3D_c4_t(const DC3D_c4_t &)             = delete;
                      
                     inline ~DC3D_c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap<std::complex<float>>(this->mx,this->mnx);
                             gms_unmap<std::complex<float>>(this->my,this->mny);
                             gms_unmap<std::complex<float>>(this->mz,this->mnz);
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                              gms_mm_free(this->mz);
                          }
                      }
                      
                     DC3D_c4_t & operator=(const DC3D_c4_t &) = delete;
                      
                     inline DC3D_c4_t & operator=(DC3D_c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) {
                             gms_unmap<std::complex<float>>(this->mx,this->mnx);
                             gms_unmap<std::complex<float>>(this->my,this->mny);
                             gms_unmap<std::complex<float>>(this->mz,this->mnz);
                           }
                           else {
                             gms_mm_free(this->mx);
                             gms_mm_free(this->my);
                             gms_mm_free(this->mz);
                           }
                           this->mnx   = rhs.mnx;
                           this->mny   = rhs.mny;
                           this->mnz   = rhs.mnz;
                           this->mx    = &rhs.mx[0];
                           this->my    = &rhs.my[0];
                           this->mz    = &rhs.mz[0];
                           rhs.mnx     = 0ULL;
                           rhs.mny     = 0ULL;
                           rhs.mnz     = 0ULL;
                           rhs.mx      = NULL;
                           rhs.my      = NULL;
                           rhs.mz      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                          this->my  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mny,64ULL);
                          this->mz  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mny,64ULL);
                      }
                      
                     
                    
              };
              
              
               struct __ATTR_ALIGN__(64) DC2D_c4_t {
                      // Complmx electric field.
                      std::complex<float> * __restrict mx
                      std::complex<float> * __restrict my;
                      std::size_t                      mmnx;
                      std::size_t                      mny;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,31)
#endif
                     inline DC2D_c4_t() nomxcept(true) {
                          
                          this->mnx  = 0ULL;
                          this->mny  = 0ULL;
                          this->mx  = NULL;
                          this->my  = NULL;
                         
                      } 
                          
                     inline DC2D_c4_t(const std::size_t _mnx,
                                   const std::size_t _mny) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC2D_c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             this->mny = _mny;
                             switch (fsize) {
                                 case:0
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->my = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline DC2D_c4_t(const std::vector<std::complex<float>> &m_x,    //shall be of the same size (no error checking implemented)
                                   const std::vector<std::complex<float>> &m_y) {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = m_x.size();
                          this->mny = m_y.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          std::memcpy(this->mx,&m_x[0],lemnx);
                          std::memcpy(this->my,&m_y[0],lemny);
                        
                     }
                     
                     inline DC2D_c4_t(const std::valarray<std::complex<float>> &m_x,
                                   const std::valarray<std::complex<float>> &m_y) {
                                 
                          using namespace gms::common;
                          this->mnx = m_x.size();
                          this->mny = m_y.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          std::memcpy(this->mx,&e_x[0],lemnx);
                          std::memcpy(this->my,&e_y[0],lemny);
                             
                     }
                             
                             
                      
                     inline DC2D_c4_t(const std::size_t _mnx,
                                   const std::size_t _mny,
                                   const std::complex<float> * __restrict m_x,
                                   const std::complex<float> * __restrict m_y) {
                                
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                 
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  
#endif
                   }  
                   
                    inline  DC2D_c4_t(DC2D_c4_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mx   = &rhs.mx[0];
                          this->my   = &rhs.my[0];
                          rhs.mnx    = 0ULL;
                          rhs.mny    = 0ULL;
                          rhs.mx     = NULL;
                          rhs.my     = NULL;
                         
                      }
                                 
                     DC2D_c4_t(const DC2D_c4_t &)             = delete;
                      
                     inline ~DC2D_c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) {
                             gms_unmap<std::complex<float>>(this->mx,this->mnx);
                             gms_unmap<std::complex<float>>(this->my,this->mny);
                             
                          }
                          else {
                              gms_mm_free(this->mx);
                              gms_mm_free(this->my);
                             
                          }
                      }
                      
                     DC2D_c4_t & operator=(const DC2D_c4_t &) = delete;
                      
                     inline DC2D_c4_t & operator=(DC2D_c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) {
                             gms_unmap<std::complex<float>>(this->mx,this->mnx);
                             gms_unmap<std::complex<float>>(this->my,this->mny);
                            
                           }
                           else {
                             gms_mm_free(this->mx);
                             gms_mm_free(this->my);
                             
                           }
                           this->mnx   = rhs.mnx;
                           this->mny   = rhs.mny;
                           this->mx    = &rhs.mx[0];
                           this->my    = &rhs.my[0];
                           rhs.mnx     = 0ULL;
                           rhs.mny     = 0ULL;
                           rhs.mx      = NULL;
                           rhs.my      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                          this->my  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mny,64ULL);
                     }
                      
               };
               
               
               struct __ATTR_ALIGN__(64) DC1D_c4_t {
                      // Complmx electric field.
                      std::complex<float> * __restrict mx
                      std::size_t                      mmnx;
                      bool                             ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,47)
#endif
                     inline DC1D_c4_t() nomxcept(true) {
                          
                          this->mnx  = 0ULL;
                          this->mx  = NULL;
                                      
                      } 
                          
                     inline DC1D_c4_t(const std::size_t _mnx) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline DC1D_c4_t(const std::size_t _mnx,
                                      const int32_t prot,
                                   const int32_t flags,
                                   const int32_t fd,
                                   const int32_t offset,
                                   const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _mnx;
                             switch (fsize) {
                                 case:0
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mx = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline DC1D_c4_t(const std::vector<std::complex<float>> &m_x) {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = m_x.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          std::memcpy(this->mx,&m_x[0],lemnx);
                                                
                     }
                     
                     inline DC1D_c4_t(const std::valarray<std::complex<float>> &m_x) {
                                   
                          using namespace gms::common;
                          this->mnx = m_x.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          std::memcpy(this->mx,&e_x[0],lemnx);
                                                    
                     }
                             
                             
                      
                     inline DC1D_c4_t(const std::size_t _mnx,
                                      const std::complex<float> * __restrict m_x) {
                                 
                          using namespace gms::common;
                          this->mnx = _mnx;
                          allocate();
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                                    
#endif
                   }  
                   
                    inline  DC1D_c4_t(DC1D_c4_t && rhs) {
                          
                          this->mnx = rhs.mnx;
                          this->mx   = &rhs.mx[0];
                          rhs.mnx    = 0ULL;
                          rhs.mx     = NULL;
                                                
                      }
                                 
                     DC1D_c4_t(const DC1D_c4_t &)             = delete;
                      
                     inline ~DC1D_c4_t() {
                      
                          using namespace gms::common;
                          if(this->issmap) 
                             gms_unmap<std::complex<float>>(this->mx,this->mnx);
                          else 
                             gms_mm_free(this->mx);    
                     }    
                         
                                  
                                                    
                     DC1D_c4_t & operator=(const DC1D_c4_t &) = delete;
                      
                     inline DC1D_c4_t & operator=(DC1D_c4_t &&rhs) {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                           if(this->ismmap) 
                              gms_unmap<std::complex<float>>(this->mx,this->mnx);
                           else 
                              gms_mm_free(this->mx);
                                                                                      
                           this->mnx   = rhs.mnx;
                           this->mx    = &rhs.mx[0];
                           rhs.mnx     = 0ULL;
                           rhs.mx      = NULL;
                           return (*this);
                      }
                      
                      inline void allocate() {
                          using namespace gms::common;
                          this->mx  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                       
                     }
                      
               };
               
               
               struct __ATTR_ALIGN__(64) DC3D_r4_t {
                     
                      float * __restrict mxr;
                      float * __restrict mxi;
                      float * __restrict myr;
                      float * __restrict myi;
                      float * __restrict mzr;
                      float * __restrict mzi;
                      std::size_t        mnx;
                      std::size_t        mny;
                      std::size_t        mnz;
                      bool              ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,51)
#endif
                      inline DC3D_r4_t() {
                         this->mnx   = 0ULL;
                         this->mny   = 0ULL;
                         this->mnz   = 0ULL;
                         this->mxr  = NULL;
                         this->mxi  = NULL;
                         this->myr  = NULL;
                         this->myi  = NULL;
                         this->mzr  = NULL;
                         this->mzi  = NULL;
                      }                    
                      
                      inline DC3D_r4_t(const std::size_t _nx,
                                       const std::size_t _ny,
                                      const std::size_t _nz) {
                             this->mnx = _nx;
                             this->mny = _ny;
                             this->mnz = _nz;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline DC3D_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const std::size_t _nz,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _nx;
                             this->mny = _ny;
                             this->mnz = _nz;
                             switch (fsize) {
                                 case:0
                                      this->mxr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mxr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mxr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->mzr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->mzi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnz,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC3D_r4_t(const std::vector<float> &m_xr,
                                    const std::vector<float> &m_xi,
                                    const std::vector<float> &m_yr,
                                    const std::vector<float> &m_yi,
                                    const std::vector<float> &m_zr,
                                    const std::vector<float> &m_zi) {
                               
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               this->mnz = m_zr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               const std::size_t lenx = sizeof(float)*this->mnz;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                               std::memcpy(this->mzr,&m_zr[0],lenz);
                               std::memcpy(this->mzi,&m_zi[0],lenz);     
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC3D_r4_t(const std::valarray<float> &m_xr,
                                    const std::valarray<float> &m_xi,
                                    const std::valarray<float> &m_yr,
                                    const std::valarray<float> &m_yi,
                                    const std::valarray<float> &m_zr,
                                    const std::valarray<float> &m_zi) {
                               
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               this->mnz = m_zr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               const std::size_t lenz = sizeof(float)*this->mnz;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                               std::memcpy(this->mzr,&m_zr[0],lenz);
                               std::memcpy(this->mzi,&m_zi[0],lenz);     
                      }
                      
                      inline DC3D_r4_t(const std::size_t _mnx,
                                    const std::size_t _mny,
                                    const std::size_t _mnz,
                                    const float * __restrict m_xr,   
                                    const float * __restrict m_xi,
                                    const float * __restrict m_yr,
                                    const float * __restrict m_yi,
                                    const float * __restrict m_zr,
                                    const float * __restrict m_zi) {
                         
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  avx512_uncached_memmove(&this->myr[0],&m_yr[0],this->mny);
	                  avx512_uncached_memmove(&this->myi[0],&m_yi[0],this->mny);
	                  avx512_uncached_memmove(&this->mzr[0],&m_zr[0],this->mnz);
	                  avx512_uncached_memmove(&this->mzi[0],&m_zi[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  avx512_cached_memmove(&this->myr[0],&m_yr[0],this->mny);
	                  avx512_cached_memmove(&this->myi[0],&m_yi[0],this->mny);
	                  avx512_cached_memmove(&this->mzr[0],&m_zr[0],this->mnz);
	                  avx512_cached_memmove(&this->mzi[0],&m_zi[0],this->mnz);
#endif          
                      }
                      
                      inline DC3D_r4_t(DC3D_r4_t &&rhs) {
                           
                          this->mnx   = rhs.mnx;
                          this->mny   = rhs.mny;
                          this->mnz   = rhs.mnz;
                          this->mxr  = &rhs.mxr[0];
                          this->mxi  = &rhs.mxi[0];
                          this->myr  = &rhs.myr[0];
                          this->myi  = &rhs.myi[0];
                          this->mzr  = &rhs.mzr[0];
                          this->mzi  = &rhs.mzi[0];
                          rhs.mnx     = 0ULL;
                          rhs.mny     = 0ULL;
                          rhs.mnz     = 0ULL;
                          rhs.mxr    = NULL;
                          rhs.mxi    = NULL;
                          rhs.myr    = NULL;
                          rhs.myi    = NULL;
                          rhs.mzr    = NULL;
                          rhs.mzi    = NULL;
                      }   
                        
                      DC3D_r4_t(const DC3D_r4_t &)             = delete;
                      
                      inline ~DC3D_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap<float>(this->mxr,this->mnx);
                              gms_unmap<float>(this->mxi,this->mnx);
                              gms_unmap<float>(this->myr,this->mny); 
                              gms_unmap<float>(this->myi,this->mny);
                              gms_unmap<float>(this->mzr,this->mnz);
                              gms_unmap<float>(this->mzi,this->mnz);
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
                      
                      DC3D_r4_t & operator=(const DC3D_r4_t &) = delete;
                      
                      inline DC3D_r4_t & operator=(DC3D_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            if(this->ismmap) {
                                gms_unmap<float>(this->mxr,this->mnx);
                                gms_unmap<float>(this->mxi,this->mnx);
                                gms_unmap<float>(this->myr,this->mny); 
                                gms_unmap<float>(this->myi,this->mny);
                                gms_unmap<float>(this->mzr,this->mnz);
                                gms_unmap<float>(this->mzi,this->mnz);
                            }
                            else {
                               gms_mm_free(this->mxr);
                               gms_mm_free(this->mxi);
                               gms_mm_free(this->myr);
                               gms_mm_free(this->myi);
                               gms_mm_free(this->mzr);
                               gms_mm_free(this->mzi);
                            }
                            this->mnx   = rhs.mnx;
                            this->mny   = rhs.mny;
                            this->mnz   = rhs.mnz;
                            this->mxr   = &rhs.mxr[0];
                            this->mxi   = &rhs.mxi[0];
                            this->myr   = &rhs.myr[0];
                            this->myi   = &rhs.myi[0];
                            this->mzr   = &rhs.mzr[0];
                            this->mzi   = &rhs.mzi[0];
                            rhs.mnx     = 0ULL;
                            rhs.mny     = 0ULL;
                            rhs.mnz     = 0ULL;
                            rhs.mxr     = NULL;
                            rhs.mxi     = NULL;
                            rhs.myr     = NULL;
                            rhs.myi     = NULL;
                            rhs.mzr     = NULL;
                            rhs.mzi     = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->mxr  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->mxi  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->myr  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                        this->myi  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                        this->mzr  = (float*)gms_mm_malloc(sizeof(float)*this->mnz,64ULL);
                        this->mzi  = (float*)gms_mm_malloc(sizeof(float)*this->mnz,64ULL);
                   }    
                   
              };
              
              
             struct __ATTR_ALIGN__(64) DC2D_r4_t {
                     
                      float * __restrict mxr;
                      float * __restrict mxi;
                      float * __restrict myr;
                      float * __restrict myi;
                      std::size_t        mnx;
                      std::size_t        mny;
                      bool              ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif
                      inline DC2D_r4_t() {
                      
                         this->mnx   = 0ULL;
                         this->mny   = 0ULL;
                         this->mxr  = NULL;
                         this->mxi  = NULL;
                         this->myr  = NULL;
                         this->myi  = NULL;
                    }                    
                      
                      inline DC2D_r4_t(const std::size_t _nx,
                                       const std::size_t _ny) {
                                     
                             this->mnx = _nx;
                             this->mny = _ny;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline DC2D_r4_t(const std::size_t _nx,
                                    const std::size_t _ny,
                                    const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _nx;
                             this->mny = _ny;
                             switch (fsize) {
                                 case:0
                                      this->mxr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_4KiB<float>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mxr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_2MiB<float>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mxr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->myr = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->myi = (float*)
                                                  gms_mmap_1GiB<float>(this->mny,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC2D_r4_t(const std::vector<float> &m_xr,
                                    const std::vector<float> &m_xi,
                                    const std::vector<float> &m_yr,
                                    const std::vector<float> &m_yi) {
                                    
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                                
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC2D_r4_t(const std::valarray<float> &m_xr,
                                    const std::valarray<float> &m_xi,
                                    const std::valarray<float> &m_yr,
                                    const std::valarray<float> &m_yi) {
                                   
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                               std::memcpy(this->myr,&m_yr[0],leny);
                               std::memcpy(this->myi,&m_yi[0],leny);
                                 
                      }
                      
                      inline DC2D_r4_t(const std::size_t _mnx,
                                    const std::size_t _mny,
                                    const float * __restrict m_xr,   
                                    const float * __restrict m_xi,
                                    const float * __restrict m_yr,
                                    const float * __restrict m_yi) {                                   
                         
                          this->mnx = _mnx;
                          this->mny = _mny;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  avx512_uncached_memmove(&this->myr[0],&m_yr[0],this->mny);
	                  avx512_uncached_memmove(&this->myi[0],&m_yi[0],this->mny);
#else
	                  avx512_cached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  avx512_cached_memmove(&this->myr[0],&m_yr[0],this->mny);
	                  avx512_cached_memmove(&this->myi[0],&m_yi[0],this->mny);
	                
#endif          
                      }
                      
                      inline DC2D_r4_t(DC2D_r4_t &&rhs) {
                           
                          this->mnx   = rhs.mnx;
                          this->mny   = rhs.mny;
                          this->mxr  = &rhs.mxr[0];
                          this->mxi  = &rhs.mxi[0];
                          this->myr  = &rhs.myr[0];
                          this->myi  = &rhs.myi[0];
                          rhs.mnx     = 0ULL;
                          rhs.mny     = 0ULL;
                          rhs.mxr    = NULL;
                          rhs.mxi    = NULL;
                          rhs.myr    = NULL;
                          rhs.myi    = NULL;
                          
                      }   
                        
                      DC2D_r4_t(const DC2D_r4_t &)             = delete;
                      
                      inline ~DC2D_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap<float>(this->mxr,this->mnx);
                              gms_unmap<float>(this->mxi,this->mnx);
                              gms_unmap<float>(this->myr,this->mny); 
                              gms_unmap<float>(this->myi,this->mny);
                              
                           }
                           else {
                               gms_mm_free(this->mxr);
                               gms_mm_free(this->mxi);
                               gms_mm_free(this->myr);
                               gms_mm_free(this->myi);
                               
                           }
                      }
                      
                      DC2D_r4_t & operator=(const DC2D_r4_t &) = delete;
                      
                      inline DC2D_r4_t & operator=(DC2D_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            if(this->ismmap) {
                                gms_unmap<float>(this->mxr,this->mnx);
                                gms_unmap<float>(this->mxi,this->mnx);
                                gms_unmap<float>(this->myr,this->mny); 
                                gms_unmap<float>(this->myi,this->mny);
                                
                            }
                            else {
                               gms_mm_free(this->mxr);
                               gms_mm_free(this->mxi);
                               gms_mm_free(this->myr);
                               gms_mm_free(this->myi);
                              
                            }
                            this->mnx   = rhs.mnx;
                            this->mny   = rhs.mny;
                            this->mxr   = &rhs.mxr[0];
                            this->mxi   = &rhs.mxi[0];
                            this->myr   = &rhs.myr[0];
                            this->myi   = &rhs.myi[0];
                            rhs.mnx     = 0ULL;
                            rhs.mny     = 0ULL;
                            rhs.mxr     = NULL;
                            rhs.mxi     = NULL;
                            rhs.myr     = NULL;
                            rhs.myi     = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->mxr  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->mxi  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->myr  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                        this->myi  = (float*)gms_mm_malloc(sizeof(float)*this->mny,64ULL);
                       
                   }    
                   
              };
              
              
              
               struct __ATTR_ALIGN__(64) DC1D_r4_t {
                     
                      float * __restrict mxr;
                      float * __restrict mxi;
                      std::size_t        mnx;
                      bool              ismmap;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,39)
#endif
                      inline DC1D_r4_t() {
                      
                         this->mnx   = 0ULL;
                         this->mxr  = NULL;
                         this->mxi  = NULL;
                        
                    }                    
                      
                      inline DC1D_r4_t(const std::size_t _nx) {
                                       
                             this->mnx = _nx;
                             allocate();
                             this->ismmap = false;
                      }  
                      
                      inline DC1D_r4_t(const std::size_t _nx,
                                       const int32_t prot,
                                    const int32_t flags,
                                    const int32_t fd,
                                    const int32_t offset,
                                    const int32_t fsize) {
                             using namespace gms::common;
                             this->mnx = _nx;
                             switch (fsize) {
                                 case:0
                                      this->mxr = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:1
                                      this->mxr = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case:2
                                      this->mxr = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->mxi = (float*)
                                                  gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     } 
                       
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC1D_r4_t(const std::vector<float> &m_xr,
                                    const std::vector<float> &m_xi) {
                                    
                               this->mnx = m_xr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                                                              
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline DC1D_r4_t(const std::valarray<float> &m_xr,
                                    const std::valarray<float> &m_xi) {
                                  
                               this->mnx = m_xr.size(); 
                               allocate();
                               this->ismmap = false;
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               std::memcpy(this->mxr,&m_xr[0],lenx);
                               std::memcpy(this->mxi,&m_xi[0],lenx);
                                                               
                      }
                      
                      inline DC1D_r4_t(const std::size_t _mnx,
                                       const float * __restrict m_xr,   
                                    const float * __restrict m_xi) {
                                                             
                          this->mnx = _mnx;
                          allocate()
                          this->ismmap = false;
#if (USE_GMS_ANTENNA_MODEL_COMMON_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  
#else
	                  avx512_cached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  
	                
#endif          
                      }
                      
                      inline DC1D_r4_t(DC1D_r4_t &&rhs) {
                           
                          this->mnx   = rhs.mnx;
                          this->mxr  = &rhs.mxr[0];
                          this->mxi  = &rhs.mxi[0];
                          rhs.mnx     = 0ULL;
                          rhs.mxr    = NULL;
                          rhs.mxi    = NULL;
                                                    
                      }   
                        
                      DC1D_r4_t(const DC1D_r4_t &)             = delete;
                      
                      inline ~DC1D_r4_t() {
                           using namespace gms::common;
                           if(this->ismmap) {
                              gms_unmap<float>(this->mxr,this->mnx);
                              gms_unmap<float>(this->mxi,this->mnx);
                                                            
                           }
                           else {
                               gms_mm_free(this->mxr);
                               gms_mm_free(this->mxi);
                                                           
                           }
                      }
                      
                      DC1D_r4_t & operator=(const DC1D_r4_t &) = delete;
                      
                      inline DC1D_r4_t & operator=(DC1D_r4_t &&rhs) {
                            using namespace gms::common;
                            if(this==&rhs) return (*this);
                            if(this->ismmap) {
                                gms_unmap<float>(this->mxr,this->mnx);
                                gms_unmap<float>(this->mxi,this->mnx);
                                                               
                            }
                            else {
                               gms_mm_free(this->mxr);
                               gms_mm_free(this->mxi);
                                                             
                            }
                            this->mnx   = rhs.mnx;
                            this->mxr   = &rhs.mxr[0];
                            this->mxi   = &rhs.mxi[0];
                            rhs.mnx     = 0ULL;
                            rhs.mxr     = NULL;
                            rhs.mxi     = NULL;
                            return (*this);
                      }
                      
                   inline void allocate() {
                        using namespace gms::common;
                        this->mxr  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                        this->mxi  = (float*)gms_mm_malloc(sizeof(float)*this->mnx,64ULL);
                                             
                   }    
                   
              };
              
              

              
              




} //gms

















#endif /*__GMS_DYN_CONTAINERS_HPP__*/
