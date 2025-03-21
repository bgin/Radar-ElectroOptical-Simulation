
#ifndef __GMS_STAT_CONTAINERS_ZMM16R4_HPP__
#define __GMS_STAT_CONTAINERS_ZMM16R4_HPP__ 131220230651


namespace file_info {

     const unsigned int GMS_STAT_CONTAINERS_ZMM16R4_MAJOR = 1;
     const unsigned int GMS_STAT_CONTAINERS_ZMM16R4_MINOR = 1;
     const unsigned int GMS_STAT_CONTAINERS_ZMM16R4_MICRO = 0;
     const unsigned int GMS_STAT_CONTAINERS_ZMM16R4_FULLVER =
       1000U*GMS_STAT_CONTAINERS_ZMM16R4_MAJOR+100U*GMS_STAT_CONTAINERS_ZMM16R4_MINOR+
       10U*GMS_STAT_CONTAINERS_ZMM16R4_MICRO;
     const char * const GMS_STAT_CONTAINERS_ZMM16R4_CREATION_DATE = "13-12-2023 06:51 +00200 (WED 13 DEC 2023 GMT+2)";
     const char * const GMS_STAT_CONTAINERS_ZMM16R4_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_STAT_CONTAINERS_ZMM16R4_SYNOPSIS      = "Dynamically allocated, 64-byte aligned ZMM16R4-based dense containers."

}





#include <cstdint>
#include "GMS_config"
#include "GMS_complex_zmm8r8.hpp"
#include "GMS_simd_memops.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)
#define USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES 0
#endif


namespace gms {

           template<std::size_t Nx,
                    std::size_t Ny,
                    std::size_t Nz>
           struct  SC3D_zmm16c4_t {
                     
                    std::size_t mnx = Nx;
                    std::size_t mny = Ny;
                    std::size_t mnz = Nz;  
                    __ATTR_ALIGN__(64)  zmm16c4_t mx[(mnx == 0) ? 1 : mnx];
                    __ATTR_ALIGN__(64)  zmm16c4_t my[(mny == 0) ? 1 : mny];
                    __ATTR_ALIGN__(64)  zmm16c4_t mz[(mnz == 0) ? 1 : mnz];
                     
             
                    SC3D_zmm16c4_t() = delete;
                          
                    inline SC3D_zmm16c4_t(zmm16c4_t m_x[Nx],
                                         zmm16c4_t m_y[Ny],
                                         zmm16c4_t m_z[Nz]) {
                          using namespace gms::common;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->mnz);
#endif                         
                      }  
                      
                      
                     inline SC3D_zmm16c4_t(const std::size_t _mnx,
                                      const std::size_t _mny,
                                      const std::size_t _mnz,
                                      zmm16c4_t * __restrict m_x,
                                      zmm16c4_t * __restrict m_y,
                                      zmm16c4_t * __restrict m_z) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->mnz);
#endif                         
                      }  
                      
                   
                      
                    
                   
                              
                                                
                 inline    SC3D_zmm16c4_t(const SC3D_zmm16c4_t &rhs)  {
                     
                          using namespace gms::common;
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&rhs.mz[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&rhs.mz[0],this->mnz);
#endif                            
                     }           
                      
                    
                      
                 inline    SC3D_zmm16c4_t & operator=(const SC3D_zmm16c4_t &rhs) {
                     
                           using namespace gms::common;
                           if(this == &rhs) return (*this);
                           this->mnx = rhs.mnx;
                           this->mny = rhs.mny;
                           this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                   avx512_uncached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                   avx512_uncached_memmove(&this->mz[0],&rhs.mz[0],this->mnz);
#else
	                   avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                   avx512_cached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                   avx512_cached_memmove(&this->mz[0],&rhs.mz[0],this->mnz);
#endif                  
                          return (*this);           
                    }
                      
                    
                 inline    zmm16c4_t __restrict *  mx_ptr() const {
                          
                          return (std::__addressof(this->mx[0]));
                    } 
                    
                 inline    zmm16c4_t __restrict * my_ptr() const {
                          
                          return (std::__addressof(this->my[0]));
                    }
                    
                 inline    zmm16c4_t __restrict * mz_ptr() const {
                         
                          return (std::__addressof(this->mz[0]));
                    } 
                     
                    
              };
              
              
           template<std::size_t Nx,
                    std::size_t Ny>
           struct  SC2D_zmm16c4_t {
                     
                    std::size_t mnx = Nx;
                    std::size_t mny = Ny;
                    __ATTR_ALIGN__(64)  zmm16c4_t mx[(mnx == 0) ? 1 : mnx];
                    __ATTR_ALIGN__(64)  zmm16c4_t my[(mny == 0) ? 1 : mny];
                   
                     
             
                    SC2D_zmm16c4_t() = delete;
                          
                    inline SC2D_zmm16c4_t(zmm16c4_t m_x[Nx],
                                     zmm16c4_t m_y[Ny]) {
                                     
                          using namespace gms::common;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  
#endif                         
                      }  
                      
                      
                     inline SC2D_zmm16c4_t(const std::size_t _mnx,
                                      const std::size_t _mny,
                                      zmm16c4_t * __restrict m_x,
                                      zmm16c4_t * __restrict m_y)
                                 
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                 
#endif                         
                      }  
                      
                   
                                          
                inline  SC2D_zmm16c4_t(const SC2D_zmm16c4_t &rhs)  {
                     
                          using namespace gms::common;
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                
#else
	                  avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                 
#endif                            
                     }           
                      
                    
                      
                  inline  SC2D_zmm16c4_t & operator=(const SC2D_zmm16c4_t &rhs) {
                     
                           using namespace gms::common;
                           if(this == &rhs) return (*this);
                           this->mnx = rhs.mnx;
                           this->mny = rhs.mny;
                          
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                   avx512_uncached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  
#else
	                   avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                   avx512_cached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  
#endif                  
                          return (*this);           
                    }
                      
                    
                 inline    zmm16c4_t __restrict *  mx_ptr() const {
                          
                          return (std::__addressof(this->mx[0]));
                    } 
                    
                 inline    zmm16c4_t __restrict * my_ptr() const {
                          
                          return (std::__addressof(this->my[0]));
                    }
                    
                     
                     
                    
              };
              
              
           template<std::size_t Nx>
           struct  SC1D_zmm16c4_t {
                     
                    std::size_t mnx = Nx;
                    __ATTR_ALIGN__(64)  zmm16c4_t mx[(mnx == 0) ? 1 : mnx];
                                   
                            
                    SC1D_zmm16c4_t() = delete;
                          
                    inline SC1D_zmm16c4_t(zmm16c4_t m_x[Nx]) {
                                    
                          using namespace gms::common;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	                
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                
	                  
#endif                         
                      }  
                      
                      
                     inline SC1D_zmm16c4_t(const std::size_t _mnx,
                                      zmm16c4_t * __restrict m_x) {
                                    
                          using namespace gms::common;
                          this->mnx = _mnx;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  	                  
#else	                 
                          avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                                  
#endif                         
                      }  
                      
                   
                      
                                                 
                                                
                   inline   SC1D_zmm16c4_t(const SC1D_zmm16c4_t &rhs)  {
                     
                          using namespace gms::common;
                          this->mnx = rhs.mnx;
                                                 
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                
	               
#else
	                  avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                 
	                 
#endif                            
                     }           
                      
                    
                      
                  inline   SC1D_zmm16c4_t & operator=(const SC1D_zmm16c4_t &rhs) {
                     
                           using namespace gms::common;
                           if(this == &rhs) return (*this);
                           this->mnx = rhs.mnx;
                                                    
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  
	                  
#else
	                   avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                
	                  
#endif                  
                          return (*this);           
                    }
                      
                    
                  inline   zmm16c4_t __restrict *  mx_ptr() const {
                          
                          return (std::__addressof(this->mx[0]));
                    } 
                    
                  
                    
                     
                     
                    
              };
              
              
              
               
             
               template<std::size_t Nx,
                        std::size_t Ny,
                        std::size_t Nz>
               struct  SC3D__m512_t {
                     
                      std::size_t mnx = Nx;
                      std::size_t mny = Ny;
                      std::size_t mnz = Nz;
                     __ATTR_ALIGN__(64) __m512  mxr[(mnx == 0) ? 1 : mnx];
                     __ATTR_ALIGN__(64) __m512  mxi[(mny == 0) ? 1 : mny];
                     __ATTR_ALIGN__(64) __m512  myr[(myr == 0) ? 1 : myr];
                     __ATTR_ALIGN__(64) __m512  myi[(myi == 0) ? 1 : myi];
                     __ATTR_ALIGN__(64) __m512  mzr[(mzr == 0) ? 1 : mzr];
                     __ATTR_ALIGN__(64) __m512  mzi[(mzi == 0) ? 1 : mzi];
                     
                   
                      inline SC3D__m512_t() = delete;
                      
                      
                                        
                      
                      inline SC3D__m512_t(const std::size_t _mnx,
                                    const std::size_t _mny,
                                    const std::size_t _mnz,
                                    const __m512 * __restrict m_xr,   
                                    const __m512 * __restrict m_xi,
                                    const __m512 * __restrict m_yr,
                                    const __m512 * __restrict m_yi,
                                    const __m512 * __restrict m_zr,
                                    const __m512 * __restrict m_zi) {
                         
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
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
                      
                     
                     inline SC3D__m512_t(const SC3D__m512_t &rhs) {
                          
                         using namespace gms::common;
                         this->mnx = rhs.mnx;
                         this->mny = rhs.mny;
                         this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                  avx512_uncached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                  avx512_uncached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                  avx512_uncached_memmove(&this->mzr[0],&rhs.mzr[0],this->mnz);
	                  avx512_uncached_memmove(&this->mzi[0],&rhs.mzi[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                  avx512_cached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                  avx512_cached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                  avx512_cached_memmove(&this->mzr[0],&rhs.mzr[0],this->mnz);
	                  avx512_cached_memmove(&this->mzi[0],&rhs.mzi[0],this->mnz);
#endif
                    }            
                      
                    
                      
                   inline   SC3D__m512_t & operator=(const SC3D__m512_t &rhs) {
                      
                          using namespace gms::common;
                          if(this == &rhs) return (*this);
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                   avx512_uncached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                   avx512_uncached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                   avx512_uncached_memmove(&this->mzr[0],&rhs.mzr[0],this->mnz);
	                   avx512_uncached_memmove(&this->mzi[0],&rhs.mzi[0],this->mnz);
#else
	                   avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                   avx512_cached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                   avx512_cached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                   avx512_cached_memmove(&this->mzr[0],&rhs.mzr[0],this->mnz);
	                   avx512_cached_memmove(&this->mzi[0],&rhs.mzi[0],this->mnz);
#endif                   
                           return (*this);       
                    }
                    
                  inline  __m512 * __restrict get_mxr() const {
                         
                         return (std::__addressof(this->mxr[0]));
                   }
                   
                  inline  __m512 * __restrict get_mxi() const {
                         
                         return (std::__addressof(this->mxi[0]));
                   }
                   
                 inline   __m512 * __restrict get_myr() const {
                         
                         return (std::__addressof(this->myr[0]));
                   }
                   
                inline    __m512 * __restrict get_myi() const {
                         
                         return (std::__addressof(this->myi[0]));
                   }
                   
                inline    __m512 * __restrict get_mzr() const {
                         
                         return (std::__addressof(this->mzr[0]));
                   }
                   
                inline    __m512 * __restrict get_mzi() const {
                         
                         return (std::__addressof(this->mzi[0]));
                   }
                      
                     
                   
              };
              
              
              template<std::size_t Nx,
                       std::size_t Ny>
              struct  SC2D__m512_t {
                     
                      std::size_t  mnx = Nx;
                      std::size_t  mny = Ny;
                     __ATTR_ALIGN__(64) __m512  mxr[(mnx == 0) ? 1 : mnx];
                     __ATTR_ALIGN__(64) __m512  mxi[(mny == 0) ? 1 : mny];
                     __ATTR_ALIGN__(64) __m512  myr[(myr == 0) ? 1 : myr];
                     __ATTR_ALIGN__(64) __m512  myi[(myi == 0) ? 1 : myi];
                                       
                     inline SC2D__m512_t() = delete;
                                                        
                    
                                           
                     inline SC2D__m512_t(const std::size_t _mnx,
                                    const std::size_t _mny,
                                    const __m512 * __restrict m_xr,   
                                    const __m512 * __restrict m_xi,
                                    const __m512 * __restrict m_yr,
                                    const __m512 * __restrict m_yi) {
                                    
                          this->mnx = _mnx;
                          this->mny = _mny;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
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
                      
                     
                    inline  SC2D__m512_t(const SC2D__m512_t &rhs) {
                          
                         using namespace gms::common;
                         this->mnx = rhs.mnx;
                         this->mny = rhs.mny;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                  avx512_uncached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                  avx512_uncached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                 
#else
	                  avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                  avx512_cached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                  avx512_cached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                  
#endif
                    }            
                      
                    
                      
                    inline  SC2D__m512_t & operator=(const SC2D__m512_t &rhs) {
                      
                          using namespace gms::common;
                          if(this == &rhs) return (*this);
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                   avx512_uncached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                   avx512_uncached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                   
#else
	                   avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                   avx512_cached_memmove(&this->myr[0],&rhs.myr[0],this->mny);
	                   avx512_cached_memmove(&this->myi[0],&rhs.myi[0],this->mny);
	                   
#endif                   
                           return (*this);       
                    }
                    
                 inline   __m512 * __restrict get_mxr() const {
                         
                         return (std::__addressof(this->mxr[0]));
                   }
                   
                 inline   __m512 * __restrict get_mxi() const {
                         
                         return (std::__addressof(this->mxi[0]));
                   }
                   
                 inline   __m512 * __restrict get_myr() const {
                         
                         return (std::__addressof(this->myr[0]));
                   }
                   
                 inline   __m512 * __restrict get_myi() const {
                         
                         return (std::__addressof(this->myi[0]));
                   }
                   
                   
            
         
             
              };
              
              
             
            template<std::size_t Nx>
            struct  SC1D__m512_t {
                     
                      std::size_t  mnx = Nx;
                     __ATTR_ALIGN__(64) __m512  mxr[(mnx == 0) ? 1 : mnx];
                     __ATTR_ALIGN__(64) __m512  mxi[(mny == 0) ? 1 : mny];
                    
                                       
                     inline SC1D__m512_t() = delete;
                                                        
                       
                    
                      
                      inline SC1D__m512_t(const std::size_t _mnx,
                                       const __m512 * __restrict m_xr,   
                                       const __m512 * __restrict m_xi) {
                                      
                          this->mnx = _mnx;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  
	                 
#else
	                  avx512_cached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                 
	                  
#endif          
                      }
                      
                     
                    inline  SC1D__m512_t(const SC1D__m512_t &rhs) {
                          
                         using namespace gms::common;
                         this->mnx = rhs.mnx;
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                 
	                 
#else
	                  avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                  
	                  
#endif
                    }            
                      
                    
                      
                   inline   SC1D__m512_t & operator=(const SC1D__m512_t &rhs) {
                      
                          using namespace gms::common;
                          if(this == &rhs) return (*this);
                          this->mnx = rhs.mnx;
                        
#if (USE_GMS_STAT_CONTAINERS_ZMM16R4_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                   
	                   
#else
	                   avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                 
	                   
#endif                   
                           return (*this);       
                    }
                    
                 inline   __m512 * __restrict get_mxr() const {
                         
                         return (std::__addressof(this->mxr[0]));
                   }
                   
                 inline   __m512 * __restrict get_mxi() const {
                         
                         return (std::__addressof(this->mxi[0]));
                   }
                   
                    
             
              };
              
              
              
            
           
              
              

              
              




} //gms

















#endif /*__GMS_STAT_CONTAINERS_ZMM16R4_HPP__*/
