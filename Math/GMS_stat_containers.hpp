
#ifndef __GMS_STAT_CONTAINERS_HPP__
#define __GMS_STAT_CONTAINERS_HPP__ 051220230457


namespace file_info {

     const unsigned int GMS_STAT_CONTAINERS_MAJOR = 1;
     const unsigned int GMS_STAT_CONTAINERS_MINOR = 1;
     const unsigned int GMS_STAT_CONTAINERS_MICRO = 0;
     const unsigned int GMS_STAT_CONTAINERS_FULLVER =
       1000U*GMS_STAT_CONTAINERS_MAJOR+100U*GMS_STAT_CONTAINERS_MINOR+
       10U*GMS_STAT_CONTAINERS_MICRO;
     const char * const GMS_STAT_CONTAINERS_CREATION_DATE = "05-12-2023 04:57 +00200 (TUE 05 DEC 2023 GMT+2)";
     const char * const GMS_STAT_CONTAINERS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_STAT_CONTAINERS_SYNOPSIS      = "Dynamically allocated, 64-byte aligned dense containers."

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
#include "GMS_simd_memops.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_STAT_CONTAINERS_NT_STORES)
#define USE_GMS_STAT_CONTAINERS_NT_STORES 0
#endif


namespace gms {

           template<std::size_t Nx,
                    std::size_t Ny,
                    std::size_t Nz>
           struct  SC3D_c4_t {
                     
                    std::size_t mnx = Nx;
                    std::size_t mny = Ny;
                    std::size_t mnz = Nz;  
                    __ATTR_ALIGN__(64)  std::complex<float> mx[(mnx == 0) ? 4 : mnx];
                    __ATTR_ALIGN__(64)  std::complex<float> my[(mny == 0) ? 4 : mny];
                    __ATTR_ALIGN__(64)  std::complex<float> mz[(mnz == 0) ? 4 : mnz];
                     
             
                    SC3D_c4_t() = delete;
                          
                    inline SC3D_c4_t(std::complex<float> m_x[Nx],
                                     std::complex<float> m_y[Ny],
                                     std::complex<float> m_z[Nz]) {
                          using namespace gms::common;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->mnz);
#endif                         
                      }  
                      
                      
                     inline SC3D_c4_t(const std::size_t _mnx,
                                      const std::size_t _mny,
                                      const std::size_t _mnz,
                                      std::complex<float> * __restrict m_x,
                                      std::complex<float> * __restrict m_y,
                                      std::complex<float> * __restrict m_z) {
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
                          this->mnz = _mnz;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&m_z[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&m_z[0],this->mnz);
#endif                         
                      }  
                      
                   
                      
                     inline SC3D_c4_t(const std::vector<std::complex<float>> &m_x,    //shall be of the same size (no error checking implemented)
                                      const std::vector<std::complex<float>> &m_y,    //shall be of the same size (no error checking implemented)
                                      const std::vector<std::complex<float>> &m_z) {  //shall be of the same size (no error checking implemented)
                          
                          this->mnx = m_x.size();
                          this->mny = m_y.size();
                          this->mnz = m_z.size();
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          const std::size_t lemnz = sizeof(std::complex<float>)*this->mnz;
                          std::memcpy(&this->mx[0],&m_x[0],lemnx);
                          std::memcpy(&this->my[0],&m_y[0],lemny);
                          std::memcpy(&this->mz[0],&m_z[0],lemnz);       
                     }
                     
                     inline SC3D_c4_t(const std::valarray<std::complex<float>> &m_x,
                                      const std::valarray<std::complex<float>> &m_y,
                                      const std::valarray<std::complex<float>> &m_z) {
                         
                          
                          this->mnx = m_x.size();
                          this->mny = m_y.size();
                          this->mnz = m_z.size();
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          const std::size_t lemnz = sizeof(std::complex<float>)*this->mnz;
                          std::memcpy(&this->mx[0],&m_x[0],lemnx);
                          std::memcpy(&this->my[0],&m_y[0],lemny);
                          std::memcpy(&this->mz[0],&m_z[0],lemnz);           
                     }
                             
                   
                              
                                                
                     SC3D_c4_t(const SC3D_c4_t &rhs)  {
                     
                          using namespace gms::common;
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  avx512_uncached_memmove(&this->mz[0],&rhs.mz[0],this->mnz);
#else
	                  avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  avx512_cached_memmove(&this->mz[0],&rhs.mz[0],this->mnz);
#endif                            
                     }           
                      
                    
                      
                     SC3D_c4_t & operator=(const SC3D_c4_t &rhs) {
                     
                           using namespace gms::common;
                           if(this == &rhs) return (*this);
                           this->mnx = rhs.mnx;
                           this->mny = rhs.mny;
                           this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
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
                      
                    
                     std::complex<float> __restrict *  mx_ptr() const {
                          
                          return (std::__addressof(this->mx[0]));
                    } 
                    
                     std::complex<float> __restrict * my_ptr() const {
                          
                          return (std::__addressof(this->my[0]));
                    }
                    
                     std::complex<float> __restrict * mz_ptr() const {
                         
                          return (std::__addressof(this->mz[0]));
                    } 
                     
                    
              };
              
              
           template<std::size_t Nx,
                    std::size_t Ny>
           struct  SC2D_c4_t {
                     
                    std::size_t mnx = Nx;
                    std::size_t mny = Ny;
                    __ATTR_ALIGN__(64)  std::complex<float> mx[(mnx == 0) ? 4 : mnx];
                    __ATTR_ALIGN__(64)  std::complex<float> my[(mny == 0) ? 4 : mny];
                   
                     
             
                    SC2D_c4_t() = delete;
                          
                    inline SC2D_c4_t(std::complex<float> m_x[Nx],
                                     std::complex<float> m_y[Ny]) {
                                     
                          using namespace gms::common;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                  
#endif                         
                      }  
                      
                      
                     inline SC2D_c4_t(const std::size_t _mnx,
                                      const std::size_t _mny,
                                      std::complex<float> * __restrict m_x,
                                      std::complex<float> * __restrict m_y)
                                      
                          using namespace gms::common;
                          this->mnx = _mnx;
                          this->mny = _mny;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&m_y[0],this->mny);
	                  
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&m_y[0],this->mny);
	                 
#endif                         
                      }  
                      
                   
                      
                     inline SC2D_c4_t(const std::vector<std::complex<float>> &m_x,    //shall be of the same size (no error checking implemented)
                                      const std::vector<std::complex<float>> &m_y){    //shall be of the same size (no error checking implemented)
                           
                          this->mnx = m_x.size();
                          this->mny = m_y.size();          
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          std::memcpy(&this->mx[0],&m_x[0],lemnx);
                          std::memcpy(&this->my[0],&m_y[0],lemny);
                           
                     }
                     
                     inline SC2D_c4_t(const std::valarray<std::complex<float>> &m_x,
                                      const std::valarray<std::complex<float>> &m_y) {
                              
                          this->mnx = m_x.size();
                          this->mny = m_y.size();                              
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          const std::size_t lemny = sizeof(std::complex<float>)*this->mny;
                          std::memcpy(&this->mx[0],&m_x[0],lemnx);
                          std::memcpy(&this->my[0],&m_y[0],lemny);
                             
                     }
                             
                   
                              
                                                
                    inline SC2D_c4_t(const SC2D_c4_t &rhs)  {
                     
                          using namespace gms::common;
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_uncached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                
#else
	                  avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  avx512_cached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                 
#endif                            
                     }           
                      
                    
                      
                    inline SC2D_c4_t & operator=(const SC2D_c4_t &rhs) {
                     
                           using namespace gms::common;
                           if(this == &rhs) return (*this);
                           this->mnx = rhs.mnx;
                           this->mny = rhs.mny;
                          
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                   avx512_uncached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  
#else
	                   avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                   avx512_cached_memmove(&this->my[0],&rhs.my[0],this->mny);
	                  
#endif                  
                          return (*this);           
                    }
                      
                    
                  inline   std::complex<float> __restrict *  mx_ptr() const {
                          
                          return (std::__addressof(this->mx[0]));
                    } 
                    
                  inline   std::complex<float> __restrict * my_ptr() const {
                          
                          return (std::__addressof(this->my[0]));
                    }
                    
                     
                     
                    
              };
              
              
           template<std::size_t Nx>
           struct  SC1D_c4_t {
                     
                    std::size_t mnx = Nx;
                    __ATTR_ALIGN__(64)  std::complex<float> mx[(mnx == 0) ? 4 : mnx];
                                   
                            
                    SC1D_c4_t() = delete;
                          
                    inline SC2D_c4_t(std::complex<float> m_x[Nx]) {
                                    
                          using namespace gms::common;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                 
	                
#else
	                  avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                
	                  
#endif                         
                      }  
                      
                      
                     inline SC1D_c4_t(const std::size_t _mnx,
                                      std::complex<float> * __restrict m_x) {
                                    
                          using namespace gms::common;
                          this->mnx = _mnx;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                  	                  
#else	                 
                          avx512_cached_memmove(&this->mx[0],&m_x[0],this->mnx);
	                                  
#endif                         
                      }  
                      
                   
                      
                     inline SC1D_c4_t(const std::vector<std::complex<float>> &m_x) {    //shall be of the same size (no error checking implemented)
                              
                          this->mnx = m_x.size();       
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          std::memcpy(&this->mx[0],&m_x[0],lemnx);
                                            
                     }
                     
                     inline SC1D_c4_t(const std::valarray<std::complex<float>> &m_x) {
                           
                          this->mnx = m_x.size();           
                          const std::size_t lemnx = sizeof(std::complex<float>)*this->mnx;
                          std::memcpy(&this->mx[0],&m_x[0],lemnx);
                     }
                             
                   
                              
                                                
                   inline  SC1D_c4_t(const SC1D_c4_t &rhs)  {
                     
                          using namespace gms::common;
                          this->mnx = rhs.mnx;
                                                 
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                
	               
#else
	                  avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                 
	                 
#endif                            
                     }           
                      
                    
                      
                  inline   SC1D_c4_t & operator=(const SC1D_c4_t &rhs) {
                     
                           using namespace gms::common;
                           if(this == &rhs) return (*this);
                           this->mnx = rhs.mnx;
                                                    
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                  
	                  
#else
	                   avx512_cached_memmove(&this->mx[0],&rhs.mx[0],this->mnx);
	                
	                  
#endif                  
                          return (*this);           
                    }
                      
                    
                  inline   std::complex<float> __restrict *  mx_ptr() const {
                          
                          return (std::__addressof(this->mx[0]));
                    } 
                    
                  
                    
                     
                     
                    
              };
              
              
              
               
             
               template<std::size_t Nx,
                        std::size_t Ny,
                        std::size_t Nz>
               struct  SC3D_r4_t {
                     
                      std::size_t   mnx = Nx;
                      std::size_t   mny = Ny;
                      std::size_t   mnz = Nz;
                     __ATTR_ALIGN__(64) float  mxr[(mnx == 0) ? 16 : mnx];
                     __ATTR_ALIGN__(64) float  mxi[(mny == 0) ? 16 : mny];
                     __ATTR_ALIGN__(64) float  myr[(myr == 0) ? 16 : myr];
                     __ATTR_ALIGN__(64) float  myi[(myi == 0) ? 16 : myi];
                     __ATTR_ALIGN__(64) float  mzr[(mzr == 0) ? 16 : mzr];
                     __ATTR_ALIGN__(64) float  mzi[(mzi == 0) ? 16 : mzi];
                     
                   
                      inline SC3D_r4_t() = delete;
                      
                      
                      
                     
                      //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline SC3D_r4_t(const std::vector<float> &m_xr,
                                    const std::vector<float> &m_xi,
                                    const std::vector<float> &m_yr,
                                    const std::vector<float> &m_yi,
                                    const std::vector<float> &m_zr,
                                    const std::vector<float> &m_zi) {
                               
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               this->mnz = m_zr.size();
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               const std::size_t lenx = sizeof(float)*this->mnz;
                               std::memcpy(&this->mxr[0],&m_xr[0],lenx);
                               std::memcpy(&this->mxi[0],&m_xi[0],lenx);
                               std::memcpy(&this->myr[0],&m_yr[0],leny);
                               std::memcpy(&this->myi[0],&m_yi[0],leny);
                               std::memcpy(&this->mzr[0],&m_zr[0],lenz);
                               std::memcpy(&this->mzi[0],&m_zi[0],lenz);     
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline SC3D_r4_t(const std::valarray<float> &m_xr,
                                    const std::valarray<float> &m_xi,
                                    const std::valarray<float> &m_yr,
                                    const std::valarray<float> &m_yi,
                                    const std::valarray<float> &m_zr,
                                    const std::valarray<float> &m_zi) {
                               
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               this->mnz = m_zr.size();
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               const std::size_t lenz = sizeof(float)*this->mnz;
                               std::memcpy(&this->mxr[0],&m_xr[0],lenx);
                               std::memcpy(&this->mxi[0],&m_xi[0],lenx);
                               std::memcpy(&this->myr[0],&m_yr[0],leny);
                               std::memcpy(&this->myi[0],&m_yi[0],leny);
                               std::memcpy(&this->mzr[0],&m_zr[0],lenz);
                               std::memcpy(&this->mzi[0],&m_zi[0],lenz);     
                      }
                      
                      inline SC3D_r4_t(const std::size_t _mnx,
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
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
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
                      
                     
                    inline  SC3D_r4_t(const SC3D_r4_t &rhs) {
                          
                         using namespace gms::common;
                         this->mnx = rhs.mnx;
                         this->mny = rhs.mny;
                         this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
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
                      
                    
                      
                   inline   SC3D_r4_t & operator=(const SC3D_r4_t &rhs) {
                      
                          using namespace gms::common;
                          if(this == &rhs) return (*this);
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          this->mnz = rhs.mnz;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
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
                    
                 inline   float * __restrict get_mxr() const {
                         
                         return (std::__addressof(this->mxr[0]));
                   }
                   
                 inline   float * __restrict get_mxi() const {
                         
                         return (std::__addressof(this->mxi[0]));
                   }
                   
                 inline   float * __restrict get_myr() const {
                         
                         return (std::__addressof(this->myr[0]));
                   }
                   
                inline    float * __restrict get_myi() const {
                         
                         return (std::__addressof(this->myi[0]));
                   }
                   
                 inline   float * __restrict get_mzr() const {
                         
                         return (std::__addressof(this->mzr[0]));
                   }
                   
                 inline   float * __restrict get_mzi() const {
                         
                         return (std::__addressof(this->mzi[0]));
                   }
                      
                     
                   
              };
              
              
              template<std::size_t Nx,
                        std::size_t Ny>
              struct  SC2D_r4_t {
                     
                      std::size_t        mnx = Nx;
                      std::size_t        mny = Ny;
                     __ATTR_ALIGN__(64) float  mxr[(mnx == 0) ? 16 : mnx];
                     __ATTR_ALIGN__(64) float  mxi[(mny == 0) ? 16 : mny];
                     __ATTR_ALIGN__(64) float  myr[(myr == 0) ? 16 : myr];
                     __ATTR_ALIGN__(64) float  myi[(myi == 0) ? 16 : myi];
                                       
                     inline SC2D_r4_t() = delete;
                                                        
                         //The length of arguments must be of the same size (no error checking is implemented)!!
                     inline SC2D_r4_t(const std::vector<float> &m_xr,
                                      const std::vector<float> &m_xi,
                                      const std::vector<float> &m_yr,
                                      const std::vector<float> &m_yi) {
                                  
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               std::memcpy(&this->mxr[0],&m_xr[0],lenx);
                               std::memcpy(&this->mxi[0],&m_xi[0],lenx);
                               std::memcpy(&this->myr[0],&m_yr[0],leny);
                               std::memcpy(&this->myi[0],&m_yi[0],leny);
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline SC2D_r4_t(const std::valarray<float> &m_xr,
                                    const std::valarray<float> &m_xi,
                                    const std::valarray<float> &m_yr,
                                    const std::valarray<float> &m_yi) {
                                    
                               
                               this->mnx = m_xr.size(); 
                               this->mny = m_yr.size();
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               const std::size_t leny = sizeof(float)*this->mny;
                               std::memcpy(&this->mxr[0],&m_xr[0],lenx);
                               std::memcpy(&this->mxi[0],&m_xi[0],lenx);
                               std::memcpy(&this->myr[0],&m_yr[0],leny);
                               std::memcpy(&this->myi[0],&m_yi[0],leny);
                                  
                      }
                      
                      inline SC2D_r4_t(const std::size_t _mnx,
                                    const std::size_t _mny,
                                    const float * __restrict m_xr,   
                                    const float * __restrict m_xi,
                                    const float * __restrict m_yr,
                                    const float * __restrict m_yi) {
                                    
                          this->mnx = _mnx;
                          this->mny = _mny;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
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
                      
                     
                    inline  SC2D_r4_t(const SC2D_r4_t &rhs) {
                          
                         using namespace gms::common;
                         this->mnx = rhs.mnx;
                         this->mny = rhs.mny;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
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
                      
                    
                      
                   inline   SC2D_r4_t & operator=(const SC2D_r4_t &rhs) {
                      
                          using namespace gms::common;
                          if(this == &rhs) return (*this);
                          this->mnx = rhs.mnx;
                          this->mny = rhs.mny;
                          
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
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
                    
                 inline   float * __restrict get_mxr() const {
                         
                         return (std::__addressof(this->mxr[0]));
                   }
                   
                 inline   float * __restrict get_mxi() const {
                         
                         return (std::__addressof(this->mxi[0]));
                   }
                   
                 inline   float * __restrict get_myr() const {
                         
                         return (std::__addressof(this->myr[0]));
                   }
                   
                 inline   float * __restrict get_myi() const {
                         
                         return (std::__addressof(this->myi[0]));
                   }
                   
                   
            
         
             
              };
              
              
             
            template<std::size_t Nx>
            struct  SC1D_r4_t {
                     
                      std::size_t        mnx = Nx;
                     __ATTR_ALIGN__(64) float  mxr[(mnx == 0) ? 16 : mnx];
                     __ATTR_ALIGN__(64) float  mxi[(mny == 0) ? 16 : mny];
                    
                                       
                     inline SC1D_r4_t() = delete;
                                                        
                         //The length of arguments must be of the same size (no error checking is implemented)!!
                     inline SC1D_r4_t(const std::vector<float> &m_xr,
                                      const std::vector<float> &m_xi) {
                                     
                               this->mnx = m_xr.size(); 
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               std::memcpy(&this->mxr[0],&m_xr[0],lenx);
                               std::memcpy(&this->mxi[0],&m_xi[0],lenx);
                               
                      }
                      
                      
                       //The length of arguments must be of the same size (no error checking is implemented)!!
                      inline SC1D_r4_t(const std::valarray<float> &m_xr,
                                    const std::valarray<float> &m_xi) {
                                    
                               this->mnx = m_xr.size(); 
                               const std::size_t lenx = sizeof(float)*this->mnx;
                               std::memcpy(&this->mxr[0],&m_xr[0],lenx);
                               std::memcpy(&this->mxi[0],&m_xi[0],lenx);
                                                               
                      }
                      
                      inline SC1D_r4_t(const std::size_t _mnx,
                                       const float * __restrict m_xr,   
                                       const float * __restrict m_xi) {
                                      
                                   
                          this->mnx = _mnx;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                  
	                 
#else
	                  avx512_cached_memmove(&this->mxr[0],&m_xr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&m_xi[0],this->mnx);
	                 
	                  
#endif          
                      }
                      
                     
                    inline  SC1D_r4_t(const SC1D_r4_t &rhs) {
                          
                         using namespace gms::common;
                         this->mnx = rhs.mnx;
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                  avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                 
	                 
#else
	                  avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                  avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                  
	                  
#endif
                    }            
                      
                    
                      
                    inline  SC1D_r4_t & operator=(const SC1D_r4_t &rhs) {
                      
                          using namespace gms::common;
                          if(this == &rhs) return (*this);
                          this->mnx = rhs.mnx;
                          
                          
#if (USE_GMS_STAT_CONTAINERS_NT_STORES)  == 1
	                   avx512_uncached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_uncached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                   
	                   
#else
	                   avx512_cached_memmove(&this->mxr[0],&rhs.mxr[0],this->mnx);
	                   avx512_cached_memmove(&this->mxi[0],&rhs.mxi[0],this->mnx);
	                 
	                   
#endif                   
                           return (*this);       
                    }
                    
                 inline   float * __restrict get_mxr() const {
                         
                         return (std::__addressof(this->mxr[0]));
                   }
                   
                 inline  float * __restrict get_mxi() const {
                         
                         return (std::__addressof(this->mxi[0]));
                   }
                   
                    
             
              };
              
              
              
            
           
              
              

              
              




} //gms

















#endif /*__GMS_STAT_CONTAINERS_HPP__*/
