
#ifndef __GMS_DYN_ARRAY_H__
#define __GMS_DYN_ARRAY_H__ 180820251017



namespace file_info 
{

     static const unsigned int GMS_DYN_ARRAY_MAJOR = 1;
     static const unsigned int GMS_DYN_ARRAY_MINOR = 1;
     static const unsigned int GMS_DYN_ARRAY_MICRO = 0;
     static const unsigned int GMS_DYN_ARRAY_FULLVER =
       1000U*GMS_DYN_ARRAY_MAJOR+100U*GMS_DYN_ARRAY_MINOR+
       10U*GMS_DYN_ARRAY_MICRO;
     static const char GMS_DYN_ARRAY_CREATION_DATE[] = "18-08-2025 10:17 +00200 (MON  18 AUG 2025 GMT+2)";
     static const char GMS_DYN_ARRAY_BUILD_DATE[]    = __DATE__; 
     static const char GMS_DYN_ARRAY_BUILD_TIME[]  = __TIME__;
     static const char GMS_DYN_ARRAY_SYNOPSIS[]    = "Dynamically allocated, 64-byte aligned array.";

}





#include <cstdint>
#include <complex>
#include <vector>
#include <exception> //std::terminate
#include <valarray>
#include <cstring> // std::memcpy
#include <iostream>
#include "GMS_config.h"
#include "GMS_malloc.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (USE_GMS_DYN_ARRAY_NT_STORES)
#define USE_GMS_DYN_ARRAY_NT_STORES 0
#endif

namespace gms 
{

               // Simple dynamically allocated arrays without vector-like functionality (growing,resizing, ...etc)
               // Used mainly for the baseband/narrowband signal representation and henceforth processing.

               struct alignas(64) darray_c4_t final 
               {
                      
                      std::complex<float> * __restrict m_data;
                      std::size_t                      mnx;
                      bool                             ismmap;

                     inline darray_c4_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                      
                     } 
                          
                     inline explicit darray_c4_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  darray_c4_t(const std::size_t nx,
                                         const int32_t prot,
                                         const int32_t flags,
                                         const int32_t fd,
                                         const long offset,
                                         const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->mnx = nx;
                             switch (fsize) {
                                 case 0:
                                      this->m_data = (std::complex<float>*)
                                                 gms_mmap_4KiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_data = (std::complex<float>*)
                                                 gms_mmap_2MiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_data = (std::complex<float>*)
                                                 gms_mmap_1GiB<std::complex<float>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline darray_c4_t(const std::vector<std::complex<float>> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                
                     }
                     
                     inline darray_c4_t(const std::valarray<std::complex<float>> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                    
                     }
                             
                             
                      
                     inline  darray_c4_t(const std::size_t nx,
                                         const std::complex<float> * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<float>)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
                     }  
                   
                     inline  darray_c4_t(darray_c4_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                                
                      }
                                 
                     darray_c4_t(const darray_c4_t &)             = delete;
                      
                     inline ~darray_c4_t() noexcept(false)
                     {
                          using namespace gms::common;
                          if(this->ismmap) 
                          { 
                             int32_t err1{};
                             err1 = gms_ummap<std::complex<float>>(this->m_data,this->mnx); this->m_data = NULL;
                             if(__builtin_expect(err1==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                               
                             }
                          }
                          else 
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
                     }    
                         
                                                                 
                     darray_c4_t & operator=(const darray_c4_t &) = delete;
                      
                     inline darray_c4_t & operator=(darray_c4_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                                                                                                                
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<float>)*this->mnx);
                    }

                    inline std::complex<float> * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline std::complex<float> * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (std::complex<float>*)
                                         gms_mm_malloc( sizeof(std::complex<float>)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_c4_t) = " << alignof(darray_c4_t) << '\n'
                                    << "sizeof(struct  darray_c4_t) = " << sizeof(darray_c4_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n"
                                       "&this->ismmap               ="  << &this->ismmap         << "\n";
                     }
                      
               };


               struct alignas(64) darray_c8_t final 
               {
                      
                      std::complex<double> * __restrict m_data;
                      std::size_t                      mnx;
                      bool                             ismmap;

                     inline darray_c8_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                      
                     } 
                          
                     inline explicit darray_c8_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  darray_c8_t(const std::size_t nx,
                                         const int32_t prot,
                                         const int32_t flags,
                                         const int32_t fd,
                                         const long offset,
                                         const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->mnx = nx;
                             switch (fsize) {
                                 case 0:
                                      this->m_data = (std::complex<double>*)
                                                 gms_mmap_4KiB<std::complex<double>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_data = (std::complex<double>*)
                                                 gms_mmap_2MiB<std::complex<double>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_data = (std::complex<double>*)
                                                 gms_mmap_1GiB<std::complex<double>>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline darray_c8_t(const std::vector<std::complex<double>> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                
                     }
                     
                     inline darray_c8_t(const std::valarray<std::complex<double>> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                    
                     }
                             
                             
                      
                     inline  darray_c8_t(const std::size_t nx,
                                         const std::complex<double> * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(std::complex<double>)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
                     }  
                   
                     inline  darray_c8_t(darray_c8_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                                
                      }
                                 
                     darray_c8_t(const darray_c8_t &)             = delete;
                      
                     inline ~darray_c8_t() noexcept(false)
                     {
                          using namespace gms::common;
                          if(this->ismmap) 
                          { 
                             int32_t err1{};
                             err1 = gms_ummap<std::complex<double>>(this->m_data,this->mnx); this->m_data = NULL;
                             if(__builtin_expect(err1==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                               
                             }
                          }
                          else 
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
                     }    
                         
                                                                 
                     darray_c8_t & operator=(const darray_c8_t &) = delete;
                      
                     inline darray_c8_t & operator=(darray_c8_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                                                                                                          
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(std::complex<double>)*this->mnx);
                    }

                    inline std::complex<double> * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline std::complex<double> * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (std::complex<double>*)
                                         gms_mm_malloc( sizeof(std::complex<double>)*this->mnx,64ULL);
                       
                    }
                    public: 
                    inline void info_size_alignment() const 
                    {
                          std::cout << "alignof(struct darray_c8_t) = " << alignof(darray_c8_t) << '\n'
                                    << "sizeof(struct  darray_c8_t) = " << sizeof(darray_c8_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n"
                                       "&this->ismmap               ="  << &this->ismmap         << "\n";
                    }
                      
               };

               struct alignas(64) darray_r4_t final 
               {
                      
                      float * __restrict               m_data;
                      std::size_t                      mnx;
                      bool                             ismmap;

                     inline darray_r4_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                      
                     } 
                          
                     inline explicit darray_r4_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  darray_r4_t(const std::size_t nx,
                                         const int32_t prot,
                                         const int32_t flags,
                                         const int32_t fd,
                                         const long offset,
                                         const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->mnx = nx;
                             switch (fsize) {
                                 case 0:
                                      this->m_data = (float*)
                                                 gms_mmap_4KiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_data = (float*)
                                                 gms_mmap_2MiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_data = (float*)
                                                 gms_mmap_1GiB<float>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline darray_r4_t(const std::vector<float> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                
                     }
                     
                    inline darray_r4_t(const std::valarray<float> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                    
                     }
                             
                             
                      
                     inline  darray_r4_t(const std::size_t nx,
                                         const float * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(float)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
                     }  
                   
                     inline  darray_r4_t(darray_r4_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                                
                      }
                                 
                     darray_r4_t(const darray_r4_t &)             = delete;
                      
                     inline ~darray_r4_t() noexcept(false)
                     {
                          using namespace gms::common;
                          if(this->ismmap) 
                          { 
                             int32_t err1{};
                             err1 = gms_ummap<float>(this->m_data,this->mnx); this->m_data = NULL;
                             if(__builtin_expect(err1==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                               
                             }
                          }
                          else 
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
                     }    
                         
                                                                 
                     darray_r4_t & operator=(const darray_r4_t &) = delete;
                      
                     inline darray_r4_t & operator=(darray_r4_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                                                                           
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(float)*this->mnx);
                    }

                    inline float * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline float * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (float*)
                                         gms_mm_malloc( sizeof(float)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_r4_t) = " << alignof(darray_r4_t) << '\n'
                                    << "sizeof(struct  darray_r4_t) = " << sizeof(darray_r4_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n"
                                       "&this->ismmap               ="  << &this->ismmap         << "\n";
                     }
                      
               };


               struct alignas(64) darray_r8_t final 
               {
                      
                      double * __restrict               m_data;
                      std::size_t                      mnx;
                      bool                             ismmap;

                     inline darray_r8_t() noexcept(true)
                     {
                          
                          this->mnx     = 0ULL;
                          this->m_data  = NULL;
                                      
                     } 
                          
                     inline explicit darray_r8_t(const std::size_t nx) noexcept(false)
                     {
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                      }  
                      
                     inline  darray_r8_t(const std::size_t nx,
                                         const int32_t prot,
                                         const int32_t flags,
                                         const int32_t fd,
                                         const long offset,
                                         const int32_t fsize) noexcept(false)
                     {
                             using namespace gms::common;
                             this->mnx = nx;
                             switch (fsize) {
                                 case 0:
                                      this->m_data = (double*)
                                                 gms_mmap_4KiB<double>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 1:
                                      this->m_data = (double*)
                                                 gms_mmap_2MiB<double>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 case 2:
                                      this->m_data = (double*)
                                                 gms_mmap_1GiB<double>(this->mnx,prot,flags,fd,offset);
                                      this->ismmap = true;
                                 break;
                                 default :
                                      allocate();
                                      this->ismmap = false; // do not call mmap!!                        
                             }          
                     }
                      
                     inline darray_r8_t(const std::vector<double> &data) noexcept(false)
                     {    //shall be of the same size (no error checking implemented)
                                  
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                
                     }
                     
                    inline darray_r8_t(const std::valarray<double> &data) noexcept(false)
                    {
                                   
                          using namespace gms::common;
                          this->mnx = data.size();
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx{bytes_mnx()};
                          std::memcpy(this->m_data,&data[0],lenx);
                                                    
                     }
                             
                             
                      
                     inline  darray_r8_t(const std::size_t nx,
                                         const float * __restrict data) noexcept(false)
                     {
                                 
                          using namespace gms::common;
                          this->mnx = nx;
                          allocate();
                          this->ismmap = false;
                          const std::size_t lenx = sizeof(float)*this->mnx;
                          std::memcpy(this->m_data,&data[0],lenx);
                     }  
                   
                     inline  darray_r8_t(darray_r8_t && rhs) noexcept(true)
                     {
                          
                          this->mnx    =  rhs.mnx;
                          this->m_data =  &rhs.m_data[0];
                          rhs.mnx      =  0ULL;
                          rhs.m_data   =  NULL;
                                                
                      }
                                 
                     darray_r8_t(const darray_r8_t &)             = delete;
                      
                     inline ~darray_r8_t() noexcept(false)
                     {
                          using namespace gms::common;
                          if(this->ismmap) 
                          { 
                             int32_t err1{};
                             err1 = gms_ummap<double>(this->m_data,this->mnx); this->m_data = NULL;
                             if(__builtin_expect(err1==-1,0))
                             {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                               
                             }
                          }
                          else 
                          {
                             gms_mm_free(this->m_data); this->m_data = NULL;  
                          }  
                     }    
                         
                                                                 
                     darray_r8_t & operator=(const darray_r8_t &) = delete;
                      
                     inline darray_r8_t & operator=(darray_r8_t &&rhs) noexcept(true)
                     {
                           using namespace gms::common;
                           if(this==&rhs) return (*this);
                                                                           
                           gms_swap(this->mnx,rhs.mnx);
                           gms_swap(this->m_data,rhs.m_data);
                           rhs.mnx       = 0ULL;
                           rhs.m_data    = NULL;
                           return (*this);
                      }

                    inline std::size_t bytes_mnx() const noexcept(true)
                    {
                         return (sizeof(double)*this->mnx);
                    }

                    inline double * begin() const noexcept
                    {
                         return this->m_data;
                    }

                    inline double * end() const noexcept
                    {
                         return this->m_data+this->mnx;
                    }

                    inline constexpr std::size_t size() const noexcept
                    {
                           return this->mnx;
                    }

                    inline bool empty() const noexcept
                    {
                         return (this->size() == 0ull);
                    }

                    private:
                    inline void allocate() noexcept(false)
                    {
                          using namespace gms::common;
                          this->m_data  = (double*)
                                         gms_mm_malloc( sizeof(double)*this->mnx,64ULL);
                       
                     }
                     public: 
                     inline void info_size_alignment() const 
                     {
                          std::cout << "alignof(struct darray_r8_t) = " << alignof(darray_r8_t) << '\n'
                                    << "sizeof(struct  darray_r8_t) = " << sizeof(darray_r8_t)  << '\n'
                                    << std::hex << std::showbase        << '\n'
                                    << "&this->m_data               ="  << (void*)&this->m_data  << "\n"
                                       "&this->mnx                  ="  << &this->mnx            << "\n"
                                       "&this->ismmap               ="  << &this->ismmap         << "\n";
                     }
                      
               };







} //gms 















#endif /*__GMS_DYN_ARRAY_H__*/