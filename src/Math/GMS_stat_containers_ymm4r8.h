
#ifndef __GMS_STAT_CONTAINERS_YMM4R8_HPP__
#define __GMS_STAT_CONTAINERS_YMM4R8_HPP__ 131220230809


namespace file_info 
{

     const unsigned int GMS_STAT_CONTAINERS_YMM4R8_MAJOR = 2;
     const unsigned int GMS_STAT_CONTAINERS_YMM4R8_MINOR = 0;
     const unsigned int GMS_STAT_CONTAINERS_YMM4R8_MICRO = 0;
     const unsigned int GMS_STAT_CONTAINERS_YMM4R8_FULLVER =
       1000U*GMS_STAT_CONTAINERS_YMM4R8_MAJOR+100U*GMS_STAT_CONTAINERS_YMM4R8_MINOR+
       10U*GMS_STAT_CONTAINERS_YMM4R8_MICRO;
     const char * const GMS_STAT_CONTAINERS_YMM4R8_CREATION_DATE = "13-12-2023 08:09 +00200 (WED 13 DEC 2023 GMT+2)";
     const char * const GMS_STAT_CONTAINERS_YMM4R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_STAT_CONTAINERS_YMM4R8_SYNOPSIS      = "Statically allocated, 32-byte aligned YMM4R8-based dense containers.";

}


#include <cstdint>
#include <valarray>
#include <algorithm>
#include "GMS_config.h"
#include "GMS_avx_memcpy.h"
#include "GMS_avx_uncached_memcpy.h"
#include "GMS_complex_ymm4r8.h"


namespace gms
{
              
              
           template<std::size_t Nx,uint32_t stream_store>
           struct  alignas(32) SC1D_ymm4c8_t final 
           {
                     
                typedef gms::math::ymm4c8_t&       reference;
                typedef const gms::math::ymm4c8_t& const_reference;
                std::size_t mnx = Nx;
                __ATTR_ALIGN__(32)  gms::math::ymm4c8_t mx[Nx ? Nx : 1ull];
                   
                   
                inline void copy_fill(gms::math::ymm4c8_t * __restrict m_x,
                                      const std::size_t n)
                {
                        using namespace gms::common;
                        if(__builtin_expect(this->mnx!=n,0)) return;
                        if constexpr (stream_store) // case: non-temporal store is true 
                        {

                             avx_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x[0]),8ull*n);
	                  }
                        else 
                        {
                             avx_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x[0]),8ull*n);
                        }

                     
                  }  
                      
                                        
                      
                  inline void copy_fill(const std::vector<gms::math::ymm4c8_t> &m_x)    //shall be of the same size (no error checking implemented)
                  {                 
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {

                             avx_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x.data()),8ull*m_x.size());
	                    }
                          else 
                          {
                             avx_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x.data()),8ull*m_x.size());
                          }
                          
                  }
                     
                  inline void copy_fill(const std::valarray<gms::math::ymm4c8_t> &m_x)
                  {                 
                         
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {

                             avx_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x[0]),8ull*m_x.size());
	                    }
                          else 
                          {

                             avx_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x[0]),8ull*m_x.size());
                          }
                                 
                  }
                             
                   
                                           
                  gms::math::ymm4c8_t  *  begin() const 
                  {
                          
                          return (std::__addressof(this->mx[0]));
                  } 
                    
                  gms::math::ymm4c8_t  * end() const 
                  {
                           return (std::__addressof(this->mx[this->mnx]));
                  }

                  inline void swap(SC1D_ymm4c8_t &other)
                  {
                        std::swap_ranges(begin(),end(),other.begin());
                  }

                  inline constexpr std::size_t size() const 
                  {
                        return this->mnx;
                  }

                  inline constexpr bool empty() const 
                  {
                        return (this->size() == 0ull);
                  }
                     
                  inline reference front()
                  {
                        return (*this->begin());
                  }  

                  inline const_reference front() const 
                  {
                        return (*this->begin());
                  }

                  inline reference back() 
                  {
                        return (this->mnx ? *(this->begin()-1ull) : *this->begin());
                  }

                  inline const_reference back() const 
                  {
                        return (this->mnx ? *(this->begin()-1ull) : *this->begin());
                  }              
                      
                    
                 inline    gms::math::ymm4c8_t   mx_ptr() const 
                 {
                          
                          return (std::__addressof(this->mx[0]));
                 } 
                    
                
                     
                    
      };
           
              
              

              
              




} //gms

















#endif /*__GMS_STAT_CONTAINERS_YMM4R8_HPP__*/
