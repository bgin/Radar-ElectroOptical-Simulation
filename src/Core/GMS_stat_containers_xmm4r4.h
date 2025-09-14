
#ifndef __GMS_STAT_CONTAINERS_XMM4R4_H__
#define __GMS_STAT_CONTAINERS_XMM4R4_H__ 131220230139


namespace file_info 
{

     const unsigned int GMS_STAT_CONTAINERS_XMM4R4_MAJOR = 1;
     const unsigned int GMS_STAT_CONTAINERS_XMM4R4_MINOR = 1;
     const unsigned int GMS_STAT_CONTAINERS_XMM4R4_MICRO = 0;
     const unsigned int GMS_STAT_CONTAINERS_XMM4R4_FULLVER =
       1000U*GMS_STAT_CONTAINERS_XMM4R4_MAJOR+100U*GMS_STAT_CONTAINERS_XMM4R4_MINOR+
       10U*GMS_STAT_CONTAINERS_XMM4R4_MICRO;
     const char * const GMS_STAT_CONTAINERS_XMM4R4_CREATION_DATE = "13-12-2023 01:39 +00200 (WED 13 DEC 2023 GMT+2)";
     const char * const GMS_STAT_CONTAINERS_XMM4R4_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_STAT_CONTAINERS_XMM4R4_SYNOPSIS      = "Statically allocated, 16-byte aligned XMM4R4-based dense containers.";

}





#include <cstdint>
#include <valarray>
#include <algorithm>
#include "GMS_config.h"
#include "GMS_sse_memcpy.h"
#include "GMS_sse_uncached_memcpy.h"
#include "GMS_complex_xmm4r4.h"





namespace gms 
{

          
              
           template<std::size_t Nx,uint32_t stream_store>
           struct  alignas(16) SC1D_xmm4c4_t final 
           {
                     
                typedef gms::math::xmm4c4_t&       reference;
                typedef const gms::math::xmm4c4_t& const_reference;
                std::size_t mnx = Nx;
                __ATTR_ALIGN__(16)  gms::math::xmm4c4_t mx[Nx ? Nx : 1ull];
                   
                   
                inline void copy_fill(gms::math::xmm4c4_t * __restrict m_x,
                                      const std::size_t n)
                {
                        using namespace gms::common;
                        if(__builtin_expect(this->mnx!=n,0)) return;
                        if constexpr (stream_store) // case: non-temporal store is true 
                        {

                             sse_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),8ull*this->mnx);
	                  }
                        else 
                        {
                             sse_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                    reinterpret_cast<float  *>(&m_x[0]),8ull*this->mnx);
                        }

                     
                  }  
                      
                                        
                      
                  inline void copy_fill(const std::vector<gms::math::xmm4c4_t> &m_x)    //shall be of the same size (no error checking implemented)
                  {                 
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {

                             sse_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x.data()),8ull*this->mnx);
	                    }
                          else 
                          {
                             sse_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                    reinterpret_cast<float  *>(&m_x.data()),8ull*this->mnx);
                          }
                          
                  }
                     
                  inline void copy_fill(const std::valarray<gms::math::xmm4c4_t> &m_x)
                  {                 
                         
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {

                             sse_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),8ull*this->mnx);
	                    }
                          else 
                          {

                             sse_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                    reinterpret_cast<float  *>(&m_x[0]),8ull*this->mnx);
                          }
                                 
                  }
                             
                   
                                           
                  gms::math::xmm4c4_t  *  begin() const 
                  {
                          
                          return (std::__addressof(this->mx[0]));
                  } 
                    
                  gms::math::xmm4c4_t  * end() const 
                  {
                           return (std::__addressof(this->mx[this->mnx]));
                  }

                  inline void swap(SC1D_xmm4c4_t &other)
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
                      
                    
                 inline    gms::math::xmm4c4_t   mx_ptr() const 
                 {
                          
                          return (std::__addressof(this->mx[0]));
                 } 
                    
                
                     
                    
              };
              
            
           
              
              

              
              




} //gms

















#endif /*__GMS_STAT_CONTAINERS_XMM4R4_H__*/
