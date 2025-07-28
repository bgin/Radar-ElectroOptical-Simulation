
#ifndef __GMS_STAT_CONTAINERS_H__
#define __GMS_STAT_CONTAINERS_H__ 051220230457


namespace file_info
{

     const unsigned int GMS_STAT_CONTAINERS_MAJOR = 2;
     const unsigned int GMS_STAT_CONTAINERS_MINOR = 0;
     const unsigned int GMS_STAT_CONTAINERS_MICRO = 0;
     const unsigned int GMS_STAT_CONTAINERS_FULLVER =
       1000U*GMS_STAT_CONTAINERS_MAJOR+100U*GMS_STAT_CONTAINERS_MINOR+
       10U*GMS_STAT_CONTAINERS_MICRO;
     const char * const GMS_STAT_CONTAINERS_CREATION_DATE = "05-12-2023 04:57 +00200 (TUE 05 DEC 2023 GMT+2)";
     const char * const GMS_STAT_CONTAINERS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_STAT_CONTAINERS_SYNOPSIS      = "Dynamically allocated, 64-byte aligned dense containers.";

}





#include <cstdint>
#include <complex>
#include <vector>
#include <valarray>
#include <algorithm>
#include "GMS_config.h"
#if defined (__AVX512F__)
#include "GMS_avx512_memcpy.h"
#include "GMS_avx512_uncached_memcpy.h"
#elif defined (__AVX__) || defined (__AVX2__)
#include "GMS_avx_memcpy.h"
#include "GMS_avx_uncached_memcpy.h"
#elif defined (SSE)
#include "GMS_sse_memcpy.h"
#include "GMS_sse_uncached_memcpy.h"
#endif 




namespace gms 
{

           template<std::size_t Nx,uint32_t stream_store>
           struct  alignas(64) StatC1D_c4_t 
           {
                     
                   typedef std::complex<float>&       reference;
                   typedef const std::complex<float>& const_reference;

                   std::size_t mnx = Nx;
#if defined (__AVX512F__)
                     
                    __ATTR_ALIGN__(64)  std::complex<float> mx[Nx ? Nx : 8ull];
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32)  std::complex<float> mx[Nx ? Nx : 4ull];
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16)  std::complex<float> mx[Nx ? Nx : 2ull];
#endif 
                  // mstream_store = stream_store;
             
                    
                          
                  inline void copy_fill(std::complex<float> * __restrict m_x,
                                         const std::size_t n)
                  {
                        using namespace gms::common;
                        if(__builtin_expect(this->mnx!=n,0)) return;
                        if constexpr (stream_store) // case: non-temporal store is true 
                        {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                                reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                                reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#endif                             
                        }

                     
                      }  
                      
                                        
                      
                  inline void copy_fill(const std::vector<std::complex<float>> &m_x)    //shall be of the same size (no error checking implemented)
                  {                 
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                                reinterpret_cast<float  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x.data()),2ull*this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                       reinterpret_cast<float  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                    reinterpret_cast<float  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                    reinterpret_cast<float  *>(&m_x.data()),2ull*this->mnx);
#endif                             
                        }
                          
                  }
                     
                  inline void copy_fill(const std::valarray<std::complex<float>> &m_x)
                  {                 
                         
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                                reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                             reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                       reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                    reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_ps(reinterpret_cast<float  *>(&this->mx[0]),
                                                    reinterpret_cast<float  *>(&m_x[0]),2ull*this->mnx);
#endif                             
                        }
                                 
                  }
                             
                   
                                           
                  std::complex<float>  *  begin() const 
                  {
                          
                          return (std::__addressof(this->mx[0]));
                  } 
                    
                  std::complex<float>  * end() const 
                  {
                           return (std::__addressof(this->mx[this->mnx]));
                  }

                  inline void swap(StatC1D_c4_t &other)
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
              };
              
              
           template<std::size_t Nx,uint32_t stream_store>
           struct  alignas(64)  StatC1D_c8_t 
           {
                     
                   typedef std::complex<double>&       reference;
                   typedef const std::complex<double>& const_reference;

                   std::size_t mnx = Nx;
#if defined (__AVX512F__)
                      __ATTR_ALIGN__(64)  std::complex<double> mx[Nx ? Nx : 4ull];
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32)  std::complex<double> mx[Nx ? Nx : 2ull];
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16)  std::complex<double> mx[Nx ? Nx : 1ull];
#endif 
                  //bool mstream_store = stream_store;
             
                    
                          
                  inline void copy_fill(std::complex<double> * __restrict m_x,
                                         const std::size_t n)
                  {
                        using namespace gms::common;
                        if(__builtin_expect(this->mnx!=n,0)) return;
                        if constexpr (stream_store) // case: non-temporal store is true 
                        {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                                reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                                reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#endif                             
                        }

                     
                      }  
                      
                  
                      
                  inline void copy_fill(const std::vector<std::complex<double>> &m_x)    //shall be of the same size (no error checking implemented)
                  {                 
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                                reinterpret_cast<double  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x.data()),2ull*this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                       reinterpret_cast<double  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x.data()),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x.data()),2ull*this->mnx);
#endif                             
                        }
                          
                  }
                     

                  inline void copy_fill(const std::valarray<std::complex<double>> &m_x)
                  {                 
                         
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                                reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_pd(reinterpret_cast<double  *>(&this->mx[0]),
                                                             reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_ps(reinterpret_cast<double  *>(&this->mx[0]),
                                                       reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_ps(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_ps(reinterpret_cast<double  *>(&this->mx[0]),
                                                    reinterpret_cast<double  *>(&m_x[0]),2ull*this->mnx);
#endif                             
                        }
                                 
                  }
                             
                   
                                           
                  std::complex<double>  *  begin() const 
                  {
                          
                          return (std::__addressof(this->mx[0]));
                  } 
                    
                  std::complex<double>  * end() const 
                  {
                           return (std::__addressof(this->mx[this->mnx]));
                  }

                  inline void swap(StatC1D_c8_t &other)
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
              };
              
               
             
            template<std::size_t Nx,uint32_t stream_store>
            struct  alignas(64) StatV1D_r4_t 
            {
                     
                   typedef float&       reference;
                   typedef const float& const_reference;

                   std::size_t mnx = Nx;
#if defined (__AVX512F__)
                     
                    __ATTR_ALIGN__(64)  float mx[Nx ? Nx : 16ull];
#elif defined (__AVX__) || defined (__AVX2__)
                    __ATTR_ALIGN__(32)  float mx[Nx ? Nx : 8ull];
#elif defined (__SSE__)
                    __ATTR_ALIGN__(16)  float mx[Nx ? Nx : 4ull];
#endif 
                  // mstream_store = stream_store;
             
                    
                          
                  inline void copy_fill(float * __restrict m_x,
                                        const std::size_t n)
                  {
                        using namespace gms::common;
                        if(__builtin_expect(this->mnx!=n,0)) return;
                        if constexpr (stream_store) // case: non-temporal store is true 
                        {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                                &m_x[0],this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                             &m_x[0],this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                             &m_x[0],this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_ps(&this->mx[0]),
                                                       &m_x[0],this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_ps(&this->mx[0],
                                                    &m_x[0],this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_ps(&this->mx[0],
                                                    &m_x[0],this->mnx);
#endif                             
                        }

                     
                      }  
                      
                                        
                      
                  inline void copy_fill(const std::vector<float> &m_x)    //shall be of the same size (no error checking implemented)
                  {                 
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                                &m_x.data(),this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                             &m_x.data(),this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                             &m_x.data(),this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_ps(&this->mx[0],
                                                       &m_x.data(),this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_ps(&this->mx[0],
                                                    &m_x.data(),this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_ps(&this->mx[0],
                                                    &m_x.data(),this->mnx);
#endif                             
                        }
                          
                  }
                     
                  inline void copy_fill(const std::valarray<float> &m_x)
                  {                 
                         
                          using namespace gms::common;
                          if(__builtin_expect(this->mnx!=m_x.size(),0)) return;
                          if constexpr (stream_store) // case: non-temporal store is true 
                          {
#if defined (__AVX512F__)
	                       avx512_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                                &m_x[0],this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                             &m_x[0],this->mnx);
#elif defined (__SSE__)
                             sse_uncached_memcpy_unroll8x_ps(&this->mx[0],
                                                             &m_x[0],this->mnx);
	                      
#endif
                        }
                        else 
                        {
#if defined (__AVX512F__)
	                       avx512_memcpy_unroll8x_ps(&this->mx[0],
                                                       &m_x[0],this->mnx);
#elif defined (__AVX__) || defined (__AVX2__) 
                             avx_memcpy_unroll8x_ps(&this->mx[0],
                                                    &m_x[0],this->mnx);
#elif defined (__SSE__)
                             sse_memcpy_unroll8x_ps(&this->mx[0],
                                                    &m_x[0],this->mnx);
#endif                             
                        }
                                 
                  }
                             
                   
                                           
                  float  *  begin() const 
                  {
                          
                          return (std::__addressof(this->mx[0]));
                  } 
                    
                  float  * end() const 
                  {
                           return (std::__addressof(this->mx[this->mnx]));
                  }

                  inline void swap(StatV1D_r4_t &other)
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
              };
              
              
            
           
              
              

              
              




} //gms

















#endif /*__GMS_STAT_CONTAINERS_H__*/
