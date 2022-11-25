
#ifndef __GMS_ANTENNA_TYPES_V2_HPP__
#define __GMS_ANTENNA_TYPES_V2_HPP__


namespace file_info {

     const unsigned int GMS_ANTENNA_TYPES_V2_MAJOR = 1;
     const unsigned int GMS_ANTENNA_TYPES_V2_MINOR = 1;
     const unsigned int GMS_ANTENNA_TYPES_V2_MICRO = 0;
     const unsigned int GMS_ANTENNA_TYPES_V2_FULLVER =
       1000U*GMS_ANTENNA_TYPES_V2_MAJOR+100U*GMS_ANTENNA_TYPES_V2_MINOR+
       10U*GMS_ANTENNA_TYPES_V2_MICRO;
     const char * const GMS_ANTENNA_TYPES_V2_CREATION_DATE = "24-11-2022 16:08 +00200 (THR 24 NOV 2022 GMT+2)";
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

#if defined (__INTEL_COMPILER)
#include <../perf_headers/c++/valarray>
#else
#include <valarray>
#endif
#include <cstdint>
#include <complex>
#include "GMS_config"


namespace gms {

         namespace radiolocation {

                  template<typename T>
                  struct alignas(64) Exyz_t {

                     int32_t            m_npts;
                     std::valarray<T>   m_ex;
                     std::valarray<T>   m_ey;
                     std::valarray<T>   m_ez;

                     Exyz_t()          = default;
 
                     Exyz_t(const int32_t npts) {

                           m_npts  = npts;
                           m_ex    = std::valarray<T>(m_npts);
                           m_ey    = std::valarray<T>(m_npts);
                           m_ez    = std::valarray<T>(m_npts);
                     }

                     Exyz_t(const Exyz_t &x) {

                           m_npts   = x.m_npts;
                           m_ex     = x.m_ex;
                           m_ey     = x.m_ey;
                           m_ez     = x.m_ez;
                     }

                     Exyz_t(Exyz_t &&x) noexcept(true) {

                          m_npts     = std::move(m_npts);
                          m_ex       = std::move(x.m_ex);
                          m_ey       = std::move(x.m_ey);
                          m_ez       = std::move(x.m_ez); 
                     }

                    ~Exyz()             = default;

                     Exyz_t & operator=(const Exyz_t &x) {

                          if(this == &x) return (*this);
                          Exyz_t<T> tmp(x);
                          std::swap(*this,tmp);
                          return (*this);
                     }

                     Exyz_t & operator=(Exyz_t &&x) noexcept(true) {
                          
                          if(this == &x) return (*this);
                          *this = std::move(x);
                          return (*this);
                     }
                     
               };

                  
                  template<typename T>
                  struct alignas(64) Hxyz_t {

                     int32_t            m_npts;
                     std::valarray<T>   m_hx;
                     std::valarray<T>   m_hy;
                     std::valarray<T>   m_hz;

                     Hxyz_t()          = default;
 
                     Hxyz_t(const int32_t npts) {

                           m_npts  = npts;
                           m_hx    = std::valarray<T>(m_npts);
                           m_hy    = std::valarray<T>(m_npts);
                           m_hz    = std::valarray<T>(m_npts);
                     }

                     Hxyz_t(const Hxyz_t &x) {

                           m_npts   = x.m_npts;
                           m_hx     = x.m_hx;
                           m_hy     = x.m_hy;
                           m_hz     = x.m_hz;
                     }

                     Hxyz_t(Hxyz_t &&x) noexcept(true) {

                          m_npts     = std::move(m_npts);
                          m_hx       = std::move(x.m_hx);
                          m_hy       = std::move(x.m_hy);
                          m_hz       = std::move(x.m_hz); 
                     }

                    ~Hxyz()             = default;

                     Hxyz_t & operator=(const Hxyz_t &x) {

                          if(this == &x) return (*this);
                          Hxyz_t<T> tmp(x);
                          std::swap(*this,tmp);
                          return (*this);
                     }

                     Hxyz_t & operator=(Hxyz_t &&x) noexcept(true) {
                          
                          if(this == &x) return (*this);
                          *this = std::move(x);
                          return (*this);
                     }
                     
               };
 


                  /*
                      // Complex magnetic current
                   */  
                  template<typename T>
                  struct alignas(64) JMxyz_t {

                     int32_t            m_npts;
                     std::valarray<T>   m_jmx;
                     std::valarray<T>   m_jmy;
                     std::valarray<T>   m_jmz;

                     JMxyz_t()          = default;
 
                     JMxyz_t(const int32_t npts) {

                           m_npts   = npts;
                           m_jmx    = std::valarray<T>(m_npts);
                           m_jmy    = std::valarray<T>(m_npts);
                           m_jmz    = std::valarray<T>(m_npts);
                     }

                     JMxyz_t(const JMxyz_t &x) {

                           m_npts    = x.m_npts;
                           m_jmx     = x.m_jmx;
                           m_jmy     = x.m_jmy;
                           m_jmz     = x.m_jmz;
                     }

                     JMxyz_t(JMxyz_t &&x) noexcept(true) {

                          m_npts      = std::move(m_npts);
                          m_jmx       = std::move(x.m_jmx);
                          m_jmy       = std::move(x.m_jmy);
                          m_jmz       = std::move(x.m_jmz); 
                     }

                    ~JMxyz()             = default;

                     JMxyz_t & operator=(const JMxyz_t &x) {

                          if(this == &x) return (*this);
                          JMxyz_t<T> tmp(x);
                          std::swap(*this,tmp);
                          return (*this);
                     }

                     JMxyz_t & operator=(JMxyz_t &&x) noexcept(true) {
                          
                          if(this == &x) return (*this);
                          *this = std::move(x);
                          return (*this);
                     }
                     
               }; 


                    /*
                      // Complex electric current
                   */  
                  template<typename T>
                  struct alignas(64) JExyz_t {

                     int32_t            m_npts;
                     std::valarray<T>   m_jex;
                     std::valarray<T>   m_jey;
                     std::valarray<T>   m_jez;

                     JExyz_t()          = default;
 
                     JExyz_t(const int32_t npts) {

                           m_npts   = npts;
                           m_jex    = std::valarray<T>(m_npts);
                           m_jey    = std::valarray<T>(m_npts);
                           m_jez    = std::valarray<T>(m_npts);
                     }

                     JExyz_t(const JExyz_t &x) {

                           m_npts    = x.m_npts;
                           m_jex     = x.m_jex;
                           m_jey     = x.m_jey;
                           m_jez     = x.m_jez;
                     }

                     JExyz_t(JExyz_t &&x) noexcept(true) {

                          m_npts      = std::move(m_npts);
                          m_jex       = std::move(x.m_jex);
                          m_jey       = std::move(x.m_jey);
                          m_jez       = std::move(x.m_jez); 
                     }

                    ~JExyz()             = default;

                     JExyz_t & operator=(const JExyz_t &x) {

                          if(this == &x) return (*this);
                          JExyz_t<T> tmp(x);
                          std::swap(*this,tmp);
                          return (*this);
                     }

                     JExyz_t & operator=(JExyz_t &&x) noexcept(true) {
                          
                          if(this == &x) return (*this);
                          *this = std::move(x);
                          return (*this);
                     }
                     
               }; 


                  /*
                        // Time-Harmonic complex exponential
                   */
                  template<typename T1, typename T2>
                  struct alignas(64) eikr_t {

                         int32_t             m_npts;
                         T1                  m_k;
                         std::valarray<T1>   m_R;
                         std::valarray<T2>   m_ce;

                         eikr_t()            = default;
                         eikr_t(const int32_t npts,
                                const T1      k) {

                              m_npts   = npts;
                              m_k      = k;
                              m_R      = std::valarray<T1>(m_npts);
                              m_ce     = std::valarray<T2>(m_npts);
                         }

                         eikr_t(const eikr_t &x) {
                               
                              m_npts   = x.m_npts;
                              m_k      = x.m_k;
                              m_R      = x.m_R;
                              m_ce     = x.m_ce;
                         }

                         eikr_t(eikr_t &&x) noexcept(true) {
                             
                              m_npts    = std::move(x.m_npts);
                              m_k       = std::move(x.m_k);
                              m_R       = std::move(x.m_R);
                              m_ce      = std::move(x.m_ce);
                         }

                        ~eikr_t()              = default;

                         eikr_t & operator=(const eikr_t &x) {
                               
                              if(this == &x) return (*this);
                              eikr_t<T1,T2> tmp(x);
                              std::swap(*this,tmp);
                              return (*this);
                         }

                         eikr_t & operator=(eikr_t &&x) noexcept(true) {
                              
                              if(this == &x) return (*this);
                              *this = std::move(x);
                              return (*this);
                         }
               };


                  /*
                       // ! Formula (1-37)
                       //! Average level of side lobes
                   */
                  template<typename T>
                  struct alignas(64) f137_t {
                         
                       int32_t          m_nth;
                       int32_t          m_nph;
                       T                m_ith;
                       T                m_iph;
                       T                m_ifac;
                       T                m_omega;
                       T                m_avsl;
                       std::valarray<T> m_sinth;
                       std::valarray<T> m_F;

                       f137_t()         = default;

                       f137_t(const int32_t nth,
                              const int32_t nph,
                              const T       ith,
                              const T       iph,
                              const T       ifac,
                              const T       omega) {
                           
                           m_nth       = nth;
                           m_nph       = nph;
                           m_ith       = ith;
                           m_iph       = iph;
                           m_ifac      = ifac;
                           m_omega     = omega;
                           m_avsl      = static_cast<T>(0.0);
                           m_sinth     = std::valarray<T>(m_nth);
                           m_F         = std::valarray<T>(m_nth*m_nph);
                      }


                      f137_t(const f137_t &x) {
                           
                           m_nth       = x.m_nth;
                           m_nph       = x.m_nph;
                           m_ith       = x.m_ith;
                           m_iph       = x.m_iph;
                           m_ifac      = x.m_ifac;
                           m_omega     = x.m_omega;
                           m_avsl      = x.m_avsl;
                           m_sinth     = x.m_sinth;
                           m_F         = x.m_F;
                      }

                      f137_t(f137_t &&x) noexcept(true) {
                           
                           m_nth       = std::move(x.m_nth);
                           m_nph       = std::move(x.m_nph);
                           m_ith       = std::move(x.m_ith);
                           m_iph       = std::move(x.m_iph);
                           m_ifac      = std::move(x.m_ifac);
                           m_omega     = std::move(x.m_omega);
                           m_avsl      = std::move(x.m_avsl);
                           m_sinth     = std::move(x.m_sinth);
                           m_F         = std::move(x.m_F);
                      }

                     ~f137_t()                  = default;

                      f137_t & operator=(const f137_t &x) {
                          
                           if(this == &x) return (*this);
                           f137_t<T> tmp(x);
                           std::swap(*this,tmp);
                           return (*this);
                      }

                      f137_t & operator=(f137_t &&x) noexcept(true) {
                           
                           if(this == &x) return (*this);
                           *this = std::move(x);
                            return (*this);
                      }
               };
               

                 
       } // radiolocation

} // gms












#endif /*__GMS_ANTENNA_TYPES_V2_HPP__*/
