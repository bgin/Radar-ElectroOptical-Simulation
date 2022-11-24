
#ifndef __GMS_ANTENNA_TYPES_V2_H__
#define __GMS_ANTENNA_TYPES_V2_H__


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

                   using vr8 = std::valarray<double>;
                   using vr4 = std::valarray<float>;
                   using vc4 = std::valarray<std::complex<float>>;
                   using vc8 = std::valarray<std::complex<double>>;

                   // col: 20
                   struct alignas(64) E_c4_t {

                          public:
                          
                          int32_t        m_npts;
                          vc4            m_ex;
                          vc4            m_ey;
                          vc4            m_ez;

                          E_c4_t()      = default;

                          E_c4_t(const int32_t);

                          E_c4_t(const E_c4_t &);

                          E_c4_t(E_c4_t &&) noexcept(true);

                         ~E_c4_t()     = default;

                          E_c4_t & operator=(const E_c4_t &);

                          E_c4_t & operator=(E_c4_t &&) noexcept(true);
               };
               // col: 16

                 
       } // radiolocation

} // gms








#endif /*__GMS_ANTENNA_TYPES_V2_H__*/
