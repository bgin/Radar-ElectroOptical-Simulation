

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_DERIV_AERODYN_ANGLES_R4_HPP__
#define __GMS_DERIV_AERODYN_ANGLES_R4_HPP__ 240420250951

namespace file_info {

     const unsigned int GMS_DERIV_AERODYN_ANGLES_R4_MAJOR = 1;
     const unsigned int GMS_DERIV_AERODYN_ANGLES_R4_MINOR = 1;
     const unsigned int GMS_DERIV_AERODYN_ANGLES_R4_MICRO = 0;
     const unsigned int GMS_DERIV_AERODYN_ANGLES_R4_FULLVER =
       1000U*GMS_DERIV_AERODYN_ANGLES_R4_MAJOR+100U*GMS_DERIV_AERODYN_ANGLES_R4_MINOR+
       10U*GMS_DERIV_AERODYN_ANGLES_R4_MICRO;
     const char * const GMS_DERIV_AERODYN_ANGLES_R4_CREATION_DATE = "24-04-2025 09:51  +00200 (THR 24 APR 2025 GMT+2)";
     const char * const GMS_DERIV_AERODYN_ANGLES_R4_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_DERIV_AERODYN_ANGLES_R4_SYNOPSIS      = "Dynamically allocated, 64-byte aligned 1s Derivative of Aerodynamic Angles (single-precision) type.";

}

#include <cstdint>
#include <cassert>
#include "GMS_config.h"
#include "GMS_malloc.h"
#if (USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"
#endif

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_DERIV_AERODYN_ANGLES_R4_NT_STORES)
#define USE_GMS_DERIV_AERODYN_ANGLES_R4_NT_STORES 0
#endif

namespace gms {

        namespace fdm {
               
                   /* Bernard Etkin, "Dynamics of Atmospheric Flight, formulae: 5.2.13, 5.2.14, page: 128"*/
                   struct __ATTR_ALIGN__(64) DerivAerodynAngles_r4_t final 
                   {
                          float * __restrict mdPw;
                          float * __restrict mdQw;
                          float * __restrict mdRw;
                          std::size_t        mn;
                          bool               mmalloc_type; // true stands for gms_mm_malloc, false stands for gms_tbb_malloc
                          bool               mfree_type;   // true stands for gms_mm_free,   false stands for gms_tbb_free
#if (USE_STRUCT_PADDING) == 1
                        PAD_TO(0,30)
#endif                 

                          DerivAerodynAngles_r4_t() = delete;

                          inline explicit DerivAerodynAngles_r4_t(const std::size_t n,
                                                                  const bool malloc_type,
                                                                  const bool free_type) noexcept(false)
                          {
                                 assert(n>0);
                                 this->mn = n;
                                 this->mmalloc_type = malloc_type;
                                 this->mfree_type   = free_type;
#if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_PROLOGE_BODY   
#endif 
                                   this->allocate();
#if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_EPILOGE_BODY
#endif 
                          }

                          inline DerivAerodynAngles_r4_t(const DerivAerodynAngles_r4_t &) = delete;

                          inline DerivAerodynAngles_r4_t(DerivAerodynAngles_r4_t && rhs) noexcept(true)
                          {
                                 this->mn    = rhs.mn;
                                 this->mmalloc_type = rhs.mmalloc_type;
                                 this->mfree_type   = rhs.mfree_type;
                                 this->mdPw = &rhs.mdPw[0];
                                 this->mdQw = &rhs.mdQw[0];
                                 this->mdRw = &rhs.mdRw[0];
                          }

                          inline ~DerivAerodynAngles_r4_t() noexcept(true)
                          {
                                 using namespace gms::common;
                                 if(this->mfree_type==true)
                                 {
                                     gms_mm_free(this->mdPw); this->mdPw = NULL;
                                     gms_mm_free(this->mdQw); this->mdQw = NULL;
                                     gms_mm_free(this->mdRw); this->mdRw = NULL;
                                 }
                                 else
                                 {
                                     gms_tbb_free(this->mdPw); this->mdPw = NULL;
                                     gms_tbb_free(this->mdQw); this->mdQw = NULL;
                                     gms_tbb_free(this->mdRw); this->mdRw = NULL;
                                 }
                          }

                          DerivAerodynAngles_r4_t & 
                          operator=(const DerivAerodynAngles_r4_t &) = delete;

                          
                          inline DerivAerodynAngles_r4_t &
                          operator=(DerivAerodynAngles_r4_t && rhs) noexcept(true)
                          {
                                using namespace gms::common;
                                if(this==&rhs) return (*this);
                                gms_swap(this->mn,   rhs.mn);
                                gms_swap(this->mmalloc_type, rhs.mmalloc_type);
                                gms_swap(this->mfree_type,   rhs.mfree_type);
                                gms_swap(this->mdPw, rhs.mdPw);
                                gms_swap(this->mdQw, rhs.mdQw);
                                gms_swap(this->mdRw, rhs.mdRw);

                                return (*this);
                          }

                          inline std::size_t size_mnbytes() const noexcept(true)
                          {
                               return static_cast<std::size_t>(this->mn*sizeof(float));
                          }

                          inline void allocate() noexcept(false) 
                          {
                               using namespace gms::common;
                               const std::size_t mnbytes{size_mnbytes()};
                               if(this->mmalloc_type==true)
                               {
                                     this->mdPw{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                                     this->mdQw{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                                     this->mdRw{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               }
                               else
                               {
                                     this->mdPw{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                                     this->mdQw{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                                     this->mdRw{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                               }
                          }

                        
                   };
        }
}


#endif /*__GMS_DERIV_AERODYN_ANGLES_R4_HPP__*/