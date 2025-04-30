

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

#ifndef __GMS_DERIV_EULER_ANGLES_R8_HPP__
#define __GMS_DERIV_EULER_ANGLES_R8_HPP__ 240420250910

namespace file_info {

     const unsigned int GMS_DERIV_EULER_ANGLES_R8_MAJOR = 1;
     const unsigned int GMS_DERIV_EULER_ANGLES_R8_MINOR = 1;
     const unsigned int GMS_DERIV_EULER_ANGLES_R8_MICRO = 0;
     const unsigned int GMS_DERIV_EULER_ANGLES_R8_FULLVER =
       1000U*GMS_DERIV_EULER_ANGLES_R8_MAJOR+100U*GMS_DERIV_EULER_ANGLES_R8_MINOR+
       10U*GMS_DERIV_EULER_ANGLES_R8_MICRO;
     const char * const GMS_DERIV_EULER_ANGLES_R8_CREATION_DATE = "24-04-2025 09:10  +00200 (THR 24 APR 2025 GMT+2)";
     const char * const GMS_DERIV_EULER_ANGLES_R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_DERIV_EULER_ANGLES_R8_SYNOPSIS      = "Dynamically allocated, 64-byte aligned 1s Derivative of Euler Angles (double-precision) type.";

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
#if !defined (USE_GMS_DERIV_EULER_ANGLES_R8_NT_STORES)
#define USE_GMS_DERIV_EULER_ANGLES_R8_NT_STORES 0
#endif

namespace gms {

        namespace fdm {
               
                   struct __ATTR_ALIGN__(64) DerivEulerAngles_r8_t final 
                   {
                          double * __restrict mdPhi;
                          double * __restrict mdTht;
                          double * __restrict mdPsi;
                          std::size_t        mn;
                          bool               mmalloc_type; // true stands for gms_mm_malloc, false stands for gms_tbb_malloc
                          bool               mfree_type;   // true stands for gms_mm_free,   false stands for gms_tbb_free
 #if (USE_STRUCT_PADDING) == 1
                        PAD_TO(0,30)
#endif                 

                          DerivEulerAngles_r8_t() = delete;

                          inline  DerivEulerAngles_r8_t(const std::size_t n,
                                                        const bool malloc_type,
                                                        const bool free_type) noexcept(false)
                          {
                                 assert(n>0ULL);
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

                          inline DerivEulerAngles_r8_t(const DerivEulerAngles_r8_t &) = delete;

                          inline DerivEulerAngles_r8_t(DerivEulerAngles_r8_t && rhs) noexcept(true)
                          {
                                 this->mn    = rhs.mn;
                                 this->mmalloc_type = rhs.mmalloc_type;
                                 this->mfree_type   = rhs.mfree_type;
                                 this->mdPhi = &rhs.mdPhi[0];
                                 this->mdTht = &rhs.mdTht[0];
                                 this->mdPsi = &rhs.mdPsi[0];
                          }

                          inline ~DerivEulerAngles_r8_t() noexcept(true)
                          {
                                 using namespace gms::common;
                                 if(this->mfree_type==true)
                                 {
                                      gms_mm_free(this->mdPsi); this->mdPsi = NULL;
                                      gms_mm_free(this->mdTht); this->mdTht = NULL;
                                      gms_mm_free(this->mdPhi); this->mdPhi = NULL;
                                 }
                                 else 
                                 {
                                      gms_tbb_free(this->mdPsi); this->mdPsi = NULL;
                                      gms_tbb_free(this->mdTht); this->mdTht = NULL;
                                      gms_tbb_free(this->mdPhi); this->mdPhi = NULL;
                                 }
                          }

                          DerivEulerAngles_r8_t & operator=(const DerivEulerAngles_r4_t &) = delete;

                          inline DerivEulerAngles_r8_t & operator=(DerivEulerAngles_r4_t && rhs) noexcept(true)
                          {
                                 using namespace gms::common;
                                 if(this==&rhs) return (*this);
                                 gms_swap(this->mn,    rhs.mn);
                                 gms_swap(this->mmalloc_type, rhs.mmalloc_type);
                                 gms_swap(this->mfree_type,   rhs.mfree_type);
                                 gms_swap(this->mdPhi, rhs.mdPhi);
                                 gms_swap(this->mdTht, rhs.mdTht);
                                 gms_swap(this->mdPsi, rhs.mdPsi);

                                 return (*this);
                          }


                          inline std::size_t size_mnbytes() const noexcept(true)
                          {
                               return static_cast<std::size_t>(this->mn*sizeof(double));
                          }

                          inline void allocate() noexcept(false) 
                          {
                               using namespace gms::common;
                               const std::size_t mnbytes{size_mnbytes()};
                               if(this->mmalloc_type==true)
                               {
                                   this->mdPhi{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                                   this->mdTht{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                                   this->mdPsi{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               }
                               else 
                               {
                                   this->mdPhi{reinterpret_cast<double * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                                   this->mdTht{reinterpret_cast<double * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                                   this->mdPsi{reinterpret_cast<double * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                               }
                          }

                        
                   };
        }
}


#endif /*__GMS_DERIV_EULER_ANGLES_R8_HPP__*/