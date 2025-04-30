
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

#ifndef __GMS_EARTH_MOTION_R4_HPP__
#define __GMS_EARTH_MOTION_R4_HPP__ 200420250959


namespace file_info {

     const unsigned int GMS_EARTH_MOTION_R4_MAJOR = 1;
     const unsigned int GMS_EARTH_MOTION_R4_MINOR = 1;
     const unsigned int GMS_EARTH_MOTION_R4_MICRO = 0;
     const unsigned int GMS_EARTH_MOTION_R4_FULLVER =
       1000U*GMS_EARTH_MOTION_R4_MAJOR+100U*GMS_EARTH_MOTION_R4_MINOR+
       10U*GMS_EARTH_MOTION_R4_MICRO;
     const char * const GMS_EARTH_MOTION_R4_CREATION_DATE = "20-04-2025 09:59 AM +00200 (SUN 20 APR 2025 GMT+2)";
     const char * const GMS_EARTH_MOTION_R4_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_EARTH_MOTION_R4_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Earth Motion (single-precision) type.";

}

#include <cstdint>
#include <immintrin.h>
#include <cstdlib>
#include <cassert>
#include "GMS_config.h"
#include "GMS_malloc.h"
#if (USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"
#endif

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_EARTH_MOTION_R4_NT_STORES)
#define USE_GMS_EARTH_MOTION_R4_NT_STORES 0
#endif

namespace gms {

    namespace fdm {

               /*
                    Earth Motion angular velocity
                    Bernard Etkin, "Dynamics of Atmospheric Flight", page: 124.
               */
              struct __ATTR_ALIGN__(64) EarthMotion_r4_t 
              {
                     // Only non-zero vector-elements are defined
                     float * __restrict mOE; //rate of rotation 7.27*10^-5 rad/sec
                     float * __restrict mcosLat; // cosine of Latitude
                     float * __restrict msinLat; // sine of Latitude
                     float * __restrict mcosLon; // cosine of Longitude
                     float * __restrict msinLon; // sine of Longtitude
                     std::size_t        mnsec;   // number of seconds of Earth's rotation
                     bool               mmalloc_type; // true for gms_mm_malloc, false for gms_tbb_malloc
                     bool               mfree_type;   // true for gms_mm_free,   false for gms_tbb_free
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,14)
#endif
                     EarthMotion_r4_t() = delete;

                     inline  EarthMotion_r4_t(const std::size_t nsec,
                                              const bool        malloc_type,
                                              const bool        free_type) noexcept(false)
                     {
                          assert(nsec>0ULL);
                          this->mnsec = nsec;
                          this->mmalloc_type = malloc_type;
                          this->mfree_type   = free_type;
#if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_PROLOGE_BODY   
#endif                           
                          allocate();
#if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_EPILOGE_BODY
#endif
                     }

                     EarthMotion_r4_t(const EarthMotion_r4_t &) = delete;

                     inline EarthMotion_r4_t(EarthMotion_r4_t && rhs) noexcept(true)
                     {
                          this->mnsec        = rhs.mnsec;
                          this->mmalloc_type = rhs.mmalloc_type;
                          this->mfree_type   = rhs.mfree_type;
                          this->mOE          = &rhs.mOE[0];
                          this->mcosLat      = &rhs.mcosLat[0];
                          this->msinLat      = &rhs.msinLat[0];
                          this->mcosLon      = &rhs.mcosLon[0];
                          this->msinLon      = &rhs.msinLon[0];
                     }

                     inline ~EarthMotion_r4_t() noexcept(true)
                     {
                         using namespace gms::common;
                         if(this->mfree_type==true)
                         {
                             gms_mm_free(this->msinLon); this->msinLon = NULL;
                             gms_mm_free(this->mcosLon); this->mcosLon = NULL;
                             gms_mm_free(this->msinLat); this->msinLat = NULL;
                             gms_mm_free(this->mcosLat); this->mcosLat = NULL;
                             gms_mm_free(this->mOE);     this->mOe     = NULL;
                         }
                         else 
                         {
                             gms_tbb_free(this->msinLon); this->msinLon = NULL;
                             gms_tbb_free(this->mcosLon); this->mcosLon = NULL;
                             gms_tbb_free(this->msinLat); this->msinLat = NULL;
                             gms_tbb_free(this->mcosLat); this->mcosLat = NULL;
                             gms_tbb_free(this->mOE);     this->mOe     = NULL;
                         }
                     }

                     EarthMotion_r4_t & operator=(const EarthMotion_r4_t &) = delete;

                     EarthMotion_r4_t & operator=(EarthMotion_r4_t && rhs) noexcept(true)
                     {
                          using namespace gms::commom;
                          if(this==&rhs) return (*this)
                          gms_swap(this->mnsec,  rhs.mnsec);
                          gms_swap(this->mmalloc_type, rhs.mmalloc_type);
                          gms_swap(this->mfree_type,   rhs.mfree_type);
                          gms_swap(this->mOE,    rhs.mOE);
                          gms_swap(this->mcosLat,rhs.mcosLat);
                          gms_swap(this->msinLat,rhs.msinLat);
                          gms_swap(this->mcosLat,rhs.mcosLat);
                          return (*this);
                     }

                     inline std::size_t size_nsecbytes() const noexcept(true)
                     {
                         return static_cast<std::size_t>(this->mnsec*sizeof(float));
                     }
                     
                     inline void allocate() noexcept(false)
                     {
                          using namespace gms::common;
                          const std::size_t mnsecbytes{size_nsecbytes()};
                          if(this->mmalloc_type==true)
                          {
                              this->mOE{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnsecbytes,64ULL))};
                              this->mcosLat{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnsecbytes,64ULL))};
                              this->msinLat{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnsecbytes,64ULL))};
                              this->mcosLon{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnsecbytes,64ULL))};
                              this->msinLon{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnsecbytes,64ULL))};
                          }
                          else 
                          {
                              this->mOE{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnsecbytes,64ULL))};
                              this->mcosLat{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnsecbytes,64ULL))};
                              this->msinLat{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnsecbytes,64ULL))};
                              this->mcosLon{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnsecbytes,64ULL))};
                              this->msinLon{reinterpret_cast<float * __restrict>(gms_tbb_malloc(mnsecbytes,64ULL))};
                          }
                     }
                     
                     
                     
              };
    }
}











#endif /*__GMS_EARTH_MOTION_R4_HPP__*/