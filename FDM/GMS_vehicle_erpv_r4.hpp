

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

#ifndef __GMS_VEHICLE_ERPV_R4_HPP__
#define __GMS_VEHICLE_ERPV_R4_HPP__ 240420250951


namespace file_info {

     const unsigned int GMS_VEHICLE_ERPV_R4_MAJOR = 1;
     const unsigned int GMS_VEHICLE_ERPV_R4_MINOR = 1;
     const unsigned int GMS_VEHICLE_ERPV_R4_MICRO = 0;
     const unsigned int GMS_VEHICLE_ERPV_R4_FULLVER =
       1000U*GMS_VEHICLE_ERPV_R4_MAJOR+100U*GMS_VEHICLE_ERPV_R4_MINOR+
       10U*GMS_VEHICLE_ERPV_R4_MICRO;
     const char * const GMS_VEHICLE_ERPV_R4_CREATION_DATE = "24-04-2025 09:51  +00200 (THR 24 APR 2025 GMT+2)";
     const char * const GMS_VEHICLE_ERPV_R4_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_VEHICLE_ERPV_R4_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Vehicle's Earth Related Position and Velocity (single-precision) type.";

}

#include <cstdint>
#include <immintrin.h>
#include <cstdlib>
#include "GMS_config.h"
#include "GMS_malloc.h"
#if (USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"
#endif

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_VEHICLE_ERPV_R4_NT_STORES)
#define USE_GMS_VEHICLE_ERPV_R4_NT_STORES 0
#endif

namespace gms {

        namespace fdm {

                      /* Bernard Etkin, "Dynamics of Atmospheric Flight, formulae: 5.3.1, 5.3.4 page: 129"*/
                     struct __ATTR_ALIGN__(64) VehicleERPV_r4_t final 
                     {
                            float * __restrict mdR;    //  1st derivative of geocentric radius
                            float * __restrict mdLon;   // 1st derivative of longtitude, i.e. mu
                            float * __restrict mdLat;   // 1st derivative of latitude,   i.e. lambda
                            float * __restrict mVxv;    // velocity component at xv.
                            float * __restrict mVyv;    // velocity component at yv.
                            float * __restrict mVzv;    // velocity component at zv 
                            std::size_t        mn;
#if (USE_STRUCT_PADDING) == 1
                        PAD_TO(0,8)
#endif                           
                                                  
                            VehicleERPV_r4_t() = delete;

                            inline explicit (const std::size_t n) noexcept(true)
                            {
                                   assert(n>0);
                                   this->mn = n;
#if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_PROLOGE_BODY   
#endif                                
                                   this->allocate();
#if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_EPILOGE_BODY
#endif 
                            }   

                            VehicleERPV_r4_t(const VehicleERPV_r4_t &) = delete;

                            inline VehicleERPV_r4_t(VehicleERPV_r4_t && rhs) noexcept(true)
                            {
                                   this->mn     = rhs.mn;
                                   this->mdR    = &rhs.mdR[0];
                                   this->mdLon  = &rhs.mdLon[0];
                                   this->mdLat  = &rhs.mdLat[0];
                                   this->mVxv   = &rhs.mVxv[0];
                                   this->mVyv   = &rhs.mVyv[0];
                                   this->mVzv   = &rhs.mVzv[0];
                                   
                            } 

                            inline ~VehicleERPV_r4_t() noexcept(true)
                            {
                                using namespace gms::common;
                                gms_mm_free(this->mVzv);  this->mVzv  = NULL;
                                gms_mm_free(this->mVyv);  this->mVyv  = NULL;
                                gms_mm_free(this->mVxv);  this->mVxv  = NULL;
                                gms_mm_free(this->mdLat); this->mdLat = NULL;
                                gms_mm_free(this->mdLon); this->mdLon = NULL;
                                gms_mm_free(this->mdR);   this->mdR   = NULL;
                                this->mn = 0ULL;
                            }

                            VehicleERPV_r4_t & operator=(const VehicleERPV_r4_t &) = delete;

                            inline VehicleERPV_r4_t & operator=(VehicleERPV_r4_t && rhs) noexcept(true)
                            {
                                  using namespace gms::common;
                                  if(this==&rhs) return (*this);

                                  gms_swap(this->mn,    rhs.mn);
                                  gms_swap(this->mdR,   rhs.mdR);
                                  gms_swap(this->mdLon, rhs.mdLon);
                                  gms_swap(this->mdLat, rhs.mdLat);
                                  gms_swap(this->mVxv,  rhs.mVxv);
                                  gms_swap(this->mVyv,  rhs.mVyv);
                                  gms_swap(this->mVzv,  rhs.mVzv);
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
                               this->mdR{  reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mdLon{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mdLat{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mVxv{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mVyv{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mVzv{reinterpret_cast<float * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                            } 

                          
                };
        }
}












#endif /*__GMS_VEHICLE_ERPV_R4_HPP__*/