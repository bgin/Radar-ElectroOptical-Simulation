

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

#ifndef __GMS_VEHICLE_EIA_R8_HPP__
#define __GMS_VEHICLE_EIA_R8_HPP__


namespace file_info {

     const unsigned int GMS_VEHICLE_EIA_R8_MAJOR = 1;
     const unsigned int GMS_VEHICLE_EIA_R8_MINOR = 1;
     const unsigned int GMS_VEHICLE_EIA_R8_MICRO = 0;
     const unsigned int GMS_VEHICLE_EIA_R8_FULLVER =
       1000U*GMS_VEHICLE_EIA_R8_MAJOR+100U*GMS_VEHICLE_EIA_R8_MINOR+
       10U*GMS_VEHICLE_EIA_R8_MICRO;
     const char * const GMS_VEHICLE_EIA_R8_CREATION_DATE = "26-04-2025 15:15  +00200 (SAT 24 APR 2025 GMT+2)";
     const char * const GMS_VEHICLE_EIA_R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_VEHICLE_EIA_R8_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Vehicle's Inertial Acceleration Equations (double-precision) type.";

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
#if !defined (USE_GMS_VEHICLE_EIA_R8_NT_STORES)
#define USE_GMS_VEHICLE_EIA_R8_NT_STORES 0
#endif


namespace gms 
{

       namespace fdm
       {
                  /*Bernard Etkin, "Dynamics of Atmospheric Flight", eq: 5.3.18, 5.3.22, page: 132-133*/
                struct __ATTR_ALIGN__(64) VehicleEIA_r8_t final 
                {
                       double * __restrict mAxv;  //acceleration component: xv.
                       double * __restrict mAyv;  //acceleration component: yv.
                       double * __restrict mAzv;  //acceleration component: zv 
                       std::size_t        mn;
                       bool               mmalloc_type; // true for gms_mm_malloc, false for gms_tbb_malloc
                       bool               mfree_type;   // true for gms_mm_free,   false for gms_tbb_free
#if (USE_STRUCT_PADDING) == 1
                        PAD_TO(0,30)
#endif                      

                        VehicleEIA_r8_t() = delete;

                        inline VehicleEIA_r8_t(const std::size_t n,
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

                        VehicleEIA_r8_t(const VehicleEIA_r8_t &) = delete;

                        inline VehicleEIA_r8_t(VehicleEIA_r8_t && rhs) noexcept(true)
                        {
                               this->mn   = rhs.mn;
                               this->mmalloc_type = rhs.mmalloc_type;
                               this->mfree_type   = rhs.mfree_type;
                               this->mAxv = &rhs.mAxv[0];
                               this->mAyv = &rhs.mAyv[0];
                               this->mAzv = &rhs.mAzv[0];
                        }

                        inline ~VehicleEIA_r8_t() noexcept(true)
                        {
                               using namespace gms::common;
                               if(this->mfree_type==true)
                               {
                                  gms_mm_free(this->mAzv); this->mAzv = NULL;
                                  gms_mm_free(this->mAyv); this->mAyv = NULL;
                                  gms_mm_free(this->mAzv); this->mAzv = NULL;
                               }
                               else 
                               {
                                  gms_tbb_free(this->mAzv); this->mAzv = NULL;
                                  gms_tbb_free(this->mAyv); this->mAyv = NULL;
                                  gms_tbb_free(this->mAzv); this->mAzv = NULL;
                               }
                               this->mn = 0ULL;
                        }

                        VehicleEIA_r8_t & operator=(const VehicleEIA_r8_t &) = delete;

                        inline VehicleEIA_r8_t & operator=(VehicleEIA_r8_t && rhs) noexcept(true)
                        {
                               using namespace gms::common;
                               if(this==&rhs) return (*this);

                               gms_swap(this->mn,  rhs.mn);
                               gms_swap(this->mmalloc_type, rhs.mmalloc_type);
                               gms_swap(this->mfree_type,   rhs.mfree_type);
                               gms_swap(this->mAxv,rhs.mAxv);
                               gms_swap(this->mAyv,rhs.mAyv);
                               gms_swap(this->mAzv,rhs.mAzv);
                               
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
                                  this->mAxv{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                                  this->mAyv{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                                  this->mAzv{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               }
                               else 
                               {
                                  this->mAxv{reinterpret_cast<double * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                                  this->mAyv{reinterpret_cast<double * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                                  this->mAzv{reinterpret_cast<double * __restrict>(gms_tbb_malloc(mnbytes,64ULL))};
                               }
                
                        }
                        
                };
       } //fdm
       
} // gms









#endif /*__GMS_VEHICLE_EIA_R8_HPP__*/