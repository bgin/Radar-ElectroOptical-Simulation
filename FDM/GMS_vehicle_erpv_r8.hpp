


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

#ifndef __GMS_VEHICLE_ERPV_R8_HPP__
#define __GMS_VEHICLE_ERPV_R8_HPP__ 260420250647


namespace file_info {

     const unsigned int GMS_VEHICLE_ERPV_R8_MAJOR = 1;
     const unsigned int GMS_VEHICLE_ERPV_R8_MINOR = 1;
     const unsigned int GMS_VEHICLE_ERPV_R8_MICRO = 0;
     const unsigned int GMS_VEHICLE_ERPV_R8_FULLVER =
       1000U*GMS_VEHICLE_ERPV_R8_MAJOR+100U*GMS_VEHICLE_ERPV_R8_MINOR+
       10U*GMS_VEHICLE_ERPV_R8_MICRO;
     const char * const GMS_VEHICLE_ERPV_R8_CREATION_DATE = "26-04-2025 06:47  +00200 (SAT 26 APR 2025 GMT+2)";
     const char * const GMS_VEHICLE_ERPV_R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_VEHICLE_ERPV_R8_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Vehicle's Earth Related Position and Velocity (double-precision) type.";

}

#include <cstdint>
#include <immintrin.h>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "GMS_config.h"
#include "GMS_malloc.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_VEHICLE_ERPV_R8_NT_STORES)
#define USE_GMS_VEHICLE_ERPV_R8_NT_STORES 0
#endif

namespace gms {

        namespace fdm {

                      /* Bernard Etkin, "Dynamics of Atmospheric Flight, formulae: 5.3.1, page: 129"*/
                     struct __ATTR_ALIGN__(64) VehicleERPV_r8_t 
                     {
                            double * __restrict mdR;    //  1st derivative of geocentric radius
                            double * __restrict mdLon;   // 1st derivative of longtitude, i.e. mu
                            double * __restrict mdLat;   // 1st derivative of latitude,   i.e. lambda
                            std::size_t        mn;
#if (USE_STRUCT_PADDING) == 1
                        PAD_TO(0,32)
#endif                        
                            VehicleERPV_r8_t() = delete;

                            inline explicit (const std::size_t n) noexcept(true)
                            {
                                   assert(n>0);
                                   this->mn = n;
                                   this->allocate();
                            }   

                            VehicleERPV_r8_t(const VehicleERPV_r8_t &) = delete;

                            inline VehicleERPV_r8_t(VehicleERPV_r8_t && rhs) noexcept(true)
                            {
                                   this->mn     = rhs.mn;
                                   this->mdR    = &rhs.mdR[0];
                                   this->mdLon  = &rhs.mdLon[0];
                                   this->mdLat  = &rhs.mdLat[0];
                            } 

                            inline ~VehicleERPV_r8_t() noexcept(true)
                            {
                                using namespace gms::common;
                                gms_mm_free(this->mdLat); this->mdLat = NULL;
                                gms_mm_free(this->mdLon); this->mdLon = NULL;
                                gms_mm_free(this->mdR);   this->mdR   = NULL;
                            }

                            VehicleERPV_r8_t & operator=(const VehicleERPV_r8_t &) = delete;

                            inline VehicleERPV_r8_t & operator=(VehicleERPV_r8_t && rhs) noexcept(true)
                            {
                                  using namespace gms::common;
                                  if(this==&rhs) return (*this);
                                  gms_swap(this->mn,    rhs.mn);
                                  gms_swap(this->mdR,   rhs.mdR);
                                  gms_swap(this->mdLon, rhs.mdLon);
                                  gms_swap(this->mdLat, rhs.mdLat);

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
                               this->mdR{  reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mdLon{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mdLat{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                            } 

                           bool dump_state_to_file(const char * fname) const noexcept(false)
                           {
                              assert(fname != nullptr);
                              File * fp = nullptr;
                              bool b_result;

                              fp = std::fopen(fname,"w");
                              if(!fp)
                              {
                                  std::printf("[%s]: failed to open a file: %s\n",__PRETTY_FUNCTION__,fname);
                                  b_result = false;
                                  return b_result;
                              }

                              fprintf(fp,"Writing content of array mdPw[%d] to file: %s, status: started\n\n",this->mn,fname);
                              for(std::size_t i = 0ULL; i != this->mn; ++i)
                              {
                                  fprintf(fp, "mdR[%llu,%d] = %2.9f\n",i,this->mn,this->mdR[i]);
                              }
                              fprintf(fp,"\n\n");
                              
                              fprintf(fp,"Writing content of array mdLon[%d] to file: %s, status: started\n\n",this->mn,fname);
                              for(std::size_t i = 0ULL; i != this->mn; ++i)
                              {
                                  fprintf(fp, "mdLon[%llu,%d] = %2.9f\n",i,this->mn,this->mdLon[i]);
                              }
                              fprintf(fp,"\n\n");

                              fprintf(fp,"Writing content of array mdLat[%d] to file: %s, status: started\n\n",this->mn,fname);
                              for(std::size_t i = 0ULL; i != this->mn; ++i)
                              {
                                  fprintf(fp, "mdLat[%llu,%d] = %2.9f\n",i,this->mn,this->mdLat[i]);
                              }
                              fprintf(fp,"\n\n");

                              fprintf(fp, "Finished writing a content to file: %s\n", fname);
                              fclose(fp);
                              b_result = true;
                              return b_result;
                         }             
                };
        }
}












#endif /*__GMS_VEHICLE_ERPV_R8_HPP__*/