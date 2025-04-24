


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

#ifndef __GMS_DERIV_AERODYN_ANGLES_R8_HPP__
#define __GMS_DERIV_AERODYN_ANGLES_R8_HPP__ 240420251031

namespace file_info {

     const unsigned int GMS_DERIV_AERODYN_ANGLES_R8_MAJOR = 1;
     const unsigned int GMS_DERIV_AERODYN_ANGLES_R8_MINOR = 1;
     const unsigned int GMS_DERIV_AERODYN_ANGLES_R8_MICRO = 0;
     const unsigned int GMS_DERIV_AERODYN_ANGLES_R8_FULLVER =
       1000U*GMS_DERIV_AERODYN_ANGLES_R8_MAJOR+100U*GMS_DERIV_AERODYN_ANGLES_R8_MINOR+
       10U*GMS_DERIV_AERODYN_ANGLES_R8_MICRO;
     const char * const GMS_DERIV_AERODYN_ANGLES_R8_CREATION_DATE = "24-04-2025 10:31  +00200 (THR 24 APR 2025 GMT+2)";
     const char * const GMS_DERIV_AERODYN_ANGLES_R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_DERIV_AERODYN_ANGLES_R8_SYNOPSIS      = "Dynamically allocated, 64-byte aligned 1s Derivative of Aerodynamic Angles (double-precision) type.";

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
#if !defined (USE_GMS_DERIV_AERODYN_ANGLES_R8_NT_STORES)
#define USE_GMS_DERIV_AERODYN_ANGLES_R8_NT_STORES 0
#endif

namespace gms {

        namespace fdm {
               
                   struct __ATTR_ALIGN__(64) DerivAerodynAngles_r8_t 
                   {
                          double * __restrict mdPw;
                          double * __restrict mdQw;
                          double * __restrict mdRw;
                          std::size_t        mn;
 #if (USE_STRUCT_PADDING) == 1
                        PAD_TO(0,32)
#endif                 

                          DerivAerodynAngles_r8_t() = delete;

                          inline explicit DerivAerodynAngles_r8_t(const std::size_t n) noexcept(false)
                          {
                                 this->mn = n;
                                 this->allocate();
                          }

                          inline DerivAerodynAngles_r8_t(const DerivAerodynAngles_r8_t &) = delete;

                          inline DerivAerodynAngles_r8_t(DerivAerodynAngles_r8_t && rhs) noexcept(true)
                          {
                                 this->mn    = rhs.mn;
                                 this->mdPw = &rhs.mdPw[0];
                                 this->mdQw = &rhs.mdQw[0];
                                 this->mdRw = &rhs.mdRw[0];
                          }

                          inline ~DerivAerodynAngles_r8_t() noexcept(true)
                          {
                                 using namespace gms::common;
                                 gms_mm_free(this->mdPw); this->mdPw = NULL;
                                 gms_mm_free(this->mdQw); this->mdQw = NULL;
                                 gms_mm_free(this->mdRw); this->mdRw = NULL;
                          }

                          DerivAerodynAngles_r8_t & 
                          operator=(const DerivAerodynAngles_r8_t &) = delete;
                          
                          
                          inline DerivAerodynAngles_r8_t &
                          operator=(DerivAerodynAngles_r8_t && rhs) noexcept(true)
                          {
                                using namespace gms::common;
                                if(this==&rhs) return (*this);
                                gms_swap(this->mn,   rhs.mn);
                                gms_swap(this->mdPw, rhs.mdPw);
                                gms_swap(this->mdQw, rhs.mdQw);
                                gms_swap(this->mdRw, rhs.mdRw);

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
                               this->mdPw{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mdQw{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                               this->mdRw{reinterpret_cast<double * __restrict>(gms_mm_malloc(mnbytes,64ULL))};
                          }

                         bool dump_state_to_file(const char * fname) const noexcept(false)
                         {
                              assert(fname != nullptr, "Null-pointer!!");
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
                                  fprintf(fp, "mdPw[%llu,%d] = %2.9f\n",i,this->mn,this->mdPw[i]);
                              }
                              fprintf(fp,"\n\n");
                              
                              fprintf(fp,"Writing content of array mdQw[%d] to file: %s, status: started\n\n",this->mn,fname);
                              for(std::size_t i = 0ULL; i != this->mn; ++i)
                              {
                                  fprintf(fp, "mdQw[%llu,%d] = %2.9f\n",i,this->mn,this->mdQw[i]);
                              }
                              fprintf(fp,"\n\n");

                              fprintf(fp,"Writing content of array mdRw[%d] to file: %s, status: started\n\n",this->mn,fname);
                              for(std::size_t i = 0ULL; i != this->mn; ++i)
                              {
                                  fprintf(fp, "mdRw[%llu,%d] = %2.9f\n",i,this->mn,this->mdRw[i]);
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


#endif /*__GMS_DERIV_AERODYN_ANGLES_R8_HPP__*/