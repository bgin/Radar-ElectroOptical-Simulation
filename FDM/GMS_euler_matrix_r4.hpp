
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

#ifndef __GMS_EULER_MATRIX_R4_HPP__
#define __GMS_EULER_MATRIX_R4_HPP__ 220420251232

namespace file_info {

     const unsigned int GMS_EULER_MATRIX_R4_MAJOR = 1;
     const unsigned int GMS_EULER_MATRIX_R4_MINOR = 1;
     const unsigned int GMS_EULER_MATRIX_R4_MICRO = 0;
     const unsigned int GMS_EULER_MATRIX_R4_FULLVER =
       1000U*GMS_EULER_MATRIX_R4_MAJOR+100U*GMS_EULER_MATRIX_R4_MINOR+
       10U*GMS_EULER_MATRIX_R4_MICRO;
     const char * const GMS_EULER_MATRIX_R4_CREATION_DATE = "22-04-2025 12:32  +00200 (TUE 22 APR 2025 GMT+2)";
     const char * const GMS_EULER_MATRIX_R4_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_EULER_MATRIX_R4_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Euler Matrix (generic representation) (single-precision) type.";

}

#include <cstdint>
#include <immintrin.h>
#include <cstdlib>
#include <array>
#include "GMS_config.h"
#include "GMS_malloc.h"
#if (USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"
#endif

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_EULER_MATRIX_R4_NT_STORES)
#define USE_GMS_EULER_MATRIX_R4_NT_STORES 0
#endif

namespace gms {

        namespace fdm {

                        /*
                              Generic represntation of the Euler Matrices (the 12 matrices). 
                              During this class instantiation the proper type will be named and
                              the standalone vectorized procedures will compute the results.
                        */
                       struct __ATTR_ALIGN__(64) EulerMatrix_r4_t
                       {
                              float * __restrict mr1;
                              float * __restrict mr2;
                              float * __restrict mr3;
                              float * __restrict mr4;
                              float * __restrict mr5;
                              float * __restrict mr6;
                              float * __restrict mr7;
                              float * __restrict mr8;
                              float * __restrict mr9;
                              std::size_t        mn;
                              std::array<char,6> mseq; // Axis rotation sequence, eg. "1,2,3"
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,48)
#endif                  
                              EulerMatrix_r4_t() = delete;

                              inline EulerMatrix_r4_t(const std::size_t n,
                                                      const std::array<char,6> seq) noexcept(false)
                              {
                                    this->mn   = n;
#if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_PROLOGE_BODY   
#endif                                       
                                    this->allocate();
   #if (USE_PMC_INSTRUMENTATION) == 1
                                   HW_PMC_COLLECTION_EPILOGE_BODY
#endif                                  
                                    this->mseq = seq;
                              }   

                              EulerMatrix_r4_t(const EulerMatrix_r4_t &) = delete;

                              inline EulerMatrix_r4_t(EulerMatrix_r4_t && rhs) noexcept(true)
                              {
                                     this->mn  = rhs.mn;
                                     this->mr1 = &rhs.mr1[0];
                                     this->mr2 = &rhs.mr2[0];
                                     this->mr3 = &rhs.mr3[0];
                                     this->mr4 = &rhs.mr4[0];
                                     this->mr5 = &rhs.mr5[0];
                                     this->mr6 = &rhs.mr6[0];
                                     this->mr7 = &rhs.mr7[0];
                                     this->mr8 = &rhs.mr8[0];
                                     this->mr9 = &rhs.mr9[0];
                                     this->mseq= rhs.mseq;
                              }

                              inline ~EulerMatrix_r4_t() noexcept(true)
                              {
                                  using namespace gms::common;
                                  gms_mm_free(this->mr9); this->mr9 = NULL;
                                  gms_mm_free(this->mr8); this->mr8 = NULL;
                                  gms_mm_free(this->mr7); this->mr7 = NULL;
                                  gms_mm_free(this->mr6); this->mr6 = NULL;
                                  gms_mm_free(this->mr5); this->mr5 = NULL;
                                  gms_mm_free(this->mr4); this->mr4 = NULL;
                                  gms_mm_free(this->mr3); this->mr3 = NULL;
                                  gms_mm_free(this->mr2); this->mr2 = NULL;
                                  gms_mm_free(this->mr1); this->mr1 = NULL;
                              }

                              EulerMatrix_r4_t & operator=(const EulerMatrix_r4_t &) = delete;

                              inline EulerMatrix_r4_t & operator=(EulerMatrix_r4_t && rhs) noexcept(true)
                              {
                                    using namespace gms::common;
                                    if(this==&rhs) return (*this);
                                    gms_swap(this->mn,  rhs.mn);
                                    gms_swap(this->mr1, rhs.mr1);
                                    gms_swap(this->mr2, rhs.mr2);
                                    gms_swap(this->mr3, rhs.mr3);
                                    gms_swap(this->mr4, rhs.mr4);
                                    gms_swap(this->mr5, rhs.mr5);
                                    gms_swap(this->mr6, rhs.mr6);
                                    gms_swap(this->mr7, rhs.mr7);
                                    gms_swap(this->mr8, rhs.mr8);
                                    gms_swap(this->mr9, rhs.mr9);
                                    gms_swap(this->mseq,rhs.mseq);

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
                                   this->mr1{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr2{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr3{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr4{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr5{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr6{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr7{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr8{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   this->mr9{reinterpret_cast<float * __restrict>(gmns_mm_malloc(mnbytes,64ULL))};
                                   
                              }  
                                 
                       };
        }
}







#endif /*__GMS_EULER_MATRIX_R4_HPP__*/