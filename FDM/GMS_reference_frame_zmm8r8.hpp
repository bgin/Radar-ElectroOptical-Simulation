

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

#ifndef __GMS_REFERENCE_FRAME_ZMM8R8_HPP__
#define __GMS_REFERENCE_FRAME_ZMM8R8_HPP__ 160420250121


namespace file_info {

     const unsigned int GMS_REFERENCE_FRAME_ZMM8R8_MAJOR = 1;
     const unsigned int GMS_REFERENCE_FRAME_ZMM8R8_MINOR = 1;
     const unsigned int GMS_REFERENCE_FRAME_ZMM8R8_MICRO = 0;
     const unsigned int GMS_REFERENCE_FRAME_ZMM8R8_FULLVER =
       1000U*GMS_REFERENCE_FRAME_ZMM8R8_MAJOR+100U*GMS_REFERENCE_FRAME_ZMM8R8_MINOR+
       10U*GMS_REFERENCE_FRAME_ZMM8R8_MICRO;
     const char * const GMS_REFERENCE_FRAME_ZMM8R8_CREATION_DATE = "16-04-2025 01:21 PM +00200 (WED 16 APR 2025 GMT+2)";
     const char * const GMS_REFERENCE_FRAME_ZMM8R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_REFERENCE_FRAME_ZMM8R8_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Reference Frame (AVX512-double).";

}

#include <cstdint>
#include <immintrin.h>
#include <cstdlib> //std::exit
#include "GMS_config.h"
#include "GMS_malloc.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_REFERENCE_FRAME_ZMM8R8_NT_STORES)
#define USE_GMS_REFERENCE_FRAME_ZMM8R8_NT_STORES 0
#endif

namespace gms {

       namespace fdm {

                

                /* **Reference Frame AVX-double** */
                struct __ATTR_ALIGN__(64) ReferenceFrame_zmm8r8_t 
                {       
                       __m512d * __restrict mFI_x;    //displacement component: x
                       __m512d * __restrict mFI_y;    //displacement component: y
                       __m512d * __restrict mFI_z;    //displacement component: z
                       __m512d * __restrict mdFI_x;   //velocity component: x
                       __m512d * __restrict mdFI_y;   //velocity component: y
                       __m512d * __restrict mdFI_z;   //velocity component: z
                       __m512d * __restrict mddFI_x;  //acceleration component: x
                       __m512d * __restrict mddFI_y;  //acceleration component: y
                       __m512d * __restrict mddFI_z;  //acceleration component: z
                       std::size_t              mnx;
                       std::size_t              mny;
                       std::size_t              mnz;
                       double                   morig_x;
                       double                   morig_y;
                       double                   morig_z;
                       double                   mdt; //time increment
                       bool                     mismmap; 
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,63)
#endif                      
                       // Unit vectors
                       constexpr static __m512d mx_hat[3] = {1.0,0.0,0.0};
                       constexpr static __m512d my_hat[3] = {0.0,1.0,0.0};
                       constexpr static __m512d mz_hat[3] = {0.0,0.0,1.0};
                     
                      
                 
                       inline ReferenceFrame_zmm8r8_t()
                       {
                            this->mnx      = 0ULL;
                            this->mny      = 0ULL;
                            this->mnz      = 0ULL;
                            this->mFI_x    = NULL;
                            this->mFI_y    = NULL;
                            this->mFI_z    = NULL;
                            this->mdFI_x   = NULL;
                            this->mdFI_y   = NULL;
                            this->mdFI_z   = NULL;
                            this->mddFI_x  = NULL;
                            this->mddFI_y  = NULL;
                            this->mddFI_z  = NULL;
                            this->morig_x  = 0.0;
                            this->morig_y  = 0.0;
                            this->morig_z  = 0.0;
                            this->mdt      = 0.0;
                            this->mismmap  = false;
                       }     

                       inline ReferenceFrame_zmm8r8_t(const std::size_t nx,
                                                     const std::size_t ny,
                                                     const std::size_t nz,
                                                     const double orig_x,
                                                     const double orig_y,
                                                     const double orig_z,
                                                     const double dt) noexcept(false)
                       {
                             this->mnx     = nx;
                             this->mny     = ny;
                             this->mnz     = nz;
                             allocate();
                             this->morig_x = orig_x;
                             this->morig_y = orig_y;
                             this->morig_z = orig_z;
                             this->mdt     = dt;
                             this->mismmap = false;
                       }                   
                     
                       inline ReferenceFrame_zmm8r8_t(const std::size_t nx,
                                                     const std::size_t ny,
                                                     const std::size_t nz,
                                                     const double orig_x,
                                                     const double orig_y,
                                                     const double orig_z,
                                                     const double      dt,
                                                     const int32_t prot,
                                                     const int32_t flags,
                                                     const int32_t fd,
                                                     const int32_t offset,
                                                     const int32_t fsize) noexcept(false)
                        {
                             using namespace gms::common;
                             this->mnx = nx;
                             this->mny = ny;
                             this->mnz = nz;
                             this->morig_x = orig_x;
                             this->morig_y = orig_y;
                             this->morig_z = orig_z;
                             this->mdt     = dt;
                             switch (fsize) 
                             {
                                case 0: 
                                    this->mFI_x   = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (__m512d __restrict*)
                                                   gms_mmap_4KiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true;
                                    break;
                                case 1:
                                    this->mFI_x   = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (__m512d __restrict*)
                                                   gms_mmap_2MiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true; 
                                    break;
                                case 2:
                                    this->mFI_x   = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (__m512d __restrict*)
                                                   gms_mmap_1GiB<__m512d>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true;  
                                    break;
                                default : 
                                    allocate();
                                    this->mismmap = false;
                             }
                        }

                       
                        inline ReferenceFrame_zmm8r8_t(ReferenceFrame_zmm8r8_t && rhs) noexcept(true)
                        {
                            this->mnx      = rhs.mnx;
                            this->mny      = rhs.mny;
                            this->mnz      = rhs.mnz;
                            this->mFI_x    = &rhs.mFI_x[0];
                            this->mFI_y    = &rhs.mFI_y[0];
                            this->mFI_z    = &rhs.mFI_z[0];
                            this->mdFI_x   = &rhs.mdFI_x[0];
                            this->mdFI_y   = &rhs.mdFI_y[0];
                            this->mdFI_z   = &rhs.mdFI_z[0];
                            this->mddFI_x  = &rhs.mddFI_x[0];
                            this->mddFI_y  = &rhs.mddFI_y[0];
                            this->mddFI_z  = &rhs.mddFI_z[0];
                            this->morig_x  = rhs.morig_x;
                            this->morig_y  = rhs.morig_y;
                            this->morig_z  = rhs.morig_z;
                            this->mdt      = rhs.mdt;
                            this->mismmap  = rhs.mismmap;
                        }

                        ReferenceFrame_zmm8r8_t(const ReferenceFrame_zmm8r8_t &) = delete;

                        inline ~ReferenceFrame_zmm8r8_t() noexcept(false)
                        {
                            using namespace gms::common;
                            if(this->mismmap)
                            {   // In the case of 'munmap' failure end the program execution!!
                                volatile int32_t e1,e2,e3,e4,
                                                 e5,e6,e7,e8,e9;
                                e1 =  gms_unmap<__m512d>(this->mddFI_z,this->mnz);
                                e2 =  gms_unmap<__m512d>(this->mddFI_y,this->mny);
                                e3 =  gms_unmap<__m512d>(this->mddFI_x,this->mnx);
                                e4 =  gms_unmap<__m512d>(this->mdFI_z,this->mnz);
                                e5 =  gms_unmap<__m512d>(this->mdFI_y,this->mny);
                                e6 =  gms_unmap<__m512d>(this->mdFI_x,this->mnx);
                                e7 =  gms_unmap<__m512d>(this->mFI_z,this->mnz);
                                e8 =  gms_unmap<__m512d>(this->mFI_y,this->mny);
                                e9 =  gms_unmap<__m512d>(this->mFI_x,this->mnx);
                                if(e1 || e2 || e3 || e4 || e5 || e6 || e7 || e8 || e9)
                                { 
                                    std::exit(EXIT_FAILURE);
                                }
                            }
                            else
                            {
                                 gms_mm_free(this->mddFI_z);
                                 gms_mm_free(this->mddFI_y);
                                 gms_mm_free(this->mddFI_x);
                                 gms_mm_free(this->mdFI_z);
                                 gms_mm_free(this->mdFI_y);
                                 gms_mm_free(this->mdFI_x);
                                 gms_mm_free(this->mFI_z);
                                 gms_mm_free(this->mFI_y);
                                 gms_mm_free(this->mFI_x);
                            }
                        }

                        ReferenceFrame_zmm8r8_t & operator=(const ReferenceFrame_zmm8r8_t &) = delete;

                        inline ReferenceFrame_zmm8r8_t & operator=(ReferenceFrame_zmm8r8_t && rhs) noexcept(true)
                        {
                             using namespace gms::common;
                             if(this==&rhs) return (*this);

                             gms_swap(this->mnx,rhs.mnx);
                             gms_swap(this->mny,rhs.mny);
                             gms_swap(this->mnz,rhs.mnz);
                             gms_swap(this->mFI_x,rhs.mFI_x);
                             gms_swap(this->mFI_y,rhs.mFI_y);
                             gms_swap(this->mFI_z,rhs.mFI_z);
                             gms_swap(this->mdFI_x,rhs.mdFI_x);
                             gms_swap(this->mdFI_y,rhs.mdFI_y);
                             gms_swap(this->mdFI_z,rhs.mdFI_z);
                             gms_swap(this->mddFI_x,rhs.mddFI_x);
                             gms_swap(this->mddFI_y,rhs.mddFI_y);
                             gms_swap(this->mddFI_z,rhs.mddFI_z);
                             gms_swap(this->morig_x,rhs.morig_x);
                             gms_swap(this->morig_y,rhs.morig_y);
                             gms_swap(this->morig_z,rhs.morig_z);
                             gms_swap(this->mdt,rhs.mdt);
                             gms_swap(this->mismmap,rhs.mismmap);
                             
                             return (*this);
                        }

                        inline void allocate() noexcept(false)
                        {
                              using namespace gms::common;

                              std::size_t nxbytes{size_nxbytes()};
                              std::size_t nybytes{size_nybytes()};
                              std::size_t nzbytes{size_nzbytes()};
                              this->mFI_x         = (__m512d __restrict*)
                                              gms_mm_malloc( nxbytes,64ULL);
                              this->mFI_y         = (__m512d __restrict*)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mFI_z         = (__m512d __restrict*)
                                              gms_mm_malloc( nzbytes,64ULL);
                              this->mdFI_x        = (__m512d __restrict*)
                                              gms_mm_malloc( nxbytes,64ULL);
                              this->mdFI_y        = (__m512d __restrict*)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mdFI_z        = (__m512d __restrict*)
                                              gms_mm_malloc( nzbytes,64ULL);
                              this->mddFI_x       = (__m512d __restrict*)
                                              gms_mm_malloc( nxbytes,64ULL);     
                              this->mddFI_y       = (__m512d __restrict*)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mddFI_z       = (__m512d __restrict*)
                                              gms_mm_malloc( nzbytes,64ULL);
                      }

                      inline std::size_t size_nxbytes() const noexcept(true) 
                      {
                             return static_cast<std::size_t>(sizeof(__m512d)*this->mnx);
                      }

                      inline std::size_t size_nybytes() const noexcept(true)
                      {
                             return static_cast<std::size_t>(sizeof(__m512d)*this->mny);
                      }

                      inline std::size_t size_nzbytes() const noexcept(true)
                      {
                             return static_cast<std::size_t>(sizeof(__m512d)*this->mnz);
                      }

                };

             

       } // fdm
} // gms































#endif /*__GMS_REFERENCE_FRAME_ZMM8R8_HPP__*/