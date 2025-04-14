

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

#ifndef __GMS_INERT_REF_FRAME_XMM2R8_HPP__
#define __GMS_INERT_REF_FRAME_XMM2R8_HPP__


namespace file_info {

     const unsigned int GMS_INERT_REF_FRAME_XMM2R8_MAJOR = 1;
     const unsigned int GMS_INERT_REF_FRAME_XMM2R8_MINOR = 1;
     const unsigned int GMS_INERT_REF_FRAME_XMM2R8_MICRO = 0;
     const unsigned int GMS_INERT_REF_FRAME_XMM2R8_FULLVER =
       1000U*GMS_INERT_REF_FRAME_XMM2R8_MAJOR+100U*GMS_INERT_REF_FRAME_XMM2R8_MINOR+
       10U*GMS_INERT_REF_FRAME_XMM2R8_MICRO;
     const char * const GMS_INERT_REF_FRAME_XMM2R8_CREATION_DATE = "12-04-2025 08:25 AM +00200 (SAT 12 APR 2025 GMT+2)";
     const char * const GMS_INERT_REF_FRAME_XMM2R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_INERT_REF_FRAME_XMM2R8_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Inertial Reference Frame.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_INERT_REF_FRAME_XMM2R8_NT_STORES)
#define USE_GMS_INERT_REF_FRAME_XMM2R8_NT_STORES 0
#endif

namespace gms {

       namespace fdm {

                

                /* **Inertial Reference Frame XMM2R8-double** */
                struct __ATTR_ALIGN__(64) InertRefFrame_xmm2r8_t 
                {       
                       __m128d * __restrict mFI_x;    //displacement component: x
                       __m128d * __restrict mFI_y;    //displacement component: y
                       __m128d * __restrict mFI_z;    //displacement component: z
                       __m128d * __restrict mdFI_x;   //velocity component: x
                       __m128d * __restrict mdFI_y;   //velocity component: y
                       __m128d * __restrict mdFI_z;   //velocity component: z
                       __m128d * __restrict mddFI_x;  //acceleration component: x
                       __m128d * __restrict mddFI_y;  //acceleration component: y
                       __m128d * __restrict mddFI_z;  //acceleration component: z
                       std::size_t              mnx;
                       std::size_t              mny;
                       std::size_t              mnz;
                       double                   mdt; //time increment
                       bool                     mismmap; 
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,7)
#endif                      
                       // Unit vectors
                       constexpr static __m128d mx_hat[3] = {1.0,0.0,0.0};
                       constexpr static __m128d my_hat[3] = {0.0,1.0,0.0};
                       constexpr static __m128d mz_hat[3] = {0.0,0.0,1.0};
                       constexpr static double  morig[3]  = {0.0,0.0,0.0};
                      
                 
                       inline InertRefFrame_xmm2r8_t()
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
                            this->mdt      = 0.0;
                            this->mismmap  = false;
                       }     

                       inline InertRefFrame_xmm2r8_t(const std::size_t nx,
                                                     const std::size_t ny,
                                                     const std::size_t nz,
                                                     const double dt) noexcept(false)
                       {
                             this->mnx     = nx;
                             this->mny     = ny;
                             this->mnz     = nz;
                             allocate();
                             this->mdt     = dt;
                             this->mismmap = false;
                       }                   
                     
                       inline InertRefFrame_xmm2r8_t(const std::size_t nx,
                                                     const std::size_t ny,
                                                     const std::size_t nz,
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
                             switch (fsize) 
                             {
                                case 0: 
                                    this->mFI_x   = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (__m128d __restrict*)
                                                   gms_mmap_4KiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true;
                                    break;
                                case 1:
                                    this->mFI_x   = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (__m128d __restrict*)
                                                   gms_mmap_2MiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true; 
                                    break;
                                case 2:
                                    this->mFI_x   = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (__m128d __restrict*)
                                                   gms_mmap_1GiB<__m128d>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true;  
                                    break;
                                default : 
                                    allocate();
                                    this->mismmap = false;
                             }
                        }

                        inline InertRefFrame_xmm2r8_t(const std::size_t nx,
                                                      const std::size_t ny,
                                                      const std::size_t nz,
                                                      const __m128d * __restrict FI_x,
                                                      const __m128d * __restrict FI_y,
                                                      const __m128d * __restrict FI_z,
                                                      const __m128d * __restrict dFI_x,
                                                      const __m128d * __restrict dFI_y,
                                                      const __m128d * __restrict dFI_z,
                                                      const __m128d * __restrict ddFI_x,
                                                      const __m128d * __restrict ddFI_y,
                                                      const __m128d * __restrict ddFI_z,
                                                      const double      dt) noexcept(false)
                        {
                             using namespace gms::common;
                             this->mnx     = nx;
                             this->mny     = ny;
                             this->mnz     = nz;
                             allocate();
                             this->mdt     = dt;
                             this->mismmap = false;
#if (USE_GMS_INERT_REF_FRAME_XMM2R8_NT_STORES) == 1
                             avx512_uncached_memmove(&this->mFI_x[0], &FI_x[0], this->mnx);
                             avx512_uncached_memmove(&this->mFI_y[0], &FI_y[0], this->mny);
                             avx512_uncached_memmove(&this->mFI_z[0], &FI_z[0], this->mnz);
                             avx512_uncached_memmove(&this->mdFI_x[0], &dFI_x[0], this->mnx);
                             avx512_uncached_memmove(&this->mdFI_y[0], &dFI_y[0], this->mny);
                             avx512_uncached_memmove(&this->mdFI_z[0], &dFI_z[0], this->mnz);
                             avx512_uncached_memmove(&this->mddFI_x[0], &ddFI_x[0], this->mnx);
                             avx512_uncached_memmove(&this->mddFI_y[0], &ddFI_y[0], this->mny);
                             avx512_uncached_memmove(&this->mddFI_z[0], &ddFI_z[0], this->mnz);
#else
                             avx512_cached_memmove(&this->mFI_x[0], &FI_x[0], this->mnx);
                             avx512_cached_memmove(&this->mFI_y[0], &FI_y[0], this->mny);
                             avx512_cached_memmove(&this->mFI_z[0], &FI_z[0], this->mnz);
                             avx512_cached_memmove(&this->mdFI_x[0], &dFI_x[0], this->mnx);
                             avx512_cached_memmove(&this->mdFI_y[0], &dFI_y[0], this->mny);
                             avx512_cached_memmove(&this->mdFI_z[0], &dFI_z[0], this->mnz);
                             avx512_cached_memmove(&this->mddFI_x[0], &ddFI_x[0], this->mnx);
                             avx512_cached_memmove(&this->mddFI_y[0], &ddFI_y[0], this->mny);
                             avx512_cached_memmove(&this->mddFI_z[0], &ddFI_z[0], this->mnz);
#endif
                        }

                        inline InertRefFrame_xmm2r8_t(InertRefFrame_xmm2r8_t && rhs) noexcept(true)
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
                            this->mdt      = rhs.mdt;
                            this->mismmap  = rhs.mismmap;
                        }

                        InertRefFrame_xmm2r8_t(const InertRefFrame_xmm2r8_t &) = delete;

                        inline ~InertRefFrame_xmm2r8_t() 
                        {
                            using namespace gms::common;
                            if(this->mismmap)
                            {
                                 gms_unmap<__m128d>(this->mddFI_z,this->mnz);
                                 gms_unmap<__m128d>(this->mddFI_y,this->mny);
                                 gms_unmap<__m128d>(this->mddFI_x,this->mnx);
                                 gms_unmap<__m128d>(this->mdFI_z,this->mnz);
                                 gms_unmap<__m128d>(this->mdFI_y,this->mny);
                                 gms_unmap<__m128d>(this->mdFI_x,this->mnx);
                                 gms_unmap<__m128d>(this->mFI_z,this->mnz);
                                 gms_unmap<__m128d>(this->mFI_y,this->mny);
                                 gms_unmap<__m128d>(this->mFI_x,this->mnx);
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

                        InertRefFrame_xmm2r8_t & operator=(const InertRefFrame_xmm2r8_t &) = delete;

                        inline InertRefFrame_xmm2r8_t & operator=(InertRefFrame_xmm2r8_t && rhs) 
                        {
                             using namespace gms::common;
                             if(this==&rhs) return (*this);
                             if(this->mismmap) 
                             {
                                 gms_unmap<__m128d>(this->mddFI_z,this->mnz);
                                 gms_unmap<__m128d>(this->mddFI_y,this->mny);
                                 gms_unmap<__m128d>(this->mddFI_x,this->mnx);
                                 gms_unmap<__m128d>(this->mdFI_z,this->mnz);
                                 gms_unmap<__m128d>(this->mdFI_y,this->mny);
                                 gms_unmap<__m128d>(this->mdFI_x,this->mnx);
                                 gms_unmap<__m128d>(this->mFI_z,this->mnz);
                                 gms_unmap<__m128d>(this->mFI_y,this->mny);
                                 gms_unmap<__m128d>(this->mFI_x,this->mnx);
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
                              this->mdt      = rhs.mdt;
                              this->mismmap  = rhs.mismmap;

                              rhs.mnx        = 0ULL;
                              rhs.mny        = 0ULL;
                              rhs.mnz        = 0ULL;
                              rhs.mFI_x      = NULL;
                              rhs.mFI_y      = NULL;
                              rhs.mFI_z      = NULL;
                              rhs.mdFI_x     = NULL; 
                              rhs.mdFI_y     = NULL; 
                              rhs.mdFI_z     = NULL; 
                              rhs.mddFI_x    = NULL;
                              rhs.mddFI_y    = NULL; 
                              rhs.mddFI_z    = NULL;
                              rhs.mdt        = 0.0;
                              return (*this);
                        }

                        inline void allocate() noexcept(false)
                        {
                              using namespace gms::common;

                              std::size_t nxbytes{size_nxbytes()};
                              std::size_t nybytes{size_nybytes()};
                              std::size_t nzbytes{size_nzbytes()};
                              this->mFI_x         = (__m128d __restrict*)
                                              gms_mm_malloc( nxbytes,64ULL);
                              this->mFI_y         = (__m128d __restrict*)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mFI_z         = (__m128d __restrict*)
                                              gms_mm_malloc( nzbytes,64ULL);
                              this->mdFI_x        = (__m128d __restrict*)
                                              gms_mm_malloc( nxbytes,64ULL);
                              this->mdFI_y        = (__m128d __restrict*)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mdFI_z        = (__m128d __restrict*)
                                              gms_mm_malloc( nzbytes,64ULL);
                              this->mddFI_x       = (__m128d __restrict*)
                                              gms_mm_malloc( nxbytes,64ULL);     
                              this->mddFI_y       = (__m128d __restrict*)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mddFI_z       = (__m128d __restrict*)
                                              gms_mm_malloc( nzbytes,64ULL);
                      }

                      inline std::size_t size_nxbytes() const noexcept(true) 
                      {
                             return static_cast<std::size_t>(sizeof(__m128d)*this->mnx);
                      }

                      inline std::size_t size_nybytes() const noexcept(true)
                      {
                             return static_cast<std::size_t>(sizeof(__m128d)*this->mny);
                      }

                      inline std::size_t size_nzbytes() const noexcept(true)
                      {
                             return static_cast<std::size_t>(sizeof(__m128d)*this->mnz);
                      }

                };

             

       } // fdm
} // gms































#endif /*__GMS_INERT_REF_FRAME_XMM2R8_HPP__*/