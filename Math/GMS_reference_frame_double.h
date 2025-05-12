
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

#ifndef __GMS_REFERENCE_FRAME_DOUBLE_H__
#define __GMS_REFERENCE_FRAME_DOUBLE_H__ 180420250713


namespace file_info {

     static const unsigned int GMS_REFERENCE_FRAME_DOUBLE_MAJOR = 1;
     static const unsigned int GMS_REFERENCE_FRAME_DOUBLE_MINOR = 1;
     static const unsigned int GMS_REFERENCE_FRAME_DOUBLE_MICRO = 0;
     static const unsigned int GMS_REFERENCE_FRAME_DOUBLE_FULLVER =
       1000U*GMS_REFERENCE_FRAME_DOUBLE_MAJOR+100U*GMS_REFERENCE_FRAME_DOUBLE_MINOR+
       10U*GMS_REFERENCE_FRAME_DOUBLE_MICRO;
     static const char  GMS_REFERENCE_FRAME_DOUBLE_CREATION_DATE[] = "18-04-2025 07:13 AM +00200 (FRI 18 APR 2025 GMT+2)";
     static const char  GMS_REFERENCE_FRAME_DOUBLE_BUILD_DATE[]    = __DATE__; 
     static const char  GMS_REFERENCE_FRAME_DOUBLE_BUILD_TIME[]    = __TIME__;
     static const char  GMS_REFERENCE_FRAME_DOUBLE_SYNOPSIS[]      = "Dynamically allocated, 64-byte aligned Reference Frame (double-precision) type.";

}

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <exception>
#include "GMS_config.h"
#include "GMS_malloc.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_REFERENCE_FRAME_DOUBLE_NT_STORES)
#define USE_GMS_REFERENCE_FRAME_DOUBLE_NT_STORES 0
#endif

namespace gms {

       namespace fdm {
            
                 /* Inertial Reference Frame single-precision */
                struct __ATTR_ALIGN__(64) ReferenceFrame_double_t 
                {       
                       double * __restrict mFI_x;    //displacement component: x
                       double * __restrict mFI_y;    //displacement component: y
                       double * __restrict mFI_z;    //displacement component: z
                       double * __restrict mdFI_x;   //velocity component: x
                       double * __restrict mdFI_y;   //velocity component: y
                       double * __restrict mdFI_z;   //velocity component: z
                       double * __restrict mddFI_x;  //acceleration component: x
                       double * __restrict mddFI_y;  //acceleration component: y
                       double * __restrict mddFI_z;  //acceleration component: z
                       std::size_t             mnx;
                       std::size_t             mny;
                       std::size_t             mnz;
                       double                   morig_x;
                       double                   morig_y;
                       double                   morig_z;
                       double                   mdt; //time increment
                       bool                    mismmap; 
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,15)
#endif                          
                       // Unit vectors
                       constexpr static double mx_hat[3] = {1.0f,0.0f,0.0f};
                       constexpr static double my_hat[3] = {0.0f,1.0f,0.0f};
                       constexpr static double mz_hat[3] = {0.0f,0.0f,1.0f};
                      
                        
              
                       inline ReferenceFrame_double_t () = delete;
                      

                       inline ReferenceFrame_double_t (const std::size_t nx,
                                                       const std::size_t ny,
                                                       const std::size_t nz,
                                                       const double       orig_x,
                                                       const double       orig_y,
                                                       const double       orig_z,
                                                       const double dt) noexcept(false)
                       {
                             assert(nx>0ULL && ny>0ULL && nz>0ULL);
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
                     
                       inline ReferenceFrame_double_t (const std::size_t nx,
                                                       const std::size_t ny,
                                                       const std::size_t nz,
                                                       const double       orig_x,
                                                       const double       orig_y,
                                                       const double       orig_z,
                                                       const double       dt,
                                                       const int32_t     prot,
                                                       const int32_t     flags,
                                                       const int32_t     fd,
                                                       const int32_t     offset,
                                                       const int32_t     fsize) noexcept(false)
                        {
                             using namespace gms::common;
                             assert(nx>0ULL && ny>0ULL && nz>0ULL);
                             this->mnx     = nx;
                             this->mny     = ny;
                             this->mnz     = nz;
                             this->morig_x = orig_x;
                             this->morig_y = orig_y;
                             this->morig_z = orig_z;
                             this->mdt     = dt;
                             switch (fsize) 
                             {
                                case 0: 
                                    this->mFI_x   = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (double* __restrict)
                                                   gms_mmap_4KiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true;
                                    break;
                                case 1:
                                    this->mFI_x   = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (double* __restrict)
                                                   gms_mmap_2MiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true; 
                                    break;
                                case 2:
                                    this->mFI_x   = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mFI_y   = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mFI_z   = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mdFI_x  = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mdFI_y  = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mdFI_z  = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mddFI_x = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mnx,prot,flags,fd,offset);
                                    this->mddFI_y = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mny,prot,flags,fd,offset);
                                    this->mddFI_z = (double* __restrict)
                                                   gms_mmap_1GiB<double>(this->mnz,prot,flags,fd,offset);
                                    this->mismmap = true;  
                                    break;
                                default : 
                                    allocate();
                                    this->mismmap = false;
                             }
                        }

                      

                        inline ReferenceFrame_double_t (ReferenceFrame_double_t && rhs) noexcept(true)
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

                        ReferenceFrame_double_t (const ReferenceFrame_double_t &) = delete;

                        inline ~ReferenceFrame_double_t () 
                        {
                            using namespace gms::common;
                            if(this->mismmap)
                            {
                                int32_t e1{},e2{},e3{},e4{},e5{},
                                                 e6{},e7{},e8{},e9{};
                                 e1 = gms_ummap<double>(this->mddFI_z,this->mnz);
                                 e2 = gms_ummap<double>(this->mddFI_y,this->mny);
                                 e3 = gms_ummap<double>(this->mddFI_x,this->mnx);
                                 e4 = gms_ummap<double>(this->mdFI_z,this->mnz);
                                 e5 = gms_ummap<double>(this->mdFI_y,this->mny);
                                 e6 = gms_ummap<double>(this->mdFI_x,this->mnx);
                                 e7 = gms_ummap<double>(this->mFI_z,this->mnz);
                                 e8 = gms_ummap<double>(this->mFI_y,this->mny);
                                 e9 = gms_ummap<double>(this->mFI_x,this->mnx);
                                 if(__builtin_expect(e1==-1,0) ||
                                    __builtin_expect(e2==-1,0) || 
                                    __builtin_expect(e3==-1,0) ||
                                    __builtin_expect(e4==-1,0) || 
                                    __builtin_expect(e5==-1,0) || 
                                    __builtin_expect(e6==-1,0) || 
                                    __builtin_expect(e7==-1,0) || 
                                    __builtin_expect(e8==-1,0) || 
                                    __builtin_expect(e9==-1,0))
                                 {
#if (FAST_TERMINATE_ON_CRITICAL_ERROR) == 1
                                __builtin_trap();
#else
                                std::terminate();
#endif                                     
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

                        ReferenceFrame_double_t & operator=(const ReferenceFrame_double_t &) = delete;

                        inline ReferenceFrame_double_t & operator=(ReferenceFrame_double_t && rhs) noexcept(true)
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

                        inline std::size_t size_nxbytes() const noexcept(true) 
                        {
                             return static_cast<std::size_t>(sizeof(double)*this->mnx);
                        }

                        inline std::size_t size_nybytes() const noexcept(true)
                        {
                             return static_cast<std::size_t>(sizeof(double)*this->mny);
                        }

                        inline std::size_t size_nzbytes() const noexcept(true)
                        {
                             return static_cast<std::size_t>(sizeof(double)*this->mnz);
                        }
                       
                        inline void info_size_and_alignment() const 
                        {
                            std::cout  << "alignof(struct ReferenceFrame_double_t) = " << alignof(ReferenceFrame_double_t) << '\n'
                                      << "sizeof(struct  ReferenceFrame_double_t)  = " << sizeof(ReferenceFrame_double_t)  << '\n'
                                      << std::hex << std::showbase            << '\n'
                                      << "&this->mFI_x  =" << (void*)&this->mFI_x    << "\n"
                                         "&this->mFI_y  =" << (void*)&this->mFI_y    << "\n"
                                         "&this->mFI_z  =" << (void*)&this->mFI_z    << "\n"
                                         "&this->mdFI_x =" << (void*)&this->mdFI_x   << "\n"
                                         "&this->mdFI_y =" << (void*)&this->mdFI_y   << "\n"
                                         "&this->mdFI_z =" << (void*)&this->mdFI_z   << "\n"
                                         "&this->mddFI_x=" << (void*)&this->mddFI_x  << "\n"
                                         "&this->mddFI_y=" << (void*)&this->mddFI_y  << "\n"
                                         "&this->mddFI_z=" << (void*)&this->mddFI_z  << "\n"
                                         "&this->mnx    =" << &this->mnx      << "\n"
                                         "&this->mny    =" << &this->mny      << "\n"
                                         "&this->mnz    =" << &this->mnz      << "\n"
                                         "&this->morig_x=" << &this->morig_x  << "\n"
                                         "&this->morig_y=" << &this->morig_y  << "\n"
                                         "&this->morig_z=" << &this->morig_z  << "\n"
                                         "&this->mdt    =" << &this->mdt      << "\n"
                                         "&this->mismmap="  << &this->mismmap << "\n";
                        }

                        private:
                        inline void allocate() noexcept(false)
                        {
                              using namespace gms::common;

                              std::size_t nxbytes{size_nxbytes()};
                              std::size_t nybytes{size_nybytes()};
                              std::size_t nzbytes{size_nzbytes()};
                              this->mFI_x         = (double* __restrict)
                                              gms_mm_malloc( nxbytes,64ULL);
                              this->mFI_y         = (double* __restrict)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mFI_z         = (double* __restrict)
                                              gms_mm_malloc( nzbytes,64ULL);
                              this->mdFI_x        = (double* __restrict)
                                              gms_mm_malloc( nxbytes,64ULL);
                              this->mdFI_y        = (double* __restrict)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mdFI_z        = (double* __restrict)
                                              gms_mm_malloc( nzbytes,64ULL);
                              this->mddFI_x       = (double* __restrict)
                                              gms_mm_malloc( nxbytes,64ULL);     
                              this->mddFI_y       = (double* __restrict)
                                              gms_mm_malloc( nybytes,64ULL);
                              this->mddFI_z       = (double* __restrict)
                                              gms_mm_malloc( nzbytes,64ULL);
                      }

                     


                };


       } // fdm

} // gms





#endif /*__GMS_REFERENCE_FRAME_DOUBLE_H__*/