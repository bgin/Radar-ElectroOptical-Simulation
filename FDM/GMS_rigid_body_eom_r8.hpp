

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

#ifndef __GMS_RIGID_BODY_EOM_R8_HPP__
#define __GMS_RIGID_BODY_EOM_R8_HPP__ 180420250325


namespace file_info {

     const unsigned int GMS_RIGID_BODY_EOM_R8_MAJOR = 1;
     const unsigned int GMS_RIGID_BODY_EOM_R8_MINOR = 1;
     const unsigned int GMS_RIGID_BODY_EOM_R8_MICRO = 0;
     const unsigned int GMS_RIGID_BODY_EOM_R8_FULLVER =
       1000U*GMS_RIGID_BODY_EOM_R8_MAJOR+100U*GMS_RIGID_BODY_EOM_R8_MINOR+
       10U*GMS_RIGID_BODY_EOM_R8_MICRO;
     const char * const GMS_RIGID_BODY_EOM_R8_CREATION_DATE = "18-04-2025 03:25 PM +00200 (FRI 18 APR 2025 GMT+2)";
     const char * const GMS_RIGID_BODY_EOM_R8_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_RIGID_BODY_EOM_R8_SYNOPSIS      = "Dynamically allocated, 64-byte aligned Rigid Body Equations of Motion (double-precision) type.";

}

#include <cstdint>
#include <immintrin.h>
#include <cstdlib>
#include "GMS_config.h"
#include "GMS_malloc.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
#if !defined (USE_GMS_RIGID_BODY_EOM_R8_NT_STORES)
#define USE_GMS_RIGID_BODY_EOM_R8_NT_STORES 0
#endif


namespace  gms {
   
       namespace fdm {
                   


                   /*Rigid Body EOM*/
                   /* Bernard Etkin "Dynamics of Atmospheric Flight", page: 123 */
                   struct __ATTR_ALIGN__(64) RigidBodyEOM_r8_t 
                   {
                          // Velocity vector
                          double * __restrict mVx; 
                          double * __restrict mVy;
                          double * __restrict mVz;
                          // Acceleration vector
                          double * __restrict mAx;
                          double * __restrict mAy;
                          double * __restrict mAz;
                          // Position radius-vector
                          double * __restrict mRx;
                          double * __restrict mRy;
                          double * __restrict mRz;
                          // Rotation vector
                          double * __restrict mOx;
                          double * __restrict mOy;
                          double * __restrict mOz;   
                          std::size_t        mnx;
                          std::size_t        mny;
                          std::size_t        mnz;
                          // Inital velocity
                          double              mVx0;
                          double              mVy0;
                          double              mVz0;  
                          // Inital acceleration
                          double              mAx0;
                          double              mAy0;
                          double              mAz0;
#if (USE_STRUCT_PADDING) == 1
                      PAD_TO(0,24)
#endif
                                RigidBodyEOM_r8_t () = delete;

                        inline  RigidBodyEOM_r8_t (const std::size_t nx,
                                                  const std::size_t ny,
                                                  const std::size_t nz,
                                                  const double       Vx0,
                                                  const double       Vy0,
                                                  const double       Vz0,
                                                  const double       Ax0,
                                                  const double       Ay0,
                                                  const double       Az0) noexcept(false)
                          {
                               this->mnx  = nx;
                               this->mny  = ny;
                               this->mnz  = nz;
                               allocate();
                               this->mVx0 = Vx0;
                               this->mVy0 = Vy0;
                               this->mVz0 = Vz0;
                               this->mAx0 = Ax0;
                               this->mAy0 = Ay0;
                               this->mAz0 = Az0; 
                          }

                                RigidBodyEOM_r8_t (const RigidBodyEOM_r8_t &) = delete;
 
                        inline  RigidBodyEOM_r8_t (RigidBodyEOM_r8_t && rhs) noexcept(true)
                        {
                               this->mnx  = rhs.mnx;
                               this->mny  = rhs.mny;
                               this->mnz  = rhs.mnz;

                               this->mVx  = &rhs.mVx[0];
                               this->mVy  = &rhs.mVy[0];
                               this->mVz  = &rhs.mVz[0];

                               this->mAx  = &rhs.mAx[0];
                               this->mAy  = &rhs.mAy[0];
                               this->mAz  = &rhs.mAz[0];

                               this->mRx  = &rhs.mRx[0];
                               this->mRy  = &rhs.mRy[0];
                               this->mRz  = &rhs.mRz[0];

                               this->mOx  = &rhs.mOx[0];
                               this->mOy  = &rhs.mOy[0];
                               this->mOz  = &rhs.mOz[0];

                               this->mVx0  = rhs.mVx0;
                               this->mVy0  = rhs.mVy0;
                               this->mVz0  = rhs.mVz0;

                               this->mAx0  = rhs.mAx0;
                               this->mAy0  = rhs.mAy0;
                               this->mAz0  = rhs.mAz0;
                        }

                        inline ~RigidBodyEOM_r8_t () noexcept(true)
                        {
                               using namespace gms::common;
                               gms_mm_free(this->mOz); this->mOz = NULL;
                               gms_mm_free(this->mOy); this->mOy = NULL;
                               gms_mm_free(this->mOz); this->mOz = NULL;

                               gms_mm_free(this->mRz); this->mRz = NULL;
                               gms_mm_free(this->mRy); this->mRy = NULL;
                               gms_mm_free(this->mRx); this->mRx = NULL;

                               gms_mm_free(this->mAz); this->mAz = NULL;
                               gms_mm_free(this->mAy); this->mAy = NULL;
                               gms_mm_free(this->mAx); this->mAx = NULL;

                               gms_mm_free(this->mVz); this->mVz = NULL;
                               gms_mm_free(this->mVy); this->mVy = NULL;
                               gms_mm_free(this->mVx); this->mVx = NULL;
                        }

                               RigidBodyEOM_r8_t & operator=(const RigidBodyEOM_r8_t &) = delete;

                        inline RigidBodyEOM_r8_t & operator=(RigidBodyEOM_r8_t && rhs) noexcept(true)
                        {
                               using namespace gms::commom;
                               if(this==&rhs) return (*this);

                               gms_swap(this->mnx,rhs.mnx);
                               gms_swap(this->mny,rhs.mny);
                               gms_swap(this->mnz,rhs.mnz);

                               gms_swap(this->mVx,rhs.mVx);
                               gms_swap(this->mVy,rhs.mVy);
                               gms_swap(this->mVz,rhs.mVz);

                               gms_swap(this->mAx,rhs.mAx);
                               gms_swap(this->mAy,rhs.mAy);
                               gms_swap(this->mAz,rhs.mAz);

                               gms_swap(this->mRx,rhs.mRx);
                               gms_swap(this->mRy,rhs.mRy);
                               gms_swap(this->mRz,rhs.mRz);

                               gms_swap(this->mOx,rhs.mOx);
                               gms_swap(this->mOy,rhs.mOy);
                               gms_swap(this->mOz,rhs.mOz);

                               gms_swap(this->mVx0,rhs.mVx0);
                               gms_swap(this->mVy0,rhs.mVy0);
                               gms_swap(this->mVz0,rhs.mVz0);

                               gms_swap(this->mAx0,rhs.mAx0);
                               gms_swap(this->mAy0,rhs.mAy0);
                               gms_swap(this->mAz0,rhs.mAz0);

                               return (*this);
                        }



                        inline void allocate() noexcept(false)
                          {
                                using namespace gms::common;
                                std::size_t nxbytes{size_nxbytes()};
                                std::size_t nybytes{size_nybytes()};
                                std::size_t nzbytes{size_nzbytes()};

                                this->mVx{reinterpret_cast<double __restrict*>(gms_mm_malloc(nxbytes,64ULL))};
                                this->mVy{reinterpret_cast<double __restrict*>(gms_mm_malloc(nybytes,64ULL))};
                                this->mVz{reinterpret_cast<double __restrict*>(gms_mm_malloc(nzbytes,64ULL))};

                                this->mAx{reinterpret_cast<double __restrict*>(gms_mm_malloc(nxbytes,64ULL))};
                                this->mAy{reinterpret_cast<double __restrict*>(gms_mm_malloc(nybytes,64ULL))};
                                this->mAz{reinterpret_cast<double __restrict*>(gms_mm_malloc(nzbytes,64ULL))};

                                this->mRx{reinterpret_cast<double __restrict*>(gms_mm_malloc(nxbytes,64ULL))};
                                this->mRy{reinterpret_cast<double __restrict*>(gms_mm_malloc(nybytes,64ULL))};
                                this->mRz{reinterpret_cast<double __restrict*>(gms_mm_malloc(nzbytes,64ULL))};

                                this->mOx{reinterpret_cast<double __restrict*>(gms_mm_malloc(nxbytes,64ULL))};
                                this->mOy{reinterpret_cast<double __restrict*>(gms_mm_malloc(nybytes,64ULL))};
                                this->mOz{reinterpret_cast<double __restrict*>(gms_mm_malloc(nzbytes,64ULL))};
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

                                 

                        



                   };
       } // fdm

} // gms

































#endif /*__GMS_RIGID_BODY_EOM_R8_HPP__*/
