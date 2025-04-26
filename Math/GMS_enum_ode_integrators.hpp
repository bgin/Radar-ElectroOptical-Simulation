

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

#ifndef __GMS_ENUM_ODE_INTEGRATORS_HPP__
#define __GMS_ENUM_ODE_INTEGRATORS_HPP__ 260420250859

#include <cstdint>

namespace file_info {

     const unsigned int GMS_ENUM_ODE_INTEGRATORS_MAJOR = 1;
     const unsigned int GMS_ENUM_ODE_INTEGRATORS_MINOR = 1;
     const unsigned int GMS_ENUM_ODE_INTEGRATORS_MICRO = 0;
     const unsigned int GMS_ENUM_ODE_INTEGRATORS_FULLVER =
       1000U*GMS_ENUM_ODE_INTEGRATORS_MAJOR+100U*GMS_ENUM_ODE_INTEGRATORS_MINOR+
       10U*GMS_ENUM_ODE_INTEGRATORS_MICRO;
     const char * const GMS_ENUM_ODE_INTEGRATORS_CREATION_DATE = "26-04-2025 08:59AM  +00200 (SAT 26 APR 2025 GMT+2)";
     const char * const GMS_ENUM_ODE_INTEGRATORS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_ENUM_ODE_INTEGRATORS_SYNOPSIS      = "Enum class ODE Integrators.";

}

namespace gsm {
     
      namespace math {
            
              /*ODE Integrators scalar case*/
              enum class ODE_Integrators_Scalar : int32_t 
              {
                   embedd_fehlenberg34,
                   embedd_fehlenberg78,
                   embedd_prince_dormand45,
                   runge_kutta_f77_impl,
                   stochastic_rk
              };
            ////////////////////////////////////////////////////////////////////////////////
              /*ODE Integrators vectorized*/
              enum class ODE_Integrators_Vector : int32_t 
              {
                    
                    butcher_step_sse,
                    butcher_step_avx,
                    butcher_step_avx512,
                    euler_step_sse,
                    euler_step_avx2,
                    euler_step_avx512,
                    heuns_step_sse,
                    heuns_step_avx,
                    heuns_step_avx512,
                    nystrom_step_sse,
                    nystrom_step_avx,
                    nystrom_step_avx512,
                    rk3_step_avx,
                    rk3_step_avx512,
                    rk3v2_step_avx,
                    rk3v2_step_avx512,
                    rk4_step_avx,
                    rk4_step_avx512,
                    rk38_step_avx,
                    rk38_step_avx512,
                    rkg_step_avx,
                    rkg_step_avx512,
                    rkr2_step_avx,
                    rkr2_step_avx512,
                    rkr4_step_avx,
                    rkr4_step_avx512,
                    stochastic_rk_avx2,
                    stochastic_rk_avx512,
                    verner_step_sse,
                    verner_step_avx,
                    verner_step_avx512
              };

              constexpr bool operator==(ODE_Integrators_Scalar a,
                                        ODE_Integrators_Scalar b)
              {
                    return static_cast<ODE_Integrators_Scalar>(a)==
                           static_cast<ODE_Integrators_Scalar>(b);
              }

              constexpr bool operator!=(ODE_Integrators_Scalar a,
                                        ODE_Integrators_Scalar b)
              {
                    return static_cast<ODE_Integrators_Scalar>(a)!=
                           static_cast<ODE_Integrators_Scalar>(b);
              }

              constexpr bool operator>(ODE_Integrators_Scalar a,
                                       ODE_Integrators_Scalar b)
              {
                    return static_cast<ODE_Integrators_Scalar>(a)>
                           static_cast<ODE_Integrators_Scalar>(b);
              }

              constexpr bool operator<(ODE_Integrators_Scalar a,
                                       ODE_Integrators_Scalar b)
              {
                    return static_cast<ODE_Integrators_Scalar>(a)<
                           static_cast<ODE_Integrators_Scalar>(b);
              }

              constexpr ODE_Integrators_Scalar operator|(ODE_Integrators_Scalar a,
                                                         ODE_Integrators_Scalar b)
              {
                    return static_cast<ODE_Integrators_Scalar>(static_cast<int32_t>(a) |
                                                               static_cast<int32_t>(b));
              }

              constexpr ODE_Integrators_Scalar operator&(ODE_Integrators_Scalar a,
                                                         ODE_Integrators_Scalar b)
              {
                    return static_cast<ODE_Integrators_Scalar>(static_cast<int32_t>(a) &
                                                               static_cast<int32_t>(b));
              }

              ////////////////////////////////////////////////////////////////////////////////////////////////////

              constexpr bool operator==(ODE_Integrators_Vector a,
                                        ODE_Integrators_Vector b)
              {
                    return static_cast<ODE_Integrators_Vector>(a)==
                           static_cast<ODE_Integrators_Vector>(b);
              }

              constexpr bool operator!=(ODE_Integrators_Vector a,
                                        ODE_Integrators_Vector b)
              {
                    return static_cast<ODE_Integrators_Vector>(a)!=
                           static_cast<ODE_Integrators_Vector>(b);
              }

              constexpr bool operator>(ODE_Integrators_Vector a,
                                       ODE_Integrators_Vector b)
              {
                    return static_cast<ODE_Integrators_Vector>(a)>
                           static_cast<ODE_Integrators_Vector>(b);
              }

              constexpr bool operator<(ODE_Integrators_Vector a,
                                       ODE_Integrators_Vector b)
              {
                    return static_cast<ODE_Integrators_Vector>(a)<
                           static_cast<ODE_Integrators_Vector>(b);
              }

              constexpr ODE_Integrators_Vector operator|(ODE_Integrators_Vector a,
                                                         ODE_Integrators_Vector b)
              {
                    return static_cast<ODE_Integrators_Vector>(static_cast<int32_t>(a) |
                                                               static_cast<int32_t>(b));
              }

              constexpr ODE_Integrators_Vector operator&(ODE_Integrators_Vector a,
                                                         ODE_Integrators_Vector b)
              {
                    return static_cast<ODE_Integrators_Vector>(static_cast<int32_t>(a) &
                                                               static_cast<int32_t>(b));
              }
      }
}










#endif /*__GMS_ENUM_ODE_INTEGRATORS_HPP__*/