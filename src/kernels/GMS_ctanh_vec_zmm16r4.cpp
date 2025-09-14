




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





#include <cstdint>
#include "GMS_ctanh_vec_zmm16r4.h"
#include "GMS_csinh_vec_zmm16r4.h"
#include "GMS_ccosh_vec_zmm16r4.h"
#include "GMS_cdiv_vec_zmm16r4.h"


                   void gms::math::ctanh_zmm16r4_unroll_10x_u( const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict wrkcr,
                                                   float * __restrict wrkci,
                                                   float * __restrict wrksr,
                                                   float * __restrict wrksi,
                                                   float * __restrict ctre,
                                                   float * __restrict ctim,
                                                   const int32_t n) {

                        csinhv_zmm16r4_unroll_10x_u(xre,xim,wrksr,wrksi,n);
                        ccoshv_zmm16r4_unroll_10x_u(xre,xim,wrkcr,wrkci,n);
                        cdiv_zmm16r4_unroll_10x_u(wrksr,wrksi,wrkcr,wrkci,ctre,ctri,n);
                }   


                  
                   void gms::math::ctanh_zmm16r4_unroll_10x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) wrkcr,
                                                   float * __restrict __ATTR_ALIGN__(64) wrkci,
                                                   float * __restrict __ATTR_ALIGN__(64) wrksr,
                                                   float * __restrict __ATTR_ALIGN__(64) wrksi,
                                                   float * __restrict __ATTR_ALIGN__(64) ctre,
                                                   float * __restrict __ATTR_ALIGN__(64) ctim,
                                                   const int32_t n) {

                        csinhv_zmm16r4_unroll_10x_a(xre,xim,wrksr,wrksi,n);
                        ccoshv_zmm16r4_unroll_10x_a(xre,xim,wrkcr,wrkci,n);
                        cdiv_zmm16r4_unroll_10x_a(wrksr,wrksi,wrkcr,wrkci,ctre,ctri,n);
                }   


                   void gms::math::ctanh_zmm16r4_unroll_8x_u( const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict wrkcr,
                                                   float * __restrict wrkci,
                                                   float * __restrict wrksr,
                                                   float * __restrict wrksi,
                                                   float * __restrict ctre,
                                                   float * __restrict ctim,
                                                   const int32_t n) {

                        csinhv_zmm16r4_unroll_8x_u(xre,xim,wrksr,wrksi,n);
                        ccoshv_zmm16r4_unroll_8x_u(xre,xim,wrkcr,wrkci,n);
                        cdiv_zmm16r4_unroll_8x_u(wrksr,wrksi,wrkcr,wrkci,ctre,ctri,n);
                }   


                   
                   void gms::math::ctanh_zmm16r4_unroll_8x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) wrkcr,
                                                   float * __restrict __ATTR_ALIGN__(64) wrkci,
                                                   float * __restrict __ATTR_ALIGN__(64) wrksr,
                                                   float * __restrict __ATTR_ALIGN__(64) wrksi,
                                                   float * __restrict __ATTR_ALIGN__(64) ctre,
                                                   float * __restrict __ATTR_ALIGN__(64) ctim,
                                                   const int32_t n) {

                        csinhv_zmm16r4_unroll_8x_a(xre,xim,wrksr,wrksi,n);
                        ccoshv_zmm16r4_unroll_8x_a(xre,xim,wrkcr,wrkci,n);
                        cdiv_zmm16r4_unroll_8x_a(wrksr,wrksi,wrkcr,wrkci,ctre,ctri,n);
                }   


                  
                   void gms::math::ctanh_zmm16r4_unroll_6x_u( const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict wrkcr,
                                                   float * __restrict wrkci,
                                                   float * __restrict wrksr,
                                                   float * __restrict wrksi,
                                                   float * __restrict ctre,
                                                   float * __restrict ctim,
                                                   const int32_t n) {

                        csinhv_zmm16r4_unroll_6x_u(xre,xim,wrksr,wrksi,n);
                        ccoshv_zmm16r4_unroll_6x_u(xre,xim,wrkcr,wrkci,n);
                        cdiv_zmm16r4_unroll_6x_u(wrksr,wrksi,wrkcr,wrkci,ctre,ctri,n);
                }   


                  
                   void gms::math::ctanh_zmm16r4_unroll_6x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) wrkcr,
                                                   float * __restrict __ATTR_ALIGN__(64) wrkci,
                                                   float * __restrict __ATTR_ALIGN__(64) wrksr,
                                                   float * __restrict __ATTR_ALIGN__(64) wrksi,
                                                   float * __restrict __ATTR_ALIGN__(64) ctre,
                                                   float * __restrict __ATTR_ALIGN__(64) ctim,
                                                   const int32_t n) {

                        csinhv_zmm16r4_unroll_6x_a(xre,xim,wrksr,wrksi,n);
                        ccoshv_zmm16r4_unroll_6x_a(xre,xim,wrkcr,wrkci,n);
                        cdiv_zmm16r4_unroll_6x_a(wrksr,wrksi,wrkcr,wrkci,ctre,ctri,n);
                }   



    
