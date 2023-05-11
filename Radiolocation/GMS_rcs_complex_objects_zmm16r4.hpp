

#ifndef __GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__
#define __GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__

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

namespace file_version {

    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_FULLVER =
      1000U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MAJOR+
      100U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MINOR+
      10U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MICRO;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_CREATION_DATE = "11-05-2023 10:53 PM +00200 (THR 11 MAY 2023 GMT+2)";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_DESCRIPTION   = "AVX512 optimized Complex Surfaces Radar Cross Section functionality.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"




namespace  gms {


         namespace radiolocation {
         
              
              /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Formula: 8.1-21
              */     
              
              
              
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__*/
