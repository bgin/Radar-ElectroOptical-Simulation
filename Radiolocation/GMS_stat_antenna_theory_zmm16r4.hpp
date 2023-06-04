

#ifndef __GMS_STAT_ANTENNA_THEORY_ZMM16R4_HPP__
#define __GMS_STAT_ANTENNA_THEORY_ZMM16R4_HPP__ 060620231414


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

    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_STAT_ANTENNA_THEORY_ZMM16R4_FULLVER =
      1000U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MAJOR+
      100U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MINOR+
      10U*GMS_STAT_ANTENNA_THEORY_ZMM16R4_MICRO;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_CREATION_DATE = "04-06-2023 14:14 PM +00200 (SUN 04 JUN 2023 GMT+2)";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_STAT_ANTENNA_THEORY_ZMM16R4_DESCRIPTION   = "AVX512 (single) optimized statistical antenna theory kernels.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"


















#endif /*__GMS_STAT_ANTENNA_THEORY_ZMM16R4_HPP__*/
